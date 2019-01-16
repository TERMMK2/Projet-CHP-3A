
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>

#include "Laplacian2DPara.h"
#include "getenv_var.h"
#include "Laplacian2DParaParse.h"

#define VAR_PREFIX "CHP_"

#define VAR_NX VAR_PREFIX"nx"
#define VAR_NY VAR_PREFIX"ny"

#define VAR_XMIN VAR_PREFIX"xmin"
#define VAR_XMAX VAR_PREFIX"xmax"

#define VAR_YMIN VAR_PREFIX"ymin"
#define VAR_YMAX VAR_PREFIX"ymax"

#define VAR_DT VAR_PREFIX"dt"
#define VAR_T_FINAL VAR_PREFIX"t_final"

#define VAR_CL_HAUT VAR_PREFIX"cl_haut"
#define VAR_CL_BAS VAR_PREFIX"cl_bas"
#define VAR_CL_GAUCHE VAR_PREFIX"cl_gauche"
#define VAR_CL_DROITE VAR_PREFIX"cl_droite"

#define VAR_VAL_CL_HAUT VAR_PREFIX"val_cl_haut"
#define VAR_VAL_CL_BAS VAR_PREFIX"val_cl_bas"
#define VAR_VAL_CL_GAUCHE VAR_PREFIX"val_cl_gauche"
#define VAR_VAL_CL_DROITE VAR_PREFIX"val_cl_droite"

#define VAR_SOURCE VAR_PREFIX"source"

#define VAR_SAVED_POINT VAR_PREFIX"saved_points"
#define VAR_SAVE_ALL_FILE VAR_PREFIX"save_all_file"
#define VAR_SAVE_POINTS_FILE VAR_PREFIX"save_points_file"

#define VAR_CHEVAUCHEMENT VAR_PREFIX"chevauchement"

using namespace std;


int main(int argc, char *argv[])
{
  //MPI_Status status;

  MPI_Init(NULL, NULL);

  int Me, Np;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);

  const int Nx = getenv_var<int>(VAR_NX, 100);
  const int Ny = getenv_var<int>(VAR_NY, 100);

  // Truc machin biduletest

  const double xmin = getenv_var<double>(VAR_XMIN, 0.0);
  const double xmax = getenv_var<double>(VAR_XMAX, 1.);
  const double ymin = getenv_var<double>(VAR_YMIN, 0.0);
  const double ymax = getenv_var<double>(VAR_YMAX, 1.);

  const int chevauchement = getenv_var<int>(VAR_CHEVAUCHEMENT, 10);

//---
  double a = 1.; //Mettre 1. si on fait les cas tests de l'énoncé et 1./(1500.*1000.) si on veut comparer avec notre TER.
  double deltaT = getenv_var<double>(VAR_DT, 0.1);
  double tfinal = getenv_var<double>(VAR_T_FINAL, 10.0);

  const Laplacian2D::CL CL_bas = getenv_var<Laplacian2D::CL>(VAR_CL_BAS, Laplacian2D::CL::DIRICHLET);
  const Laplacian2D::CL CL_haut = getenv_var<Laplacian2D::CL>(VAR_CL_HAUT, Laplacian2D::CL::DIRICHLET);
  const Laplacian2D::CL CL_gauche = getenv_var<Laplacian2D::CL>(VAR_CL_GAUCHE, Laplacian2D::CL::DIRICHLET); //Peut valoir "Neumann_non_constant", cette condition est celle proposée par Mme Baranger dans notre TER, un flux au bord affine par morceaux
  const Laplacian2D::CL CL_droite = getenv_var<Laplacian2D::CL>(VAR_CL_DROITE, Laplacian2D::CL::DIRICHLET);

  //Il sagit de la valeur du Flux si CL_bas == "Neumann", ou de la Température si CL_bas == "Dirichlet"

  const double Val_CL_bas = getenv_var<double>(VAR_VAL_CL_BAS, 0.0);
  const double Val_CL_haut = getenv_var<double>(VAR_VAL_CL_HAUT, 0.0);
  const double Val_CL_gauche = getenv_var<double>(VAR_VAL_CL_GAUCHE, 3000.0);
  const double Val_CL_droite = getenv_var<double>(VAR_VAL_CL_DROITE, 0.0);

  int nb_iterations = int(ceil(tfinal / deltaT));
  //Peut prendre comme valeur "non", "polynomial", "trigonometrique" ou "instationnaire".
  Laplacian2D::Source Source = getenv_var<Laplacian2D::Source>(VAR_SOURCE, Laplacian2D::Source::POLYNOMIAL);

  double CI = 290.;

  string save_all_file = getenv_var<std::string>(VAR_SAVE_ALL_FILE, "SCHWARZ"); //Mettre "non" si on ne souhaite pas enregistrer la solution globale au cours du temps sous une forme lisible par paraview

  string save_points_file = getenv_var<std::string>(VAR_SAVE_POINTS_FILE, "points_SCHWARZ"); //Mettre non si on ne veut pas sauvegarder la température au cours du temps en des points paritculiers
  vector<point> saved_points = getenv_var<vector<point> >(VAR_SAVED_POINT, {{0.0, 0.0025}, {0.002, 0.0025}, {0.004, 0.0025}});

  //  Check param

  for (const auto& p : saved_points) {
    if (p.first > xmax || p.first < xmin || p.second > ymax || p.second < ymin)
      throw std::runtime_error("saved point must be in [xmin, xmax] * [ymin, ymax]");
  }

  if (xmin >= xmax)
    throw std::runtime_error("xmin must be lower than xmax");

  if (ymin >= ymax)
    throw std::runtime_error("ymin must be lower than ymax");

  if (deltaT < 0.0)
    throw std::runtime_error("deltaT must be greater than 0.0");

  if (chevauchement + 1 > Ny/Np)
    throw std::runtime_error("Chevauchement times the numbers of processors must be lesser than the numbers of lines");


  EC_ClassiqueP Lap;

  //-------------------------------------------------------------------------
  if(Me ==0)
    cout<<"Début de l'initialisation"<<endl;
  Lap.Initialize(xmin,xmax,ymin,ymax,Nx,Ny,a,deltaT, Me, Np, Source, chevauchement, save_all_file, save_points_file, saved_points);
  Lap.InitializeCI(CI);
  Lap.InitializeCL(CL_bas, CL_haut, CL_gauche, CL_droite, Val_CL_bas, Val_CL_haut, Val_CL_gauche, Val_CL_droite);
  Lap.InitializeMatrix();
  if(Me == 0)
    cout<<"Fin de l'initialisation"<<endl;


  //-------------------------------------------------------------------------

  auto start = chrono::high_resolution_clock::now();
  Lap.IterativeSolver(nb_iterations);
  auto finish = chrono::high_resolution_clock::now();

  double t = chrono::duration_cast<chrono::microseconds>(finish - start).count();

  if (Me == 0)
  {
    cout << "Le prog a mis " << t * 0.000001 << " secondes a s'effectuer" << endl;
  }

  MPI_Finalize();

  return 0;
}
