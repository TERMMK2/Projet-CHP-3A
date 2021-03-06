#ifndef LAPLACIAN_2D_H_
#define LAPLACIAN_2D_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <stdio.h>
#include <math.h>

#include <mpi.h>
#include "fonction.h"

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279
#endif

using point = std::pair<double, double>;

class Laplacian2D // pas fini de modifier
{

  public:
    enum class Source {NON, POLYNOMIAL, TRIGONOMETRIQUE, INSTATIONNAIRE};
    enum class CL {DIRICHLET, NEUMANN, NEUMANN_NON_CONSTANT};

  protected: // Les attributs de la classe
    double _x_min, _x_max, _y_min, _y_max, _h_x, _h_y, _a, _deltaT, _a_robin, _b_robin;
    int _Nx, _Ny;
    int _Nyloc; // Voir _Nxloc après.
    std::vector<std::vector<double> > _LapMatloc; // matrice creuse du laplacien
    std::vector<double> _floc; // vecteur source _f qui prend les données de _sol(i) pour calculer _sol(i+1)
    std::vector<double> _solloc; // vecteur solution U local

    CL _CL_bas, _CL_haut, _CL_gauche, _CL_droite;
    double _Val_CL_bas, _Val_CL_haut, _Val_CL_gauche, _Val_CL_droite;

    int _Me,_Np;

    Source _Source;
    int _chevauchement;

    std::string _save_all_file;
    bool _save_all_file_enabled;

    std::string _save_points_file;
    bool _save_points_file_enabled;

    std::vector<point> _saved_points;

    int _kmax; //nombre d'itération de SCHWARZ max.

  public: // Méthodes et opérateurs de la classe

    Laplacian2D();// Constructeur : Initialiser _x_min, _x_max, _y_min; _y_max; _N; _h; _LapMat; _x; _y et _sol.
    virtual ~Laplacian2D();

    // void Initialize(
    //   double x_min, double x_max, double y_min, double y_max,
    //   int Nx, int Ny, double a, double deltaT, int Me, int Np,
    //   Source source, std::string save_all_file, std::string _save_points_file,
    //   std::vector<point> saved_points);

    void InitializeCL(CL CL_bas, CL CL_haut, CL CL_gauche, CL CL_droite, double Val_CL_bas, double Val_CL_haut, double Val_CL_gauche, double Val_CL_droite);
    void Initialize(
      double x_min, double x_max,
      double y_min, double y_max,
      int Nx, int Ny, double a,
      double deltaT, int Me, int Np,
      Source Source, int chevauchement,
      std::string save_all_file, std::string _save_points_file,
      std::vector<point> saved_points,
      int kmax);

    void UpdateCL(int num_it);

    virtual void InitializeMatrix() = 0;

    void InitializeCI(double CI);

    virtual void IterativeSolver(int nb_iterations) = 0;   // Résout le système _LapMat * _sol = _f avec un solveur itératif.

    void SaveSol(const std::string& name_file); // Écrit les solutions dans le fichier "name_file".

    virtual void UpdateSecondMembre(int num_it) = 0;

    virtual std::vector<double> UpdateSchwartzCF(std::vector<double> frontiere_haut, std::vector<double> frontiere_bas) = 0;

};

class EC_ClassiqueP : public Laplacian2D
{
  //L'utilisation d'une unique classe fille vient de notre code utilisé dans notre TER sur lequel nous nous sommes basé pour écrire celui-ci. Même si elle n'est pas utile nous ne voulions pas la suprimer pour éviter de perdre du temps sur des choses peu utiles.

  public:
    void InitializeMatrix();
    void IterativeSolver(int nb_iterations);
    void UpdateSecondMembre(int num_it);
    std::vector<double> UpdateSchwartzCF(std::vector<double> frontiere_haut, std::vector<double> frontiere_bas);
};

#endif
