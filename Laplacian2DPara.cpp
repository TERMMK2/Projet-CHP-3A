
#include "Laplacian2DPara.h"

using namespace std;

//Constructeur :
Laplacian2D::Laplacian2D()
{
}
//Destructeur :
Laplacian2D::~Laplacian2D()
{
}

void Laplacian2D::Initialize(
    double x_min, double x_max,
    double y_min, double y_max,
    int Nx, int Ny, double a,
    double deltaT, int Me, int Np,
    Source Source, int chevauchement,
    string save_all_file, string save_points_file,
    vector<point> saved_points)
{
  // // On  initialise les constantes connues de tous les processeurs.

  _x_min = x_min;
  _y_min = y_min;
  _x_max = x_max;
  _y_max = y_max;
  _Nx = Nx;
  _Ny = Ny;
  _a = a;
  _deltaT = deltaT;
  _h_y = (y_max - y_min) / (Ny + 1.);
  _h_x = (x_max - x_min) / (Nx + 1.);
  _Me = Me;
  _Np = Np;
  _Source = Source;
  _chevauchement = chevauchement;

  _save_all_file = save_all_file;
  _save_points_file = save_points_file;

  _save_all_file_enabled = (save_all_file != "non");
  _save_points_file_enabled = (save_points_file != "non");

  _saved_points = move(saved_points);

  if (_save_all_file_enabled) //On supprime l'ancien fichier qui contient les solutions au cours du temps et on en crée un nouveau
  {
    system(("rm -Rf " + _save_all_file).c_str());
    system(("mkdir -p ./" + _save_all_file).c_str());
  }

  int i1, iN;
  charge(_Ny, _Np, _Me, i1, iN);
  if (_Me != _Np-1){
    iN += _chevauchement;
    _Nyloc = iN - i1 + 1;
  }
  else {
    _Nyloc = iN - i1 + 1;
  }

  //à modifier potentiellement pour savesol et save_points
}

void Laplacian2D::InitializeCI(double CI)
{
  // // On initialise le vecteur solution ici.

  _solloc.resize(_Nyloc * _Nx);
  for (int i = 0; i < _Nyloc * _Nx; i++)
  {
    _solloc[i] = CI;
  }

  //On met _floc de la bonne taille ici  :
  _floc.resize(_Nyloc);
}

void Laplacian2D::InitializeCL(CL CL_bas, CL CL_haut, CL CL_gauche, CL CL_droite, double Val_CL_bas, double Val_CL_haut, double Val_CL_gauche, double Val_CL_droite)
{
  // // On initialise les condition limites ici. La configuration quelque peu redondante avec Laplacian2D::Initialize vient d'une ancienne version du code de notre TER dans laquelle on initialisait toutes ces valeurs dans le main et où on ne souhaitait pas avoir trop d'arguments dans la méthode Initialize.

  _CL_bas = CL_bas;
  _CL_haut = CL_haut;
  _CL_gauche = CL_gauche;
  _CL_droite = CL_droite;
  _Val_CL_bas = Val_CL_bas;
  _Val_CL_haut = Val_CL_haut;
  _Val_CL_gauche = Val_CL_gauche;
  _Val_CL_droite = Val_CL_droite;
}

void Laplacian2D::UpdateCL(int num_it)
{
  //Cette méthode de classe nous permet de faire une condition au limites de Neumann variable au cours du temps. C'est avec cela que nous avons fait notre test pour vérifier notre code.

  double t = num_it * _deltaT;
  if (t <= 50.)
  {
    _Val_CL_gauche = 10000. * t;
  }
  else
  {
    _Val_CL_gauche = -9000. * (t - 50.) + 500000.;
  }
}

void EC_ClassiqueP::InitializeMatrix()
{
  // On initialise la matrice pentadiagonale ici.

  _LapMatloc.resize(5);
  int Nloc = _Nx * _Nyloc;

  double alpha = 1 + 2 * _a * _deltaT / (_h_x * _h_x) + 2 * _a * _deltaT / (_h_y * _h_y);
  double beta = -_a * _deltaT / (_h_x * _h_x);
  double gamma = -_a * _deltaT / (_h_y * _h_y);

  for (int i = 0; i < 5; i++)
  {
    _LapMatloc[i].resize(Nloc);
  }

  for (int i = 0; i < Nloc; i++)
  {
    _LapMatloc[0][i] = gamma;
    _LapMatloc[1][i] = beta;
    _LapMatloc[2][i] = alpha;
    _LapMatloc[3][i] = beta;
    _LapMatloc[4][i] = gamma;

    if (i % _Nx == 0)
      _LapMatloc[1][i] = 0.;

    if (i % _Nx == _Nx - 1)
      _LapMatloc[3][i] = 0.;

    if (i < _Nx)
      _LapMatloc[0][i] = 0.;

    if (i > (_Nyloc - 1) * _Nx - 1)
      _LapMatloc[4][i] = 0.;
  }

  if (_CL_gauche == CL::NEUMANN or _CL_gauche == CL::NEUMANN_NON_CONSTANT)
  {
    for (int i = 0; i < _Nyloc; i++)
    {
      _LapMatloc[2][_Nx * i] += beta; //Bord gauche
    }
  }

  if (_CL_droite == CL::NEUMANN)
  {
    for (int i = 0; i < _Nyloc; i++)
    {
      _LapMatloc[2][_Nx * (i + 1) - 1] += beta; //Bord droit
    }
  }

  if (_CL_haut == CL::NEUMANN)
  {
    for (int i = 0; i < _Nx; i++)
    {
      _LapMatloc[2][i] += gamma; //Bord haut
    }
  }

  if (_CL_bas == CL::NEUMANN)
  {
    for (int i = 0; i < _Nx; i++)
    {
      _LapMatloc[2][(_Nyloc - 1) * _Nx + i] += gamma; //Bord bas
    }
  }

}

void EC_ClassiqueP::IterativeSolver(int nb_iterations)
{
  //  Cette méthode est au coeur de la résolution du problème, elle permet d'effectuer la marche en temps
  MPI_Status status;

  int Nloc = _Nx * _Nyloc;

  //initialisation schwartz
  std::vector<double> solloc_km;
  std::vector<double> solloc_k(_solloc);
  std::vector<double> floc_k;
  solloc_km.resize(Nloc);
  floc_k.resize(Nloc);

  //  Prealocate record data
  if (_Me == 0) {
    _record_data.resize(nb_iterations + 1);
    for (auto& data : _record_data)
      data.reserve(_Nx * _Ny);
  }

  //  Prealocate frontière share buffer
  std::vector<double> frontiere_haut(_Nx);
  std::vector<double> frontiere_bas(_Nx);

  double norme;
  int kmax = Nloc + 100; //Pour l'algo du GC : Pour une matrice de taille
  //n le GC met max n étapes en théorie, comme on veut être sûr qu'il converge on prend une petite marge

  //Sauvegarde d'un point ou plusieurs points particulier au cours du temps:------------------------------------------------------------------------------------

  ofstream flux_pts;
  vector<double> _sol;

  if (_save_points_file_enabled) // à refaire (besoin de savoir dans quel proc est le point pour l'enregistrer)
  {
    //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
    flux_pts.open(_save_points_file + ".txt", ios::out);
    _sol.resize(_Nx * _Ny);
  }
  //-------------------------------------------------------------------------------------------------------------------------------------------------------------

  for (int iter = 0; iter <= nb_iterations; iter++)
  {

    //-------------PARTIE CONCERNANT LA SAUVEGARDE DE LA SOLUTION------------------------------------------
    if (_save_all_file_enabled) {
      EC_ClassiqueP::SaveSol(iter);
    }

    int i1, iN;
    charge(_Ny, _Np, _Me, i1, iN);

    if ((_save_points_file_enabled) and (_Me == 0))
    {

      vector<double> sol;
      sol.resize(_Nx * _Ny);
      for (int i = 0; i < (iN+1) * _Nx; i++)
      {
        sol[i] = _solloc[i];
      }

      for (int he = 1; he < _Np; he++)
      {
        int he_i1, he_iN;
        charge(_Ny, _Np, he, he_i1, he_iN);

        vector<double> sol_temp;
        sol_temp.resize((he_iN - he_i1 + 1) * _Nx);

        MPI_Recv(&sol_temp[0], (he_iN - he_i1 + 1) * _Nx, MPI_DOUBLE, he, 100 * he, MPI_COMM_WORLD, &status);

        for (int i = he_i1 * _Nx; i < (he_iN+1) * _Nx; i++)
        {
          sol[i] = sol_temp[i - he_i1 * _Nx];
        }
      }

      flux_pts << iter * _deltaT << " ";
      for (const auto &p : _saved_points)
      {
        const int pos = floor(p.first / _h_x) + _Nx * floor(p.second / _h_y);
        flux_pts << sol[pos] << " ";
      }
      flux_pts << endl;
    }
    if ((_save_points_file_enabled) and (_Me != 0))
    {
      MPI_Send(&_solloc[0], (iN - i1 + 1) * _Nx, MPI_DOUBLE, 0, 100 * _Me, MPI_COMM_WORLD);
    }
    //----------------FIN DE LA SAUVEGARDE DE LA SOLUTION-----------------------------

    UpdateSecondMembre(iter);

    //------------------------------------------------------------------------

    double condition_arret = 1.; //juste pour rentrer une première fois dans la boucle. Changeable en sa valeur réelle mais plus long pour pas grand chose?
    const double epsilon = 0.000001; //valeur arbitraire pour l'instant

    //-------------------debut boucle schwartz--------------------------------

    int k = 0;
    while (condition_arret > epsilon) //condition d'arrêt de la boucle de Schwartz
    {
      k++;

      if (_Np > 1)
      {
        if (_Me == 0)
        {
          MPI_Send(&solloc_k[(_Nyloc - 1 - _chevauchement) * _Nx], _Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
          MPI_Recv(&frontiere_bas[0], _Nx, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, &status);
        }

        for (int he = 1; he < _Np - 1; he++)
        {
          if (_Me == he)
          {
            MPI_Send(&solloc_k[_chevauchement*_Nx], _Nx, MPI_DOUBLE, he - 1, 1000 * _Me, MPI_COMM_WORLD);
            MPI_Recv(&frontiere_haut[0], _Nx, MPI_DOUBLE, he - 1, 100 * (he - 1), MPI_COMM_WORLD, &status);

            MPI_Send(&solloc_k[(_Nyloc - 1 - _chevauchement) * _Nx], _Nx, MPI_DOUBLE, he + 1, 100 * _Me, MPI_COMM_WORLD);
            MPI_Recv(&frontiere_bas[0], _Nx, MPI_DOUBLE, he + 1, 1000 * (he + 1), MPI_COMM_WORLD, &status);
          }
        }

        if (_Me == _Np - 1)
        {
          MPI_Send(&solloc_k[_chevauchement*_Nx], _Nx, MPI_DOUBLE, _Np - 2, 1000 * _Me, MPI_COMM_WORLD);
          MPI_Recv(&frontiere_haut[0], _Nx, MPI_DOUBLE, _Np - 2, 100 * (_Np - 2), MPI_COMM_WORLD, &status);
        }
      }
      else {
        std::fill(frontiere_bas.begin(), frontiere_bas.end(), 0.0);
        std::fill(frontiere_haut.begin(), frontiere_haut.end(), 0.0);
      }

      UpdateSchwartzCF(floc_k, frontiere_haut, frontiere_bas);
      solloc_km = solloc_k;

      CG(solloc_k, _LapMatloc, floc_k, solloc_km, 0.000001, kmax, _Nx, _Nyloc);

      double condition_arret_loc;
    
      if (_Me == 0)
      {
        double s = 0.0;

        for (int j = 0; j < _Nx; j++)
        {
          const int index = _Nx * (_Nyloc - 1) + j;
          const double component = solloc_k[index] - solloc_km[index];
          s += component * component;
        }

        condition_arret_loc = sqrt(s);
      }
      else if (_Me == _Np - 1)
      {
        double s = 0.0;

        for (int j = 0; j < _Nx; j++)
        {
          const double component = solloc_k[j] - solloc_km[j];
          s += component * component;
        }

        condition_arret_loc = sqrt(s);
      }
      else
      {
        double s1 = 0.0;
        double s2 = 0.0;

        for (int j = 0; j < _Nx; j++)
        {
          const double component = solloc_k[j] - solloc_km[j];
          s1 += component * component;
        }

        for (int j = 0; j < _Nx; j++)
        {
          const int index = _Nx * (_Nyloc - 1) + j;
          const double component = solloc_k[index] - solloc_km[index];
          s2 += component * component;
        }

        condition_arret_loc = 
          sqrt(s1) + sqrt(s2);
      }

      if (_Np > 1)
        MPI_Allreduce(&condition_arret_loc, &condition_arret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      else
        condition_arret = 0.0;
    }

    //fin de boucle schwartz

    _solloc = solloc_k;

    if (_Me == 0) //Barre de chargement
    {
      int i_barre;
      const int p = floor((((double)iter) / ((double)nb_iterations)) * 100);
      printf("[");
      for (i_barre = 0; i_barre <= p; i_barre += 2)
        printf("*");
      for (; i_barre <= 100; i_barre += 2)
        printf("-");
      printf("] %3d %%", p);
      for (i_barre = 0; i_barre < 59; ++i_barre)
        printf("%c", 8);

      fflush(stdout);
    }

  }

  // if ((_save_points_file != "non") and (_Me == 0))
  //   {
  //     flux_pts->close();
  //   }
  // delete flux_pts;

  if (_Me == 0) //Barre de chargement
  {
    printf("\n");
  }
}

void Laplacian2D::SaveSol(const int iter)
{
  // À refaire complètement -> à voir

  MPI_Status status;

  int i1, iN;
  charge(_Ny, _Np, _Me, i1, iN);
  if (_Me == 0)
  {
    vector<double> sol(_Nx * _Ny);

    for (int i = 0; i < (iN+1) * _Nx; i++)
      sol[i] = _solloc[i];

    for (int he = 1; he < _Np; he++)
    {
      int he_i1, he_iN;
      charge(_Ny, _Np, he, he_i1, he_iN);
      vector<double> sol_temp;
      sol_temp.resize((he_iN - he_i1 + 1) * _Nx);

      MPI_Recv(&sol_temp[0], (he_iN - he_i1 + 1) * _Nx, MPI_DOUBLE, he, 100 * he, MPI_COMM_WORLD, &status);

      for (int i = he_i1 * _Nx; i < (he_iN+1) * _Nx; i++)
      {
        sol[i] = sol_temp[i - he_i1 * _Nx];
      }
    }

    auto &data = _record_data[iter];

    for (int i = _Ny - 1; i >= 0; i--)
    {
      for (int j = 0; j < _Nx; j++)
      {
        data.push_back(sol[j + i * _Nx]);
      }
    }
  }
  else
  {
    MPI_Send(&_solloc[0], (iN - i1 + 1) * _Nx, MPI_DOUBLE, 0, 100 * _Me, MPI_COMM_WORLD);
  }

}

void Laplacian2D::write_record_data() const
{
  if (_Me != 0)
    return;

  const int size = _record_data.size();
  
  for (int iter = 0; iter < size; ++iter) {
    ofstream mon_flux;
    mon_flux.open(_save_all_file + "/sol_it_" + to_string(iter) + ".vtk", ios::out);
    mon_flux << "# vtk DataFile Version 3.0" << endl
            << "cell" << endl
            << "ASCII" << endl
            << "DATASET STRUCTURED_POINTS" << endl
            << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << endl
            << "ORIGIN 0 0 0" << endl
            << "SPACING " + to_string((_x_max - _x_min) / _Nx) + " " + to_string((_y_max - _y_min) / _Ny) + " 1" << endl
            << "POINT_DATA " << _Nx * _Ny << endl
            << "SCALARS sample_scalars double" << endl
            << "LOOKUP_TABLE default" << endl;

    const auto& data = _record_data[iter];
    int k = 0;

    for (int i = _Ny - 1; i >= 0; i--)
    {
      for (int j = 0; j < _Nx; j++)
      {
        // double x = j*_h_y;
        // double y = i*_h_x;
        // mon_flux << sol[j + i * _Nx] - ( x*x - x )*( y*y - y ) << " ";
        // mon_flux << sol[j + i*_Nx] - ( sin(x) + cos(y) ) << " ";
        mon_flux << data[k++] << " ";
      }
      mon_flux << endl;
    }
    mon_flux.close();
  }
}

void EC_ClassiqueP::UpdateSecondMembre(int num_it)
{
  // Cette méthode nous permet de mettre à jours le terme source à chaque itération
  //pour prendre en compte les effets des conditions limites et le terme source

  //------------------
  //CL
  //------------------
  double gamma = -_a * _deltaT / (_h_y * _h_y);
  double beta = -_a * _deltaT / (_h_x * _h_x);

  int i1, iN;
  charge(_Ny, _Np, _Me, i1, iN);
  if (_Me != _Np-1){
    iN += _chevauchement;
  }

  _floc = _solloc;

  if ((_CL_haut == CL::DIRICHLET) and (_Me == 0)) //Condition de température en haut
  {
    for (int j = 0; j < _Nx; j++)
    {
      _floc[j] = _solloc[j] - gamma * _Val_CL_haut;
    }
  }
  if ((_CL_bas == CL::DIRICHLET) and (_Me == _Np - 1)) //Condition de température en bas
  {
    for (int j = 0; j < _Nx; j++)
    {
      _floc[_Nx * (_Nyloc - 1) + j] = _solloc[_Nx * (_Nyloc - 1) + j] - gamma * _Val_CL_bas;
    }
  }
  if (_CL_gauche == CL::DIRICHLET) //Condition de température à gauche
  {
    for (int i = 0; i < _Nyloc; i++)
    {
      _floc[i * _Nx] = _solloc[i * _Nx] - beta * _Val_CL_gauche;
    }
  }
  if (_CL_droite == CL::DIRICHLET) //Condition de température à droite
  {
    for (int i = 0; i < _Nyloc; i++)
    {
      _floc[(i + 1) * _Nx - 1] = _solloc[(i + 1) * _Nx - 1] - beta * _Val_CL_droite; //i*_Nx+(_Nx-1)
    }
  }
  if ((_CL_haut == CL::NEUMANN) and (_Me == 0)) //Condition de flux en haut
  {
    for (int j = 0; j < _Nx; j++)
    {
      _floc[j] = _solloc[j] - gamma * _Val_CL_haut * _h_y;
    }
  }
  if ((_CL_bas == CL::NEUMANN) and (_Me == _Np - 1)) //Condition de flux en bas
  {
    for (int j = 0; j < _Nx; j++)
    {
      _floc[_Nx * (_Nyloc - 1) + j] = _solloc[_Nx * (_Nyloc - 1) + j] - gamma * _Val_CL_bas * _h_y;
    }
  }
  if (_CL_gauche == CL::NEUMANN) //Condition de flux à gauche
  {
    for (int i = 0; i < _Nyloc; i++)
    {
      _floc[i * _Nx] = _solloc[i * _Nx] - beta * _Val_CL_gauche * _h_x;
    }
  }
  if (_CL_gauche == CL::NEUMANN_NON_CONSTANT) //Condition de flux à gauche
  {
    Laplacian2D::UpdateCL(num_it);
    for (int i = 0; i < _Nyloc; i++)
    {
      _floc[i * _Nx] = _solloc[i * _Nx] - beta * _Val_CL_gauche * _h_x;
    }
  }
  if (_CL_droite == CL::NEUMANN) //Condition de flux à droite
  {
    for (int i = 0; i < _Nyloc; i++)
    {
      _floc[(i + 1) * _Nx - 1] = _solloc[(i + 1) * _Nx - 1] - beta * _Val_CL_droite * _h_x; //i*_Nx+(_Nx-1)
    }
  }
  if (_Source == Source::TRIGONOMETRIQUE)
  {
    if (_Me == 0)
    {
      for (int j = 0; j < _Nx; j++) // En haut
      {
        double x = (j % _Nx) * _h_x;
        double y = 0.;
        _floc[j] = _solloc[j] - gamma * (sin(x) + cos(y));
      }
    }
    if (_Me == _Np - 1)
    {
      for (int j = 0; j < _Nx; j++) // En bas
      {
        double x = ((_Nx * (_Nyloc - 1) + j) % _Nx) * _h_x;
        double y = _y_max;
        _floc[_Nx * (_Nyloc - 1) + j] = _solloc[_Nx * (_Nyloc - 1) + j] - gamma * (sin(x) + cos(y));
      }
    }
    for (int i = 0; i < _Nyloc; i++) // A gauche
    {
      double x = 0.;
      double y = ((i * _Nx) / _Nx + i1) * _h_y;
      _floc[i * _Nx] = _solloc[i * _Nx] - beta * (sin(x) + cos(y));
    }
    for (int i = 0; i < _Nyloc; i++) // A droite
    {
      double x = _x_max;
      double y = (((i + 1) * _Nx - 1) / _Nx + i1) * _h_y;
      _floc[(i + 1) * _Nx - 1] = _solloc[(i + 1) * _Nx - 1] - beta * (sin(x) + cos(y));
    }
  }

  //------------------
  //Source Term
  //------------------
  for (int k = 0; k < _Nyloc * _Nx; k++)
  {
    double x = (k % _Nx) * _h_x;
    double y = (k / _Nx + i1) * _h_y;

    if (_Source == Source::POLYNOMIAL)
    {
      _floc[k] = _floc[k] + 2 * _deltaT * (y - y * y + x - x * x);
    }
    if (_Source == Source::TRIGONOMETRIQUE)
    {
      _floc[k] = _floc[k] + _deltaT * (sin(x) + cos(y));
    }
    if (_Source == Source::INSTATIONNAIRE)
    {
      _floc[k] = _floc[k] + _deltaT * exp(-(x / 2.) * (x / 2.)) * exp(-(y / 2) * (y / 2)) * cos(M_PI * num_it * _deltaT / 2.);
    }
  }
}

void EC_ClassiqueP::UpdateSchwartzCF(std::vector<double>& floc_k, const std::vector<double> &frontiere_haut, const std::vector<double> &frontiere_bas)
{
  //Cette methode permet de mettre à jour floc pour prendre en compte l'évolution des conditions à la frontiere entre les procs dans la boucle de schwartz
  double gamma = -_a * _deltaT / (_h_y * _h_y);
  floc_k = _floc;

  if ((_Me > 0) and (_Me < _Np - 1))
  {
    for (int j = 0; j < _Nx; j++)
    {
      floc_k[_Nx * (_Nyloc - 1) + j] = _floc[_Nx * (_Nyloc - 1) + j] - gamma * frontiere_bas[j];
      floc_k[j] = _floc[j] - gamma * frontiere_haut[j];
    }
  }

  if (_Me == 0)
  {
    for (int j = 0; j < _Nx; j++)
    {
      floc_k[_Nx * (_Nyloc - 1) + j] = _floc[_Nx * (_Nyloc - 1) + j] - gamma * frontiere_bas[j];
    }
  }

  if (_Me == _Np - 1)
  {
    for (int j = 0; j < _Nx; j++)
    {
      floc_k[j] = _floc[j] - gamma * frontiere_haut[j];
    }
  }

}

