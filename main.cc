
#include <iostream>
#include <cmath>
#include <vector>

#include "fonction.h"
#include "getenv_var.h"

#include "Laplacian2DParaParse.h"


int main(int argc, char * argv[])
{
  using namespace std;

  using pair = std::pair<float, float>;

  const char *str = "(1, 2), (4, 5.8), (8.0, 2.4)";
  char *end;

  std::vector<pair> vec{
    str_parser<std::vector<pair> >::parse(str, &end)
  };

  if (end == str)
    std::cout << "MAtch error" << std::endl;

  for (const auto& v : vec) {
    std::cout << "(" << v.first << " ," << v.second << ")" << std::endl;
  }
  

  /*
  double err = 1.e-3;
  int kmax = 100000;

  const int nx=50;
  const int ny=50;
  int N = nx*ny;
  //double alpha, beta, gamma, dx, dy;
  std::vector<std::vector<double> > A;
  std::vector<double> x,x0,b;
  x.resize(N);
  x0.resize(N);
  b.resize(N);
  


  for (int i=0 ; i<N; i++)
  {
    x[i] = 1.;
    x0[i] = 0.;
  }



  Diag_init(nx,ny,A); //Initialisation des diagonales de la matrice A
  
  b = prodMVC(A,x,nx,ny); //Produit Matrice-Vecteur creux


  x = CG(A, b, x0, err, kmax, nx, ny);

  std::cout << "x =" << std::endl;
  for (int i=0 ; i<N ; i++)
  {
    std::cout << x[i] <<std::endl;
  }


  return 0;
*/
}










 // for(int i=0; i<N; i++)
  //   {
  //     A[i].resize(N);
  //     x[i] = 1.;
  //     if (i<nx){ 
  // 	if (i==0){
  // 	  A[i][0] = D3[i];
  // 	  A[i][1] = D4[i];
  // 	  A[i][nx] = D5[i];
  // 	}
  // 	else{
  // 	  A[i][i-1] = D2[i];
  // 	  A[i][i] = D3[i];
  // 	  A[i][i+1] = D4[i];
  // 	  A[i][i+nx] = D5[i];
  // 	}
  //     }

  //     else if (i>((ny-1)*nx-1)){ 
  // 	if (i==N-1){
  // 	  A[i][N-2] = D2[i];
  // 	  A[i][N-1] = D3[i];
  // 	  A[i][N-1-nx] = D1[i];
  // 	}
  // 	else{
  // 	  A[i][i-1] = D2[i];
  // 	  A[i][i] = D3[i];
  // 	  A[i][i+1] = D4[i];
  // 	  A[i][i-nx] = D1[i];
  // 	}
  //     }

  //     else{
  // 	A[i][i-1] = D2[i];
  // 	A[i][i] = D3[i];
  // 	A[i][i+1] = D4[i];
  // 	A[i][i-nx] = D1[i];
  // 	A[i][i+nx] = D5[i];
  //     }

  //   }
  
  //prodMV(argc, argv, A, x); //Produit Matrice-Vecteur
