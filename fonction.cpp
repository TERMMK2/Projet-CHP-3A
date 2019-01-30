#include "fonction.h"

void charge(int N, int Np, int Me, int &i1, int &iN)
{
  const int q = N / Np;
  const int r = N - q * Np;

  if (Me < r)
  {
    i1 = (q + 1) * Me;
    iN = (q + 1) * (Me + 1) - 1;
  }
  else
  {
    i1 = q * Me + r;
    iN = q * (Me + 1) + r - 1;
  }
}

void prodMVC(std::vector<double>& output, const std::vector<std::vector<double>>& A, const std::vector<double>& x, int Nx, int Ny)
{
  int N = Nx * Ny;

  if (output.size() != N)
    output.resize(N);

  for (int i = 0; i < N; i++)
  {
    output[i] =
        A[0][i] * x[i - Nx]
      + A[1][i] * x[i - 1]
      + A[2][i] * x[i]
      + A[3][i] * x[i + 1]
      + A[4][i] * x[i + Nx];
  }

}

double dot(const std::vector<double>& u, const std::vector<double>& v)
{
  //Calcul le produit scalaire entre 2 vecteurs.

  const int N = u.size();
  double y = 0.0;

  for (int i = 0; i < N; i++)
    y += u[i] * v[i];

  return y;
}

void CG(std::vector<double> &output, const std::vector<std::vector<double>> &A, const std::vector<double> &b, const std::vector<double> &x0, double err, int kmax, int nx, int ny)
{
  // Algorithme du gradient conjugué parallèle qui prend en argument uniquement des vecteurs locaux et renvoie un vecteur local.

  int k = 0;

  double norm_r, nr_carre, nr2_carre;
  std::vector<double> w, r, r_1, p, d;//, x;

  int N = nx * ny;

  w.resize(N);
  r.resize(N);
  p.resize(N);
  d.resize(N);

  output = x0;//x = x0;

  prodMVC(w, A, output, nx, ny);

  for (int i = 0; i < N; i++)
  {
    r[i] = b[i] - w[i];
  }
  p = r;

  nr_carre = dot(r, r);
  norm_r = sqrt(nr_carre); //On stocke ces deux variables puisqu'on s'en sert toutes les deux plusieurs fois dans la suite

  while ((norm_r > err) and (k < kmax))
  {
    prodMVC(d, A, p, nx, ny);

    double alpha = nr_carre / dot(p, d);

    for (int i = 0; i < N; i++)
    {
      output[i] += alpha * p[i];
      r[i] -= alpha * d[i]; //rk devient r(k+1)
    }

    nr2_carre = dot(r, r); //norme de r(k+1)

    double beta = nr2_carre / nr_carre;

    for (int i = 0; i < N; i++)
    {
      p[i] = r[i] + beta * p[i];
    }

    nr_carre = nr2_carre;
    norm_r = sqrt(nr_carre);

    k++;
  }

  //return x;
}
