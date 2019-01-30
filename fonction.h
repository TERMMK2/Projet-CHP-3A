#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>

void charge (int N, int Np, int Me, int& i1, int& iN);

void prodMVC(std::vector<double>& output, const std::vector<std::vector<double>> &Aloc, const std::vector<double> &xloc, int nx, int ny);

double dot(const std::vector<double>& u, const std::vector<double>& v);

void CG(std::vector<double>& output, const std::vector<std::vector<double>> &A, const std::vector<double> &b, const std::vector<double> &x0, double err, int kmax, int nx, int ny);
