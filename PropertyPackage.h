#ifndef PROPERTYPACKAGE_H
#define PROPERTYPACKAGE_H
#include <vector>
#include <cmath>
#include <functional>
#include "MatrixOperations.h"

using namespace std;


class PropertyPackage {

double R=8.31;

public:
 int nc;
 vector<double> omega;
 vector<double> T_cr;
 vector<double> P_cr;
 vector<double> x_mol;

public:
 PropertyPackage();
 PropertyPackage(int nc, vector<double> omega,vector<double> T_cr, vector<double> P_cr, vector<double> x_mol);
 
 vector<double> calcKi (double temp, double press); 
 vector<double> calcPi_sat(double temp); 
 vector<double> calcPi( vector<double> A,  vector<double> B, vector<double> C, double T); 
 vector<double> a_M(double temp);
 vector<double> lambda_vec();
 vector<double> alfa_m(double temp);
 vector<vector<double>> getAij(double temp, double p);
 vector<vector<double>> getAi(vector<vector<double>> mat, vector<double> xmol);
 double alfam(double temp, vector<double> xmol);
 vector<double> b_M();
 double covolParam(vector<double> xmol);
 double attractParam(double temp, vector<double> xmol);
 vector<double> analyticalPengRobinson(double press, double temp, vector<double> xmol);
 double evalPengRobinsonEq(double press, double temp, vector<double> xmol, double Zc0);
 void calcPengRobinsonParam(double T, double press, vector<double> xmol); 
 vector<double> solvePengRobinsonEq(double T, double press, vector<double> xmol);
 vector<double> calcFi(double T, double press, vector<double> xmol,double Zalfa);
 double mixRule(vector<double> xmol, vector<double> property);
 vector<double> analyticalDerivativeZc(double press, double temp, vector<double> xmol);
 vector<vector<double>> getAijDer(double temp, double p);
 vector<double> calcFiDer(double T, double press, vector<double> xmol, double Zalfa, double ZalfaDer);
};



#endif
