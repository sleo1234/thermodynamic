#ifndef EXPRESSIONBUILDER_H
#define EXPRESSIONBUILDER_H

#include<vector>

using namespace std;

vector<double> removeAtIndex(vector<double>vec, int i);
vector<double> generatePairs(vector<double> ki, vector<double> xmol, double val);
double flashEquation(vector<double> ki, vector<double> xmol, double val);
double rachfordRice(vector<double> K, vector<double> xmol, double V);
double diffRachfordRice(vector<double> K, vector<double> xmol, double V);
double dRachfordRice(vector<double> K, vector<double> xmol, double V);
#endif
