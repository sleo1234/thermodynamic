#include<iostream>
#include<vector>
#include "MatrixOperations.h"
#include "ExpressionBuilder.h"
using namespace std;
vector<double> removeAtIndex(vector<double> vec, int index){

vector<double>  newVec;

for (int i=0; i<vec.size(); i++){


if (i != index) {

   newVec.push_back(vec[i]);
   }

 }


return newVec;
}


vector<double> generatePairs(vector<double> ki, vector<double> xmol,double val){


vector<double> newPairs;

int nc = ki.size();


vector<double> ci;

vector<double> Bij;

 for (int i=0; i<nc; i++){


 ci.push_back(xmol[i]*(ki[i]-1));

 }


for (int i=0; i<nc; i++){


Bij.push_back(1+val*(ki[i]-1));


}



for (int i=0; i<nc; i++){


newPairs.push_back(ci[i]*removeAtIndex(Bij,i)[i]);


}

return newPairs;

}

double flashEquation(vector<double> ki, vector<double> xmol,double val){

int nc=ki.size();

double eqn=0;
double nom=0;

vector<double> params=generatePairs(ki,xmol,val);

vector<double> Nomij;


for (int i=0; i<nc; i++){

eqn = eqn+params[i];

 }

for (int i=0; i< nc; i++){

Nomij.push_back(1+val*(ki[i]-1));

nom=nom+Nomij[i];
}
//double nom = vecSum(Nomij);
cout<<"denominator "<<eqn<<endl;
cout<<"nominator "<<nom<<endl;
double output = eqn/nom;

return output;
}

double rachfordRice(vector<double> K, vector<double> xmol, double V){


if (K.size() != xmol.size()){

throw invalid_argument("K and xmol must have same size");

}

int n=K.size();
double numerator =0.0;
double denominator=1.0;


 for (int i=0; i<n;i++){

 denominator *= (1.0+V*(K[i]-1.0));

}


for (int i=0; i<n; i++){

double term = xmol[i]*(K[i]-1.0);
double partial_product=1.0;




for (int j=0; j<n; j++){
   if (j!=i){
   partial_product *= (1.0+V*(K[j]-1.0));
  }


 }
numerator += term * partial_product;

}
double output = numerator/denominator;
return output;

}


double diffRachfordRice(vector<double> K, vector<double>, vector<double> xmol,double V){

double h=0.000002;


double diff= (rachfordRice(K,xmol,V+h)-rachfordRice(K,xmol,V-h))/(2*h);



return diff;
}
































