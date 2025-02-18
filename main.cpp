#include <iostream>
#include <iomanip>
#include "MatrixOperations.h"
#include "Solver.h"
#include <cmath>
#include "PropertyPackage.h"
#include "ExpressionBuilder.h"
#include "FlashCalculation.h"
using namespace std;
vector<double> sys(const vector<double>& vars){
 
 double x1 = vars[0];
 //double x2 = vars[1];

 //double x3 = vars[2];
 //double x4 = vars[3];
 //double x5 = vars[4];

 //double x6 = vars[5];
 //double x7 = vars[6];
 //double x8 = vars[7];

 //double x9 = vars[8];
// double x10 = vars[9];
// double x4 = vars[3];
 //double x5 = vars[4];

// double x3 = vars[2];
                           
vector <double> system = {
x1*x1*x1+3*x1*x1-4*x1+0.98
};

 return system;
 
}


vector<double> systemEq(const vector<double>& vars){
 
 double x1 = vars[0];
 
vector<double> T_cr={369.8,425.2,469.7};
vector<double> P_cr={4.25,3.8,3.37};
vector<double> xmol={0.5,0.25,0.25};
vector<double> omega={0.153,0.199,0.255};
 
int nc = T_cr.size();
double T =400;
double press = 1.5;

PropertyPackage pr(nc,omega,T_cr,P_cr,xmol);
pr.nc=nc;
pr.T_cr=T_cr;
pr.P_cr=P_cr;
pr.omega=omega;
pr.x_mol=xmol;
vector<double> K = pr.calcKi(T,press);

//double result= pr.evalPengRobinsonEq(press,T,xmol,x1);
double result = rachfordRice(K,xmol,0.5);

vector <double> system = {
result
};

 return system;
 
}


double eq(double x){

vector<double> T_cr={369.8,425.2,469.7};
vector<double> P_cr={4.25,3.8,3.37};
vector<double> xmol={0.5,0.25,0.25};
vector<double> omega={0.153,0.199,0.255};
 
int nc = T_cr.size();
double T =450;
double press = 1.5;

PropertyPackage pr(nc,omega,T_cr,P_cr,xmol);
pr.nc=nc;
pr.T_cr=T_cr;
pr.P_cr=P_cr;
pr.omega=omega;
pr.x_mol=xmol;
vector<double> K = pr.calcKi(T,press);
//double result= pr.evalPengRobinsonEq(press,T,xmol,x1);
double result = rachfordRice(K,xmol,x);


return result;

}
 
double test_eq(double x){




return x*x*x*x*x+x*x*x - x*x + 2;

}


void fugacityData() {
 vector<double> A={15.726,1872.46,-25.16};//propane
 vector<double> B={15.678,2154.9,-34.42};//butane
 vector<double> C={15.883,2477.07,-39.94};//pentane

 vector<double> T_cr={369.8,425.2,469.7};
 vector<double> P_cr={4.25,3.8,3.37};
 vector<double> xmol={0.5,0.25,0.25};
 vector<double> omega={0.153,0.199,0.255};

 int nc = T_cr.size();
 double T = 400;
 double press = 1.5;

 vector<double> P= {1.5,1.6};
 int data_size = P.size();

 vector<vector<double>>Zc(3,vector<double>(data_size,0));
 vector<vector<double>>ZcDer(3,vector<double>(data_size,0));

 //vector<double> fugDer(data_size);

 //vector<double> ZcDer(data_size);
 PropertyPackage pr(nc,omega,T_cr,P_cr,xmol);
 pr.nc=nc;
 pr.T_cr=T_cr;
 pr.P_cr=P_cr;
 pr.omega=omega;
 pr.x_mol=xmol;

 for (int j=0; j<3; j++) {
  for (int i=0; i<data_size; i++) {
   Zc[j][i] = pr.analyticalPengRobinson(P[i],T,xmol)[i];

   ZcDer[j][i] = pr.analyticalDerivativeZc(P[i],T,xmol)[i];
   //fug[i] = pr.calcFi(T,P[i],xmol,Zc[i])[0];
   //fugDer[i] = pr.calcFiDer(T,P[i],xmol,ZcDer[i])[0];
  }

 }
 //printVector(ZcDer);
 //cout<<"Fugacities vector:"<<endl;
//printVector(fug);
 //cout<<"Derivative of fugacities vector:"<<endl;
 //printVector(fugDer);
 //cout<<endl;
 cout<<"-----------------Zc------------------"<<endl;
 printMatrix(Zc);
 cout<<endl;
 cout<<"-----------------ZcDer------------------"<<endl;
 printMatrix(ZcDer);
 cout<<endl;
}


void testBisection(){


bisection(test_eq,-50,50,1e-4);

}

double solVapFrac(vector<double> K,vector<double> xmol, double x){

return rachfordRice(K,xmol,x);
}

int main() {
    // Example input matrix

    std::vector<std::vector<double>> matrix = {
        {4, 7},
        {2, 6}
    };
//testBisection();
 fugacityData();
vector<double> A={15.726,1872.46,-25.16};//propane
vector<double> B={15.678,2154.9,-34.42};//butane
vector<double> C={15.883,2477.07,-39.94};//pentane

vector<double> T_cr={369.8,425.2,469.7};
vector<double> P_cr={4.25,3.8,3.37};
vector<double> xmol={0.5,0.25,0.25};
vector<double> omega={0.153,0.199,0.255};
 
int nc = T_cr.size();
double T =400;
double press = 1.5;

PropertyPackage pr(nc,omega,T_cr,P_cr,xmol);
pr.nc=nc;
pr.T_cr=T_cr;
pr.P_cr=P_cr;
pr.omega=omega;
pr.x_mol=xmol;
vector<double> Ki = pr.calcKi(T,press);


 FlashCalculation flash = FlashCalculation(pr);
 double val2 = flash.solveBubblePoint(pr,T,xmol,1e-4,100);

 /*
//vector<double> result2 = pr.calcPi_sat(400);
//vector<double> result3 = pr.calcPi(A,B,C,400);

//vector<double> sols = pr.solvePengRobinsonEq(T,press,xmol);

vector<double> fug=pr.calcFi(350,press,xmol,0.5);
cout<<" fugacity"<<endl;
printVector(fug);
cout<<endl;
FlashCalculation flash = FlashCalculation(pr);


double Pinit = flash.calcPinit(pr,xmol,T);
//cout<<"Pinit: "<<Pinit<<endl;
//double value = flash.bubblePfun(pr, T, press, xmol);
//cout<<" fun eval :"<<value<<endl;

//double value = flash.solveBubbleP(pr, T, Pinit, xmol, 1e-4,1000);

//double val2 = flash.solveBisection(pr,T,xmol,1,1e-4);

//vector<double> diff = divVec({1,2,3},{2,2,2});
//printVector(diff);



vector<double> Ki0 = pr.calcKi(T,Pinit);
vector<double> yi0 = prodVec(xmol,Ki0);
cout<<"y10 = ";
printVector(yi0);
double h =2e-15;

cout<<"P init"<<Pinit<<endl;
double Kitest = vecSum(prodVec(xmol,prodScal(vecDiff(flash.updateKi(pr,T,Pinit+h,xmol,yi0), flash.updateKi(pr,T,Pinit-h,xmol,yi0)),1/(2*h))));
double Kfun = vecSum(prodVec(xmol,flash.updateKi(pr,T,Pinit,xmol,yi0)))-1;


double diff_next = Kfun/Kitest;
cout<<" fder "<<Kitest<<"f "<<Kfun<<endl;
double P1 = Pinit-diff_next;


vector<double> y1 = prodVec(xmol,flash.updateKi(pr,T,P1,xmol,yi0));
cout<<"y1 = ";
printVector(y1);

double Kitest1 = vecSum(prodVec(xmol,prodScal(vecDiff(flash.updateKi(pr,T,P1+h,xmol,y1), flash.updateKi(pr,T,P1-h,xmol,y1)),1/(2*h))));
double Kfun1 = vecSum(prodVec(xmol,flash.updateKi(pr,T,P1,xmol,y1)))-1;


double diff_next1 = Kfun1/Kitest1;
cout<<" diff2 "<<diff_next1<<endl;
double P2 = P1-diff_next1;
cout<<"P1: "<<P2<<endl;


vector<double> y2 = prodVec(xmol,flash.updateKi(pr,T,P2,xmol,y1));
double Kitest2 = vecSum(prodVec(xmol,prodScal(vecDiff(flash.updateKi(pr,T,P2+h,xmol,y2), flash.updateKi(pr,T,P2-h,xmol,y2)),1/(2*h))));
double Kfun2 = vecSum(prodVec(xmol,flash.updateKi(pr,T,P2,xmol,y2)))-1;

double diff_next2 = Kfun2/Kitest2;
cout<<" diff2 "<<diff_next2<<endl;
double P3 = P2-diff_next2;
cout<<"P3: "<<P3<<endl;


vector<double> y3 = prodVec(xmol,flash.updateKi(pr,T,P3,xmol,y2));
double Kitest3 = vecSum(prodVec(xmol,prodScal(vecDiff(flash.updateKi(pr,T,P3+h,xmol,y3), flash.updateKi(pr,T,P3-h,xmol,y3)),1/(2*h))));
double Kfun3 = vecSum(prodVec(xmol,flash.updateKi(pr,T,P3,xmol,y3)))-1;

double diff_next3 = Kfun3/Kitest3;
cout<<" diff2 "<<diff_next2<<endl;
double P4 = P3-diff_next3;
cout<<"p4: "<<P4<<endl;

printVector(y3);
//printVector(Kitest);
//double value = flash.initial(pr, T, Pinit,xmol);
//cout<<" fun eval :"<<value<<endl;
double value = flash.bisect(pr,Pinit-1.5,Pinit+1.5,T,xmol,1e-6,2000);
//double eval_fun = flash.bubblePfun(pr,T,3.58,xmol);
//cout<<"============ "<<eval_fun<<endl;

//vector<double> testRemoveAtIndexI ={1,2,3,4,5};
//vector<double> resultTest = removeAtIndex(testRemoveAtIndexI,2);
//cout<<"------------- here  ---------"<<endl;
//printVector(Ki);
//double value = rachfordRice(Ki,xmol,0.5);
//cout<<"=========0 "<<value<<endl;
//printVector(result3);

//printVector(P_cr);
//vector<double> vec = {1,2};    
//vector<double> result = vecTimesMat(matrix,vec);
//vector<double> x0 {0.1,0.1,0.10,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
//printVector(sys(x0));
double x0 = 0.5;
//double solution = nRaphson(solVapFrac,x0,1e-4,1000);
//cout<<"solution "<<solution<<endl;
//printVector(pr.calcFi(T,press,xmol,0.5));


*/
    return 0;
}

