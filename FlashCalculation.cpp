#include<iostream>
#include<cmath>
#include "Solver.h"
#include <vector>
#include <functional>
#include "ExpressionBuilder.h"
#include "PropertyPackage.h"
#include "FlashCalculation.h"
#include "MatrixOperations.h"
//double fun(vector<double> K, vector<double> xmol, double x){
#include <algorithm>
//return rachfordRice(K, xmol,x);
//}

using namespace std;

FlashCalculation::FlashCalculation(){


}

FlashCalculation::FlashCalculation(PropertyPackage pr){

PR=pr;

}

double FlashCalculation::calcPinit(PropertyPackage pr,vector<double> xmol, double T){

vector<double> tcr = pr.T_cr;
vector<double> pcr = pr.P_cr;


int n = tcr.size();

double Pcm=0;
double Tcm=0;

for (int i=0; i<n; i++){

Pcm = Pcm + xmol[i]*pcr[i];
Tcm = Tcm + xmol[i]*tcr[i];

}

double Pinit = 1.036*Pcm*(T/Tcm);

return Pinit;

}



double FlashCalculation::solVapFrac(vector<double> K, vector<double> xmol, double x0){

//double V0;
//double error =1e-6;
//int maxIter = 10000;




auto lambda = [](vector<double> K, vector<double> xmol,double V) -> double {return rachfordRice(K,xmol,V);};
//function<double(double)> fun = lambda;
//const function<double(double)>& fun_ref = fun;
//double vaporFrac = nRaphson(fun_ref,x0,1e-5,1000);


return 0;


}

double FlashCalculation::initial(PropertyPackage pr, double T, double press, vector<double> xmol){


int n = pr.nc;
double sum = 0.0;
double Pinit = calcPinit(pr,xmol,T); 
vector<double> K0 = pr.calcKi(T,Pinit);

cout<<"Ki0-------------------"<<endl;
printVector(K0);
cout<<endl;

vector<double> y0(n);
for(int i=0; i<n; i++) {

y0[i]=xmol[i]/K0[i];

  }
cout<<"yi0-------------------"<<endl;
printVector(y0);
cout<<endl;
vector<double> ZL = pr.analyticalPengRobinson(press,T,xmol);

vector<double> ZV = pr.analyticalPengRobinson(press,T,y0);

cout<<"ZL-------------------"<<endl;
printVector(ZL);
cout<<endl;

cout<<"ZV-------------------"<<endl;
printVector(ZV);
cout<<endl;

double ZiL = find_min_max(ZL,"min",0,ZL.size());

double ZiV = find_min_max(ZV,"max",0,ZV.size());



cout<<"ZiL-------------------"<<endl;
cout<<ZiL;
cout<<endl;

cout<<"ZiV-------------------"<<endl;
cout<<ZiV;
cout<<endl;

vector<double> FiL = pr.calcFi(T,press,xmol,ZiL);

vector<double> FiV = pr.calcFi(T,press,y0,ZiV);


cout<<"FiL-------------------"<<endl;
printVector(FiL);
cout<<endl;

cout<<"FiV-------------------"<<endl;
printVector(FiV);
cout<<endl;


vector<double> K1(n);
 for (int i=0; i<n; i++){
K1[i] = FiL[i]/FiV[i];


 }

cout<<"K1-------------------"<<endl;
printVector(K1);
cout<<endl;


for (int i=0; i<n; i++){

sum = sum+xmol[i]*K1[i];

}


cout<<"FUn val-------------------"<<endl;
cout<<sum-1;
cout<<endl;

return sum-1;

}




double FlashCalculation::bubblePfun(PropertyPackage pr, double T, double press,vector<double> xmol){

int n = pr.nc;
double sum = 0.0;

double Pinit = calcPinit(pr,xmol,T); 
vector<double> K0 = pr.calcKi(T,Pinit);

;

vector<double> K(n);
vector<double> y0(n);

double vapFrac = 0.0;

vector<double> x(n);
vector<double> y(n);

vector<double> FiL(n);
vector<double> FiV(n);

double ZiL=0;
double ZiV=0;



vector<double> PiSat(n);

vector<double> sols = pr.analyticalPengRobinson(press, T, xmol); 

pr.evalPengRobinsonEq(press,T,xmol,0);
if (sols.size()>1) {

ZiL = find_min_max(sols,"min",0,sols.size()); 

ZiV = find_min_max(sols, "max",0,sols.size());  

}



else {
 ZiV = sols[0];
 ZiL = sols[0];
}

FiL = pr.calcFi(T, press, xmol, ZiL);

for (int i=0; i<n ;i++){


y0[i]=xmol[i]*K0[i];

//y0[i]=K0[i]*xmol[i];


}

FiV = pr.calcFi(T, press, y0, ZiV);


for (int i=0; i<n; i++){

K[i] = FiL[i]/FiV[i];
y[i]=xmol[i]*K[i];


//sum = sum +xmol[i]*K[i]; 
}

vector<double> new_sols = pr.analyticalPengRobinson(press, T, y);

ZiV = find_min_max(new_sols,"max",0,new_sols.size());

ZiL = find_min_max(new_sols,"min",0,new_sols.size());
FiV = pr.calcFi(T,press,y,ZiV);
FiL = pr.calcFi(T,press,xmol,ZiL);

for (int i=0; i<n; i++){
K[i] = FiL[i]/FiV[i];
sum = sum +xmol[i]*K[i]; 
} 


cout<<"---------------------Results------------------------------"<<endl;

cout<<"vapor fraction "<<endl;
printVector(y);
cout<<endl;
cout<<"Fugacities liq phase============="<<endl;
printVector(FiL);
cout<<endl;
cout<<"Fugacities vap phase============="<<endl;
printVector(FiV);
cout<<endl;
cout<<"Eq. K============================"<<endl;
printVector(K);
cout<<endl;
cout<<ZiV;
cout<<endl;
cout<<"Value in bubble point "<<sum-1<<endl;
cout<<"--------------------Results-------------------------------"<<endl;
return sum-1;
}


double FlashCalculation::solveBubbleP(PropertyPackage pr, double T, double pinit, vector<double> xmol,double error,int maxIter){

double sol = 0;
int flag=-1;
double h = 0.00001;
vector<double> P(maxIter);
vector<double> der(maxIter);
vector<double> fun(maxIter);

P[0]=pinit;

for (int i=1; i<maxIter; i++){
der[i-1] = (bubblePfun(pr,T,P[i-1]+h,xmol)-bubblePfun(pr,T,P[i-1]-h,xmol))/(2*h); 
fun[i-1] = bubblePfun(pr,T,P[i-1],xmol);

;
cout<<" derivative: "<<der[i-1]<<endl;
cout<<" act. value: "<<fun[i-1]<<endl;

P[i] =  P[i-1] - fun[i-1]/der[i-1];
cout<<"Pressure at "<<i<<"th iteration "<<P[i]<<endl;
 if (abs(P[i]-P[i-1])<= error ){

sol = P[i];
flag=1;
cout<<"Bubble pressure "<<sol<<" MPa(a) found after "<<i<<" iterations"<<endl;
}

if (flag>0) break;
}

if (flag<0) {

cout<<"Solution not found."<<endl;
}




return sol;
 
}

double FlashCalculation::bisect(PropertyPackage pr, double a, double b, double T, vector<double> xmol,double error, int maxIter){


double sol =0;

int flag =0;

int iter =0;

int n = pr.nc;
double sum = 0.0;
double Pinit = calcPinit(pr,xmol,T); 
vector<double> K0 = pr.calcKi(T,Pinit);
vector<double> K(n);



if ( bubblePfun(pr,T,a,xmol) *  bubblePfun(pr,T,b,xmol) >0 ){

cout<<"Choose a and b such f(a) * f(b) negative";
flag=-1;
}

while (iter < maxIter){

if (flag<0){
cout<<"Solution was not found"<<endl;

}

sol=(a+b)/2;

if (abs(b-a)/2<error){

cout<<"---b-a---- "<<(b-a)<<endl;
cout<<"Solution "<<sol<<" after "<<iter<<" iterations"<<endl;
return sol;

}
iter++;


if (bubblePfun(pr, T,sol,xmol) *  bubblePfun(pr, T,a,xmol)>0){


a=sol;

}
 else{

b=sol;

 } 

}
cout<<"Solution "<<sol<<" found after "<<iter<<" iterations"<<endl;
return sol;
}

vector<double> FlashCalculation::updateKi(PropertyPackage pr, double T, double press, vector<double> xmol, vector<double> y){

int n = pr.nc;

vector<double> K(n);
vector<double> obj(n);
vector<double> ZL = pr.analyticalPengRobinson(press, T, xmol);
vector<double> ZV = pr.analyticalPengRobinson(press, T, y);

double ZiL = find_min_max(ZL,"min",0,ZL.size());
double ZiV = find_min_max(ZV,"max",0,ZV.size());


vector<double> FiV = pr.calcFi(T,press,y,ZiV);
vector<double> FiL = pr.calcFi(T,press,xmol,ZiL);

for(int i=0; i<n; i++){


K[i]=FiL[i]/FiV[i];
//obj[i]=FiL[i]-FiV[i];
}

return K;

}

vector<double> FlashCalculation::objective(PropertyPackage pr, double T, double press, vector<double> xmol, vector<double> y){

 int n = pr.nc;

 vector<double> obj(n);
 vector<double> ZL = pr.analyticalPengRobinson(press, T, xmol);
 vector<double> ZV = pr.analyticalPengRobinson(press, T, y);

 double ZiL = find_min_max(ZL,"min",0,ZL.size());
 double ZiV = find_min_max(ZV,"max",0,ZV.size());


 vector<double> FiV = pr.calcFi(T,press,y,ZiV);
 vector<double> FiL = pr.calcFi(T,press,xmol,ZiL);

 for(int i=0; i<n; i++){


  //K[i]=FiL[i]-FiV[i];
  obj[i]=FiL[i]-FiV[i];
 }

 return obj;

}

vector<double> FlashCalculation::dObjective(PropertyPackage pr, double T, double press, vector<double> xmol, vector<double> y){

 int n = pr.nc;

 vector<double> K(n);
 vector<double> dObj(n);
 vector<double> ZL = pr.analyticalPengRobinson(press, T, xmol);
 vector<double> ZV = pr.analyticalPengRobinson(press, T, y);

 vector<double> ZiLDer = pr.analyticalDerivativeZc(press, T, xmol);
 vector<double> ZiVDer = pr.analyticalDerivativeZc(press, T, y);


 double ZiL = find_min_max(ZL,"min",0,ZL.size());
 int indexL = indexOf(ZL,ZiL);

 double ZiV = find_min_max(ZV,"max",0,ZV.size());
 int indexV = indexOf(ZV,ZiV);

 double ZLDer = ZiLDer[indexL];
 double ZVDer = ZiLDer[indexV];
cout<<"ZLDer-------------------------------------- : "<<ZLDer<<endl;
cout<<"ZVDer-------------------------------------- : "<<ZVDer<<endl;
 vector<double> FiLDer = pr.calcFiDer(T,press,xmol,ZiL,ZLDer);
 vector<double> FiVDer = pr.calcFiDer(T,press,xmol,ZiV, ZVDer);

 for(int i=0; i<n; i++){


  //K[i]=FiL[i]-FiV[i];
  dObj[i]=FiLDer[i]-FiVDer[i];
 }


  return dObj;
  }


double FlashCalculation::solveBubblePoint(PropertyPackage pr, double T, vector<double> xmol, double error, int maxIter){
int n= pr.nc;
int flag=-1;
double sol=0;
double h=2e-7;

vector<vector<double>> obj(maxIter, vector<double>(n, 0));
vector<vector<double>> dObj(maxIter, vector<double>(n, 0));
 vector<vector<double>> K(maxIter, vector<double>(n, 0));
 //vector<vector<double>> (maxIter, vector<double>(n, 0));
vector<vector<double>> y(maxIter, vector<double>(n, 0));
//vector<double> ZL(maxIter);
//vector<double> ZV(maxIter);

//printMatrix(K);
double Pinit = calcPinit(pr,xmol,T);
vector<double> P(maxIter);
double ZiL=0;
 double ZiV=0;
vector<double> ki0 = pr.calcKi(T,Pinit);
double sum=0;

for (int i=0; i <n; i++){

 y[0][i] = xmol[i]*ki0[i];
 K[0][i] = ki0[i];

}
P[0]=Pinit;
//printVector(y[0]);
cout<<"here Pinit "<<Pinit<<endl;
for(int i=1; i<maxIter; i++){


 vector<double> ZL = pr.analyticalPengRobinson(P[i-1], T, xmol);
 vector<double> ZV = pr.analyticalPengRobinson(P[i-1], T, y[i-1]);

 cout<<" i "<<i<<endl;
printVector(ZL);

 if (ZL.size()>1){
   ZiL = find_min_max(ZL,"min",0,ZL.size());
  ZiV = find_min_max(ZV,"max",0,ZV.size());
}
else{
  ZiL=ZL[0];
 ZiV=ZV[0];
}

 cout<<"ZiL "<<ZiL<<endl;

 vector<double> FiV = pr.calcFi(T,P[i-1],y[i-1],ZiV);
 vector<double> FiL = pr.calcFi(T,P[i-1],xmol,ZiL);


// vecSum(prodVec(xmol,prodScal(vecDiff(flash.updateKi(pr,T,P1+h,xmol,y1), flash.updateKi(pr,T,P1-h,xmol,y1)),1/(2*h))));
//double fder = vecSum(prodVec(xmol,prodScal(vecDiff(updateKi(pr,T,P[i-1]+h,xmol,y[i-1]), updateKi(pr,T,P[i-1]-h,xmol,y[i-1])),1/(2*h))));
//vecSum(prodVec(xmol,(prodScal(vecDiff(updateKi(pr,T,P[i-1]+h, xmol,y[i-1]),updateKi(pr,T,P[i-1]-h, xmol,y[i-1])) ,1/(2*h)))));
//vecSum(prodVec(xmol,prodScal(vecDiff(updateKi(pr,T,P[i-1]+h, xmol,y[i-1]),updateKi(pr,T,P[i-1]-h, xmol,y[i-1])),1/(2*h))));
//double f = vecSum(prodVec(xmol,updateKi(pr,T,P[i-1],xmol,y[i-1])))-1;

 obj[i] = objective(pr,T,P[i-1],xmol,y[i-1]);
 dObj[i] = dObjective(pr,T,P[i-1],xmol,y[i-1]);

double fder = vecSum(prodVec(xmol,dObj[i]));
double f = vecSum(prodVec(xmol,obj[i]));
//double fder=vecSum(prodScal(vecDiff(updateKi(pr,T,P[i-1]+h,xmol,y[i-1]),updateKi(pr,T,P[i-1]-h,xmol,y[i-1])),1/(2*h)));
//double f= vecSum(updateKi(pr,T,P[i-1],xmol,y[i-1]));
//cout<<"fder "<<fder<<endl;
//cout<<"f "<<f<<endl;
//-vecSum((prodVec(xmol,divVec(K[i], Kder[i])))+1;

//?addScal(vecSum(prodVec(xmol,divVec(K[i],Kder[i]))),-1)

//double eval = vecSum(prodVec(xmol,divVec(K[i],Kder[i])))-1;
y[i]=prodVec(K[i-1],xmol);
K[i] = divVec(FiL,FiV);

P[i]=P[i-1]-f/fder;
cout<<"Pressure: "<<P[i]<<endl;
printVector(y[i]);
//printVector(K[i]);
 if (abs(P[i]-P[i-1])<= error ){

sol = P[i];
flag=1;
cout<<"Bubble pressure "<<sol<<" MPa(a) found after "<<i<<" iterations"<<endl;
}

if (flag>0) break;
}

if (flag<0) {

cout<<"Solution not found."<<endl;



}
//functio


return sol;
}


double FlashCalculation::solveBisection(PropertyPackage pr, double T, vector<double> xmol, double offset, double error){

double sol=0;


int flag=-1;
int maxIter =100;

int n =pr.nc;

double Pinit = calcPinit(pr,xmol,T);
double a = Pinit-offset;
double b = Pinit+offset;

vector<vector<double>> Ka(maxIter, vector<double>(n, 0));

vector<vector<double>> Kb(maxIter, vector<double>(n, 0));

vector<vector<double>> Kc(maxIter, vector<double>(n, 0));

vector<vector<double>> K(maxIter, vector<double>(n, 0));

vector<vector<double>> ya(maxIter, vector<double>(n, 0));

vector<vector<double>> yb(maxIter, vector<double>(n, 0));

vector<vector<double>> yc(maxIter, vector<double>(n, 0));

vector<vector<double>> y(maxIter, vector<double>(n, 0));
vector<double> ki0 = pr.calcKi(T,Pinit);
 
vector<double> ki0a = pr.calcKi(T,a);

vector<double> ki0b = pr.calcKi(T,b);
vector<double> ki0c = pr.calcKi(T,(a+b)/2);

vector<double> Funi(maxIter);

for (int i=0; i <n; i++){

 ya[0][i] = xmol[i]*ki0a[i];
 Ka[0][i] = ki0a[i];

 yb[0][i] = xmol[i]*ki0b[i];
 Kb[0][i] = ki0b[i];

 Kc[0][i] = ki0c[i]; 
 yc[0][i]=xmol[i]*ki0c[i];
}
double Ptest =(a+b)/2;



for (int i=1; i<maxIter; i++){



//Ka[i]=updateKi(pr,T,3.6,xmol,{0.63,15,0.12});
Ka[i]=updateKi(pr,T,a,xmol,ya[i-1]);
Kb[i]=updateKi(pr,T,b,xmol,yb[i-1]);


ya[i] = prodVec(Ka[i-1],xmol);
yb[i] = prodVec(Kb[i-1],xmol);

double Funa = vecSum(prodVec(xmol,Ka[i-1]))-1;
double Funb = vecSum(prodVec(xmol,Kb[i-1]))-1;


cout<<"Funa "<<Funa<<endl;
cout<<"Funb "<<Funb<<endl;

    if (Funa * Funb >0){
cout<<"Choose offset such as f(Pinit-offset) * f(Pinit+offset) < 0"<<endl;

    }
   double c=a;

    while( (b-a)>=error) {
     
     c=(a+b)/2;
         cout<<"b-a "<<b-a<<endl;
        
       // cout<<"Func--------------- "<<Func<<endl;
                   // if (Func == error) return c; 
          yc[i] = prodVec(Kc[i-1],xmol);
          Kc[i]=updateKi(pr,T,c,xmol,yc[i-1]);
          double Func = vecSum(prodVec(xmol,Kc[i-1]))-1;
          cout<<"Func------------- "<<Func<<endl;          
            if (Func==error) return c;
                
            if(Func * Funa <0){
             cout<<"here-------------------------"<<endl;
          ///yc[i] = prodVec(Kc[i-1],xmol);
          //Kc[i]=updateKi(pr,T,c,xmol,yc[i-1]);
          //Func = vecSum(prodVec(xmol,Kc[i]))-1;          
                 b=c;
               }
            else {
            //yc[i] = prodVec(Ka[i-1],xmol);
          //Kc[i]=updateKi(pr,T,c,xmol,yc[i-1]);
       double Func = vecSum(prodVec(xmol,Kc[i]))-1;
            cout<<"here------2------------"<<endl;
             a=c;
         sol=c;
}
break;    
}

//cout<<"Sol------------------- "<<sol<<endl;
//cout<<"Func------------------ "<<Func<<endl;
//printMatrix(ya);
if ((b-a)<=error){
cout<<"Solution: "<<sol<<endl;
 break;
}
     }

return sol;
}
