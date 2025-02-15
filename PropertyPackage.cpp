#include <iostream>
#include "PropertyPackage.h"
#include <vector>
#include "Solver.h"
#include <string.h>
double R=8.31;



PropertyPackage::PropertyPackage(){


}

PropertyPackage::PropertyPackage(int Nc, vector<double> acc, vector<double> Tcr, vector<double> Pcr, vector<double> xmol){
nc=Nc;
omega=acc;
T_cr=Tcr;
P_cr=Pcr;
x_mol=xmol;


}

vector<double> PropertyPackage::calcKi (double temp, double press) {
        vector<double> K_i(nc);
          
                
        for (int i=0; i < nc; i++) {

            K_i[i] = (P_cr[i]/press)*exp(5.37*(1+omega[i])*(1-T_cr[i]/temp));
              
        }
   return K_i;
    }


vector<double> PropertyPackage::calcPi_sat (double temp){

vector<double> P_i(nc);

for (int i=0; i < nc; i++) {

  P_i[i] = exp(log(P_cr[i])+(7.0/3.0)*log(10.0)*(1+omega[i])*(1-T_cr[i]/temp));
 


}


return P_i;


}

//Antoine equation for vapor pressure 
vector<double> PropertyPackage::calcPi(vector<double> A,vector<double> B,vector<double> C, double temp){
        vector<double> Pi (nc);

        for (int i=0; i < nc; i++) {
            Pi[i] = exp(A[i]-B[i]/(temp+C[i]))/7600.0;
        }
        return Pi;
}

vector<double> PropertyPackage::lambda_vec(){



vector<double> lamb_vec(nc);

vector<double>acc = omega;

        for (int i=0; i < nc; i++) {
            if (acc[i] < 0.5215){

              lamb_vec[i] = 0.37464 + 1.5423 * acc[i] - 0.26992 * acc[i] * acc[i];
            }

            else {
                lamb_vec[i] = 0.3796 + 1.485 * acc[i] - 0.1644 * acc[i] * acc[i] + 0.01666 * acc[i] * acc[i] * acc[i];
            }

        }
        return lamb_vec;



}



vector<double> PropertyPackage::alfa_m(double temp){
vector<double> alfa(nc);
vector<double> lambda = lambda_vec();

  for (int i=0; i < nc; i++) {
         alfa[i] = pow(1+ lambda[i] * (1 - sqrt(temp/T_cr[i])),2);
         
        }
return alfa;
}

vector<double> PropertyPackage::a_M(double temp){

   vector<double> a(nc);
        vector<double> alfa = alfa_m (temp);
        for (int i=0; i < nc; i++) {
            a[i] = 0.45724 * alfa[i] * R * R * (T_cr[i] * T_cr[i]) / P_cr[i];

        }

        return a;
}



vector<vector<double>>  PropertyPackage::getAij(double temp, double p){
double k = p/(temp*temp*R*R);
vector<double> a(nc);
vector<double> alfa = alfa_m(temp);
vector<vector<double>> Aij(nc,vector<double>(nc));

  for (int i=0; i<nc;i++){
    for (int j=0; j<nc;j++){
    
    Aij[i][j] = k * sqrt(a_M(temp)[i]*a_M(temp)[j]);
 
        }
      }

return Aij;
}


vector<vector<double>>  PropertyPackage::getAijDer(double temp, double p){
    double k = 1/(temp*temp*R*R);
    vector<double> a(nc);
    vector<double> alfa = alfa_m(temp);
    vector<vector<double>> Aij(nc,vector<double>(nc));

    for (int i=0; i<nc;i++){
        for (int j=0; j<nc;j++){

            Aij[i][j] = k * sqrt(a_M(temp)[i]*a_M(temp)[j]);

        }
    }

    return Aij;
}





vector<vector<double>> PropertyPackage::getAi(vector<vector<double>> mat,vector<double> xmol){

int n= xmol.size();

vector<vector<double>>result(n,vector<double>(n)); 

for (int i=0; i<n; i++){


result[i] =  prodVec(xmol,mat[i]);

}

return result;

}


double PropertyPackage::alfam(double temp, vector<double> xmol){

vector<double> alfa = alfa_m(temp);

double sum=0;

for(int i=0; i<nc;i++){

   for(int j=0; j<nc;j++){

if (xmol[i]==0) {
      alfa[i]=0;
      alfa[j]=0;
        sum = sum+xmol[j]*alfa[i]*alfa[j];

   }



   }

}


return sum;
}




vector<double> PropertyPackage::b_M(){

 vector<double> b(nc);

   for (int i=0; i<nc; i++){
        b[i] = 0.077796*R*(T_cr[i]/P_cr[i]);


      }




return b;
}

double PropertyPackage::covolParam(vector<double> xmol){


double B = 0;

vector<double> bm = b_M();
double sum =0;

for (int i=0; i<nc; i++){

sum = sum + bm[i] * xmol[i];


}

B=sum;

return B;
}

double PropertyPackage::attractParam(double temp, vector<double> xmol){
double A=0.0;
vector<double> am = a_M(temp);

double sum =0;
 for (int i=0; i<nc; i++){
     for (int j=0; j<nc; j++){

   sum = sum + xmol[i]*xmol[j]*sqrt(am[i]*am[j]);


    }
  }
    
A = sum;

return A;


}

        vector<double> PropertyPackage::analyticalPengRobinson(double press, double temp, vector<double> xmol)  {
        vector<double> sols;
        double pi = 4*atan(1);

        double A = (attractParam (temp,xmol) * press) / (R * R * temp * temp);
        double B = (covolParam(xmol) * press) /(R * temp);

        double C2 = B-1.0;
        double C1 = (A-3*B*B-2*B);
        double C0 =(B*B*B+B*B-A*B);
        double Q1 = C2*C1/6.0-C0/2.0-(C2*C2*C2)/27.0;
        double P1 = C2*C2/9.0-C1/3.0;
        double D = Q1*Q1-P1*P1*P1;
  
        double teta = 0.0;

        if ( D >= 0.0){
            double sign1 = (Q1+sqrt(D))/(abs(Q1+sqrt(D)));
            double sign2 = (Q1-sqrt(D))/(abs(Q1-sqrt(D)));
            double sol1 =sign1*pow(abs((Q1+sqrt(D))),1.0/3.0)+sign2*pow(abs((Q1-sqrt(D))),1.0/3.0)-C2/3.0;
            sols.push_back(sol1);
            
        }

        else {

            
            double t1 = (Q1*Q1)/(P1*P1*P1);
            double t2 = (sqrt(1-t1))/((sqrt(t1)*(Q1/abs(Q1))));
            if (atan(t2) <0){
                teta = atan(t2)+pi;
            }
            else {
                teta = atan(t2);
            }
            double sol1 = 2*sqrt(P1)*cos(teta/3.0)-C2/3.0;
            double sol2 = 2*sqrt(P1)*cos((teta+2*pi)/3.0)-C2/3.0;
            double sol3 = 2*sqrt(P1)*cos((teta+4*pi)/3.0)-C2/3.0;
            sols.push_back(sol1);
            sols.push_back(sol2);
            sols.push_back(sol3);
        }
        return sols;
    }
 
double PropertyPackage::evalPengRobinsonEq(double press, double temp, vector<double> xmol,double Zc0) {


 double A = (attractParam (temp,xmol) * press) / (R * R * temp * temp);
 double B = (covolParam(xmol) * press) /(R * temp);

 double coeff1 = 1.0-B;
 double coeff2 = (A-2*B-3*B*B);
 double coeff3 =(A*B-B*B-B*B*B);
 string exp = "Zc^3-"+to_string(coeff1)+"*Zc^2+"+to_string(coeff2)+"*Zc-"+to_string(coeff3);
 
double val = Zc0*Zc0*Zc0-coeff1*Zc0*Zc0+coeff2*Zc0-coeff3;


return val;

}


vector<double> PropertyPackage::analyticalDerivativeZc(double press, double temp, vector<double> xmol){

    double pi = 4*atan(1);
    vector<double> sols;
    double A = (attractParam (temp,xmol) * press) / (R * R * temp * temp);
    double B = (covolParam(xmol) * press) /(R * temp);

    double C2 = B-1.0;
    double C1 = (A-3*B*B-2*B);
    double C0 =(B*B*B+B*B-A*B);
    double Q1 = C2*C1/6.0-C0/2.0-(C2*C2*C2)/27.0;
    double P1 = C2*C2/9.0-C1/3.0;
    double D = Q1*Q1-P1*P1*P1;





    double Ader = (attractParam (temp,xmol)) / (R * R * temp * temp);
    double Bder = (covolParam(xmol)) /(R * temp);

    double C2der = Bder;
    double C1der = (Ader - 6*Bder*B-2*Bder);

    double C0der = 3*B*B*Bder+2*Bder*B-(Ader*B+A*Bder);
    double Q1der = (C2der*C1+C1der*C2)/6-C0der/2-3*C2*C2*C2der/27;
    double P1der = C2*C2*C2der/9-C1der/3;
    double Dder = Q1*Q1-P1*P1*P1; 2*Q1*Q1der-3*P1*P1*P1der;

    double teta = 0.0;
    double tetaDer = 0.0;


    if ( D >= 0.0){
        double sign1 = (Q1+sqrt(D))/(abs(Q1+sqrt(D)));
        double sign2 = (Q1-sqrt(D))/(abs(Q1-sqrt(D)));
        double sol1 =sign1*pow(abs((Q1+sqrt(D))),1.0/3.0)+sign2*pow(abs((Q1-sqrt(D))),1.0/3.0)-C2/3.0;

        double sign1der,sign2der=0;
        double QplusD = abs(Q1+sqrt(D));
        double QplusDder = abs(Q1der +D/(2*sqrt(D)));
        double QminusD = abs(Q1+sqrt(D));
        double QminusDder = abs(Q1der - D/(2*sqrt(D)));

        double first_term = sign1*pow(abs((Q1+sqrt(D))),1.0/3.0);
        double first_termDer = sign1der*pow(abs((Q1+sqrt(D))),1.0/3.0)+sign1*(1/3)*pow(abs((Q1+sqrt(D))),-2.0/3.0)*QplusDder;

        double second_term = sign1*pow(abs((Q1-sqrt(D))),1.0/3.0);
        double second_termDer = sign2der*pow(abs((Q1+sqrt(D))),1.0/3.0)+sign2*(1/3)*pow(abs((Q1-sqrt(D))),-2.0/3.0)*QminusDder;

        double sol1der= first_termDer+second_termDer;
        sols.push_back(sol1der);

    }

    else {


        double t1 = (Q1*Q1)/(P1*P1*P1);

        double first_term=(Q1*Q1);
        double second_term=(P1*P1*P1);

        double Q1sqDer = 2*Q1*Q1der;
        double P13Der = 3*P1*P1*P1der;

        double t1der = (1/P13Der*P13Der)*(Q1sqDer*second_term-first_term*P13Der);


        double t2 = (sqrt(1-t1))/((sqrt(t1)*(Q1/abs(Q1))));//-----(f/(g*Q1/(abs(Q1))
                                                          //-----------==
        /*

        (sqrt(1-t1))/((sqrt(t1)*(Q1/abs(Q1)))) =  f
                                                 ------
                                                  g * Q1
                                                     ----
                                                      abs(Q1)
            f = sqrt(1-t1), g = sqrt(t1)                */
        double f = sqrt(1-t1);
        double g = sqrt(t1);
        double gder = t1der/(2*sqrt(t1));
        double fder = -1*t1der /(2*sqrt(1-t1));

        double gTimesQ1absQ1 = g*Q1/abs(Q1);
        double gTimesQ1absQ1der = gder*Q1/abs(Q1);
        double t2der = (fder * gTimesQ1absQ1 - f * gTimesQ1absQ1der)/(gTimesQ1absQ1 * gTimesQ1absQ1);


        if (atan(t2) <0){
            teta = atan(t2)+pi;
            tetaDer = t2der/(atan(t2)*atan(t2)+1);

        }
        else {
            teta = atan(t2);
            tetaDer = t2der/(atan(t2)*atan(t2)+1);
        }
        //double sol1 = 2*sqrt(P1)*cos(teta/3.0)-C2/3.0;

        double f1 = 2*sqrt(P1);
        double f2 = cos(teta/3.0);
        double f1der = P1der / sqrt(P1);
        double f2der = -sin(teta/3.0) * (tetaDer / 3);

        double sol1Der = f1der * f2 + f1 * f2der-(1/3)*C2der;
        /*

         */
        double sol2 = 2*sqrt(P1)*cos((teta+2*pi)/3.0)-C2/3.0;

        double f22 = cos((teta+2*pi)/3.0);
        double f22der = -sin((teta+2*pi)/3.0) * (tetaDer / 3);

        double sol2Der = f1der * f22 + f22der * f1 - (1/3)*C2der;


        double sol3 = 2*sqrt(P1)*cos((teta+4*pi)/3.0)-C2/3.0;

        double f33 = cos((teta+4*pi)/3.0);
        double f33der = -sin((teta+4*pi)/3.0) * (tetaDer / 3);

        double sol3Der = f1der * f33 + f33der * f1 - (1/3)*C2der;
        sols.push_back(sol1Der);
        sols.push_back(sol2);
        sols.push_back(sol3);
        }


return sols;
}





void PropertyPackage::calcPengRobinsonParam(double T, double press, vector<double> xmol){


b_M();
alfa_m(T);
a_M(T);
lambda_vec();
attractParam(T,xmol);
covolParam(xmol);


}

vector<double> PropertyPackage::solvePengRobinsonEq(double T, double press, vector<double> xmol) {

calcPengRobinsonParam(T,press,xmol);


vector<double> sols =analyticalPengRobinson(press, T, xmol);

int noOfSOlutions = sols.size();



return sols;

}


//check for mutliple real solutions
//highest compresibility factor - vapor phas, lowest - liquid phase;
 

// calculate fugacity coefficients
vector<double> PropertyPackage::calcFi(double T, double press, vector<double> xmol, double Zalfa){


double Vm =(Zalfa*R*T)/press;//cm3/mol
double bm = covolParam(xmol);
double aa = attractParam(T, xmol);
vector<double> bi = b_M();
double rad2 = sqrt(2);
vector<double> ffi(nc);
vector<double> fug(nc);
double Bm = bm*press/(R*T);
double Am = aa*press/(R*R*T*T);
vector<double> Bi(nc);
vector<double> vecSum = getVecSum(xmol,getAij(T,press));


for (int i=0; i<nc;i++){

Bi[i]=press*bi[i]/(R*T);

 ffi[i] =(Bi[i]/Bm)*(Zalfa-1.0)-log(Zalfa-Bm)-
                    (Am/(2*rad2*Bm))*( (2.0/Am)*vecSum[i] - Bi[i]/Bm)*log((Zalfa+(1+rad2)*Bm)/(Zalfa+(1-rad2)*Bm));

fug[i] = exp(ffi[i]);

cout<<"Bi------------ "<<Bi[i]<<endl;

cout<<"Bm------------ "<<Bm<<endl;

cout<<"Am------------ "<<Am<<endl;
  }

return fug;

}


vector<double> PropertyPackage::calcFiDer(double T, double press, vector<double> xmol, double Zalfa){

    double ZalfaDer = analyticalDerivativeZc(T,press,xmol)[0];

    double Vm =(Zalfa*R*T)/press;//cm3/mol
    double VmDer = ZalfaDer*R*T/press-(1/press*press)*Zalfa;

    double bm = covolParam(xmol);

    double aa = attractParam(T, xmol);

    vector<double> bi = b_M();

    double rad2 = sqrt(2);

    vector<double> ffi(nc);
    vector<double> fug(nc);

    double Bm = bm*press/(R*T);
    double Am = aa*press/(R*R*T*T);

    double BmDer = bm/(R*T);
    double AmDer = aa/(R*R*T*T);

    vector<double> Bi(nc);
    vector<double> BiDer(nc);
    vector<double> Bi_BmDer(nc);
    vector<double> vecSum = getVecSum(xmol,getAij(T,press));
    vector<double> vecSumDer = getVecSum(xmol,getAijDer(T,press));

    //derivative of log(Zalfa-Bm)
    double logZalfa_Bm = (ZalfaDer - BmDer)/(Zalfa-Bm);
    //derivative of (Am/(2*rad2*Bm)) -----> AmDer / (2*rad2*Bm) - (2*rad2*BmDer) * Am
    //                                                         --------------------------
    //                                                         ((2*rad2*Bm)*(2*rad2*Bm))

    double Am_rad2_Bm = AmDer/(2*rad2*Bm) - ((2*rad2*BmDer) * Am)/((2*rad2*Bm)*(2*rad2*Bm));
    //derivative (2.0/Am)*vecSum[i]
    vector<double> Am_VecSum(nc);

    double logZalfa_rad2 = (ZalfaDer+rad2*BmDer)/(Zalfa+(1+rad2*Bm)) - (ZalfaDer-rad2*BmDer)/(Zalfa+(1-rad2*Bm));
    //2*vecSum[i] '   2  * vecSum[i] + vecSumDer[i] * 2
    //-----------  = (---) '                         ------
    //    Am          Am                                Am

    /*
        2*vecSum[i] '   -2*AmDer*vecSum[i]   vecSumDer[i] * 2
        -----------  =  ------------------- + ---------------
            Am                Am*Am                  Am
     */
    //--------derivative of log((Zalfa+(1+rad2)*Bm)/(Zalfa+(1-rad2)*Bm)

   /*
       Zalfa+(1+rad2*Bm) '
   log ------------------ = (log (Zalfa + (1+rad2*Bm) - log (Zalfa + (1 - rad2*Bm))'
       Zalfa+(1-rad2*Bm)

    [Zalfa+(1+rad2*Bm)]'      [Zalfa+(1-rad2*Bm)]'
=  --------------------- --  -------------------
     Zalfa+(1+rad2*Bm)        Zalfa+(1-rad2*Bm)

     ZalfaDer + rad2*BmDer      ZalfaDer - rad2*BmDer
 =  --------------------- --  ----------------------
      Zalfa+(1+rad2*Bm)        Zalfa+(1-rad2*Bm)

    */

    /* Derivative of   (Am/(2*rad2*Bm))*( (2.0/Am)*vecSum[i] -
                       Bi[i]/Bm)*log((Zalfa+(1+rad2)*Bm)/(Zalfa+(1-rad2)*Bm));

       (f * g * h)' = f' * g * h + h * g' * f + f * g * h'
       [f * (g - h * l)]' = f'(g-h * l) + f (g' -h' * l - l * h)

       WHERE f = Am/(2*rad2*Bm)) <===============> Am_rad2_Bm ------ double f
             g = (2.0/Am)*vecSum[i] <====>  Am_VecSum[i] ----- vector<double> g(nc)
             h = Bi[i]/Bm      <====> Bi_BmDer[i] ------ vector<double> h(nc)
             l = log((Zalfa+(1+rad2)*Bm)/(Zalfa+(1-rad2)*Bm)) <===> logZalfa_rad2 ----- double l
*/

      double f = Am/(2*rad2*Bm);
      double fder = Am_rad2_Bm;
      vector<double> g(nc);
      vector<double> gder(nc);
      vector<double> h(nc);
      vector<double> hder(nc);
      double l = logZalfa_rad2;

      vector<double> second_term_derivative(nc);
      vector<double> first_term_derivative(nc);


    for (int i=0; i<nc;i++){

        Bi[i]=press*bi[i]/(R*T);
        BiDer[i]=bi[i]/(R*T);

        //BmDer = bm/(R*T);
        Bi_BmDer[i] = BiDer[i]/Bm-(1/(Bm*Bm))*Bi[i];
        Am_VecSum[i] = (-2*AmDer*vecSum[i])/(Am*Am) + 2*vecSumDer[i]/Am;

        //f'(g-h * l) + f (g' -h' * l - l * h)
    second_term_derivative[i] = fder*(g[i] - h[i] * l) + f * (gder[i] - hder[i] * l - l * h[i]);
    ffi[i] =(Bi[i]/Bm)*(Zalfa-1.0)-log(Zalfa-Bm)-
    (Am/(2*rad2*Bm))*( (2.0/Am)*vecSum[i] -
                       Bi[i]/Bm)*log((Zalfa+(1+rad2)*Bm)/(Zalfa+(1-rad2)*Bm));

     //2*vecSum[i] '   2  * vecSum[i] + vecSumDer[i] * 2
    //-----------  = (---) '                         ------
    //    Am          Am                                Am

    /*
        //2*vecSum[i] '   -2*AmDer*vecSum[i]   vecSumDer[i] * 2
        //-----------  =  ------------------- + ---------------
        //    Am                Am*Am                  Am
     */

    fug[i] = exp(ffi[i]);

    }

  return {0};
  }



//vector<double> PropertyPackage::

//vector<double> PropertyPackage::


//vector<double> PropertyPackage::
 
double PropertyPackage::mixRule(vector<double> xmol, vector<double> property){

double propMixture=0;


for (int i=0; i<nc; i++){


propMixture = propMixture + xmol[i]*property[i];

}






return propMixture;
}























