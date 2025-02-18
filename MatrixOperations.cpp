#include "MatrixOperations.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <functional>
#include <string>

const double EPSILON = 1e-9; // Threshold for considering values as zero
using namespace std;

void printMatrix(vector<vector<double>> mat) {

 for(const auto& row : mat) {
   for (double val : row) {
         cout<<setw(10)<<setprecision(10)<<val<< " ";
     
  }  
   cout<< endl;
 }

}

void printVector(vector<double> vec){


for (int i=0; i< vec.size(); i++){

    cout<<setw(10)<<setprecision(6)<<vec[i]<< " ";
 }


}

double vecSum(vector<double> vec){

double sum =0;

int n=vec.size();

for (int i=0; i<n;i++){


sum = sum+vec[i];

}

return sum;
}

vector<double> vecDiff(vector<double> vec1, vector<double> vec2){

int n1 = vec1.size();
int n2 = vec2.size();

vector<double> result(n1);

for (int i=0; i<n1; i++){

result[i]= vec1[i]-vec2[i];

}

return result;

}

double euclidean_norm(vector<double> vec){
double sum = 0.0;

for (int i=0; i<vec.size(); i++){


sum = sum+vec[i]*vec[i];

}



return sqrt(sum);
}

vector<double> prodScal(vector<double> vec, double scal){

int n= vec.size();
vector<double> result(n);

for(int i=0; i<n; i++){



result[i] = scal * vec[i]; 

}

return result;

}



vector<double> vecTimesMat(vector<vector<double>> mat, vector<double> vec){

int v_size = vec.size();
vector<double> result(v_size);

  for (int i=0; i < v_size; i++) {
    for (int j=0; j< mat[i].size(); j++){
      result[i] += mat[i][j]*vec[j];
     }

  }


return result; 


}
//Updated according to java codebase
//not used for calc the sum of a vector

vector<double> prodVec(vector<double> vec1, vector<double>vec2){

int n = vec1.size();

vector<double> result (n);

for (int i=0; i<n; i++){

    result[i] = vec1[i]*vec2[i];
  }

return result;
}

vector<double> divVec(vector<double> vec1, vector<double> vec2){


int n = vec1.size();
vector<double> result(n);


for(int i=0; i<n; i++){

result[i]=vec1[i]/vec2[i];


}


return result;

}

vector<double> addScal(vector<double> vec, double scal){


int n = vec.size();
vector<double> result(n);


for(int i=0; i<n; i++){

result[i]=vec[i]+scal;


}


return result;


}

std::vector<std::vector<double>> invertMatrix(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(2 * n, 0.0));

    // Create the augmented matrix by appending the identity matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = matrix[i][j];
        }
        augmentedMatrix[i][n + i] = 1.0;
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Find the pivot
        int pivotRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::fabs(augmentedMatrix[j][i]) > std::fabs(augmentedMatrix[pivotRow][i])) {
                pivotRow = j;
            }
        }

        // Swap rows if necessary
        if (std::fabs(augmentedMatrix[pivotRow][i]) < EPSILON) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }
        std::swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);

        // Normalize the pivot row
        double pivotValue = augmentedMatrix[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            augmentedMatrix[i][j] /= pivotValue;
        }

        // Eliminate the current column in other rows
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double factor = augmentedMatrix[j][i];
                for (int k = 0; k < 2 * n; ++k) {
                    augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                }
            }
        }
    }

    // Extract the inverse matrix from the augmented matrix
    std::vector<std::vector<double>> inverseMatrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverseMatrix[i][j] = augmentedMatrix[i][n + j];
        }
    }

    return inverseMatrix;
}


vector<vector<double>> computeJacobian(const function<vector<double>(const vector<double>&)>& equations,
const vector<double>& vars) {
double stepSize = 1.0e-6;

int numEquations = equations(vars).size();
int numVariables = vars.size();

vector<vector<double>> jacobian (numEquations, vector<double>(numVariables,0.0));

for (int j=0; j< numVariables; ++j){


vector<double> perturbedVars = vars;
perturbedVars[j] += stepSize;

vector<double> fPlus = equations(perturbedVars);
vector<double> fOriginal = equations(vars);

for (int i=0; i< numEquations; ++i){

jacobian[i][j] = (fPlus[i]-fOriginal[i])/stepSize;

  }

}

return jacobian;

}

//return either min or max from vector of double elements



double find_min_max(std::vector<double>& nums, string m,int start, int end) {

if (m == "min") {
    auto min = min_element(nums.begin(),nums.end());

    return (double)(*min);
}
    else if (m == "max") {
        auto max = max_element(nums.begin(),nums.end());
        return (double)(*max);
    }
else {
        throw std::invalid_argument("Invalid input: choose either min or max");
    }

}


vector<double> getVecSum(vector<double> xmol, vector<vector<double>> mat){
          int n = xmol.size();
          vector<double> xSum (n);
          vector<vector<double>> newMat(n,vector<double>(n, 0.0));
          double sum = 0.0;

           for (int i=0; i< n; i++){
             newMat[i] = prodVec(xmol,mat[i]);
           }
                for (int i=0; i< n; i++){
                    xSum[i] = vecSum(newMat[i]);
                }

           return xSum;
            }

int indexOf(vector<double> v, double element){
  int index = -1;
  auto it = lower_bound(v.begin(), v.end(), element);
  index = (int) (it - v.begin());
  return index;
  }
