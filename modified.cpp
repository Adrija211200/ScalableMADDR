#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <cstdlib>
using namespace Rcpp;

// [[Rcpp::export]]
bool is_column_of_matrix(const NumericVector& vector, const NumericMatrix& matrix) {
  
  for (int col = 0; col < matrix.ncol(); ++col) {
    bool isColumn = true;
    for (int row = 0; row < matrix.nrow(); ++row) {
      if (matrix(row, col) != vector[row]) {
        isColumn = false;
        break;
      }
    }
    if (isColumn) {
      return true;
    }
  }
  return false;
}

// [[Rcpp::export]]
double euclidean_distance(NumericVector x, NumericVector y) {
  if(x.size() != y.size()){
    Rprintf("Not applicable");
  }
  double dist = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    dist += (x[i] - y[i])*(x[i] - y[i]);
  }
  
  return sqrt(dist);
}

// [[Rcpp::export]]
double sum_distances(const Rcpp::NumericVector& x,const Rcpp::NumericVector& y,const Rcpp::NumericMatrix& data) {
  double sum = 0.0;
  for (int i = 0; i < data.nrow(); ++i) {
    sum += std::abs(euclidean_distance(x, data(i, Rcpp::_)) - euclidean_distance(y, data(i, Rcpp::_)));
  }
  return sum;
}


// [[Rcpp::export]]
double calc_MADD(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& data) {
  
  int status = is_column_of_matrix(x, Rcpp::transpose(data)) + is_column_of_matrix(y, Rcpp::transpose(data));
  double commonsum = sum_distances(x,y, data);
  double result = (commonsum - status*euclidean_distance(x,y))/(data.nrow() - status); 
  return result;
}

// [[Rcpp::export]]
std::vector<int> sampleWithoutReplacement(int n, int k) {
  // Initialize the vector to store the sampled indices
  std::vector<int> result(k);
  
  // Fill the vector with the first k indices
  for (int i = 0; i < k; ++i) {
    result[i] = i + 1;
  }
  
  // Randomly shuffle the vector using Fisher-Yates algorithm
  for (int i = k; i < n; ++i) {
    int j = rand() % (i + 1);
    if (j < k) {
      result[j] = i + 1;
    }
  }
  
  return result;
}