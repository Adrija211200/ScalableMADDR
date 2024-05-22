setwd("C:/Users/Adrija Saha/Desktop/Final Sem Project")
set.seed(2203)
library(Rcpp)
sourceCpp("modified.cpp")
library(proxy)
library(pracma)

# Function to find an orthonormal basis for the subspace of V orthogonal to e
orthogonalize_to_e <- function(V, e) {
  status <- T
  j <- 1; v_0 <- NULL
  while(status){
    if(abs(t(V[, j])%*%e) > 1e-10){
      status <- F
      v_0 <- V[, j]
    }
    j<-j+ 1 
  }
  v_0<-as.matrix(v_0)
  # Compute the new matrix W
  W <- apply(V, 2, function(v) {
    v <- as.matrix(v)
    cont<- (t(v) %*% (as.matrix(e)) / (t(v_0) %*% as.matrix(e)))[1,1]
    new_col <- v - cont * v_0
    return(new_col)
  })
  
  # Remove columns with all zeros
  W <- W[, colSums(W != 0) > 0]
  
  # Apply Gram-Schmidt orthogonalization
  result <- gramSchmidt(as.matrix(W))
  
  # The orthogonal basis is the first n columns of Q
  orthogonal_basis <- result$Q
  
  return(orthogonal_basis)
}

### Computing the elementary symmetric polynomials
compute_elementary_symmetric <- function(k, eigenvalues) {
  N <- length(eigenvalues)
  
  # Initialize the array to store the elementary symmetric polynomials
  e <- matrix(0, nrow <- N + 1, ncol <- k+1)
  e[,1] <- 1
  
  # Compute elementary symmetric polynomials using the given algorithm
  for (l in 2:(k+1)) {
    for (n in 2:(N+1)) {
      e[n , l ] <- e[n-1, l ] + eigenvalues[n-1] * e[n-1, l-1]
    }
  }
  
  # Output the k-th elementary symmetric polynomial
  return(e)
}

### Sampling k eigenvectors

sample_k_eigenvectors <- function(k, eigenvalues) {
  N <- length(eigenvalues)
  
  e_k_n <- compute_elementary_symmetric(k, eigenvalues)
  
  # Initialize variables
  J <- numeric(0)
  l <- k+1
  
  # Sampling eigenvectors
  for (n in (N+1):2) {
    if (l == 1) {
      break
    }
    
    if (runif(1) < eigenvalues[n-1] * e_k_n[n-1,l-1] / e_k_n[n,l]){
      J <- c(J, n-1)
      l <- l - 1
    }
  }
  
  # Output the set of sampled indices
  return(J)
}

### Sampling from a k-DPP

sample_Y_from_eigendecomposition <- function(k,eigendecomposition) {
  
  
  # Extract eigenvectors and eigenvalues
  eigenvalues <- eigendecomposition$values
  eigenvectors <- eigendecomposition$vectors
  N <- length(eigenvalues)  
  
  J <- sample_k_eigenvectors(k,eigenvalues)
  V <- as.matrix(eigenvectors[,J])
  Y <- numeric(0)
  
  # Sampling set Y
  while (ncol(V) > 0) {
    sample_prob <- rowSums(V^2)/ncol(V)
    i <- sample(1:nrow(V), 1,prob = sample_prob)
    Y <- c(Y, i)
    e_i<- rep(0,nrow(V))
    e_i[i]<-1
    if(ncol(V)==1){
      break
    }else{
      V <- orthogonalize_to_e(V,as.matrix(e_i))}
    
  }
  return(Y)
}
### Data loading

library(MASS)

# Set seed for reproducibility
set.seed(2203)

# Define the number of samples and features
num_samples <- 60
num_features <- 50
num_classes <- 3

# Create mean vectors for each class
class_means <- matrix(c(-0.5,0,0.5), nrow <- num_classes,ncol <- num_features)

# Generate data using multivariate normal distribution with varying mean
generate_data <- function(num_samples, class_means) {
  class_labels <- rep(1:num_classes, each <- num_samples/num_classes)
  data <- matrix(0, nrow <- num_samples, ncol <- num_features)
  
  for (i in 1:num_samples) {
    class_label <- class_labels[i]
    mean_vector <- class_means[class_label, ]
    data[i, ] <- mvrnorm(1, mu <- mean_vector, Sigma <- diag(num_features))
  }
  
  return(data)
}
# Custom KNN function using MADD distance metric
knn_madd <- function(train_data, test_data, train_labels, k,sampled_data_matrix) {
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  predictions <- character(n_test)
  
  for (i in 1:n_test) {
    distances <- numeric(n_train)
    
    # Calculate MADD distances between the test point and all training points
    for (j in 1:n_train) {
      distances[j] <- calc_MADD(test_data[i, ], train_data[j, ],as.matrix(sampled_data_matrix))
    }
    
    # Find the indices of the k-nearest neighbors
    nearest_neighbors <- order(distances)[1:k]
    
    # Make prediction based on majority class of k-nearest neighbors
    majority_class <- table(train_labels[nearest_neighbors])
    predictions[i] <- names(majority_class)[which.max(majority_class)]
  }
  
  return(as.factor(predictions))
}

knn_madd_optimized <- function(train_data, test_data, train_labels, k,sampled_data_matrix) {
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  predictions <- character(n_test)
  
  result_matrix <- outer(1:nrow(test_data), 1:nrow(train_data), 
                         Vectorize(function(i, j){ 
                           calc_MADD(test_data[i, ], train_data[j, ], 
                                     as.matrix(sampled_data_matrix))}))
  predictions <- apply(result_matrix,1,function(d){
    nearest_neighbors<- order(d)[1:k]
    majority_class <- table(train_labels[nearest_neighbors])
    return(names(majority_class)[which.max(majority_class)])
  })
  
  return(as.factor(predictions))
}


R<- 100

result <- replicate(R,{
  # Generate data
  simulated_data <- generate_data(num_samples, class_means)
  
  # Add class_labels to the dataset
  
  simulated_data <- data.frame(simulated_data, class <- 
                                 factor(rep(1:num_classes, each <- num_samples/num_classes)))
  
  
  # Generating Test Data
  # Generate data
  test_data <- generate_data(450, class_means)
  
  # Add class_labels to the dataset
  
  test_data <- data.frame(test_data, class <- 
                            factor(rep(1:num_classes, each <- 450/num_classes)))
  
  
  
  # Choice of L-ensemble matrix
  data_without_label <- as.matrix(simulated_data[,-ncol(simulated_data)])
  test_data_without_label <- as.matrix(test_data[,-ncol(test_data)])
  
  # Custom KNN function using MADD distance metric
  T1 <- system.time({# Sample from 1st population
    A1 <- as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]);
    L1 <- as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)])%*%t(
      as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]));
    result1<-eigen(L1);
    sample1 <- sample_Y_from_eigendecomposition(k=16,result1);
    # Sample from 2nd population
    A2 <- as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]);
    L2 <- as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)])%*%t(
      as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]));
    result2<-eigen(L2);
    sample2 <- sample_Y_from_eigendecomposition(k=16,result2);
    # Sample from 3rd population
    A3 <- as.matrix(simulated_data[simulated_data$class==3,-ncol(simulated_data)]);
    L3 <- as.matrix(simulated_data[simulated_data$class==3,-ncol(simulated_data)])%*%t(
      as.matrix(simulated_data[simulated_data$class==3,-ncol(simulated_data)]));
    result3<-eigen(L3);
    sample3 <- sample_Y_from_eigendecomposition(k=16,result3);
    sampled_data <- rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample1,],
                         simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample2,],
                         simulated_data[simulated_data$class==3,-ncol(simulated_data)][sample3,]);
    
      test_prediction <- 
                      knn_madd(data_without_label,test_data_without_label,
                               simulated_data$class,5,as.matrix(sampled_data))})[3]
  M1<- table(test_data$class,test_prediction)  
  mis_prob1 <- 1-sum(diag(M1))/sum(M1)
  
  T2<- system.time({# Sample from 1st population
    A1 = as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]);
    
    L1 = exp(-(as.matrix(dist(A1, by_rows = T,upper = T,diag = T))^2)/num_features);
    result1=eigen(L1);
    sample1 <- sample_Y_from_eigendecomposition(k=16,result1);
    # Sample from 2nd population
    A2 = as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]);
    L2 = exp(-(as.matrix(dist(A2, by_rows = T,upper = T,diag = T))^2)/num_features);
    result2=eigen(L2);
    sample2 <- sample_Y_from_eigendecomposition(k=16,result2);
    # Sample from 3rd population
    A3 = as.matrix(simulated_data[simulated_data$class==3,-ncol(simulated_data)]);
    
    L3 = exp(-(as.matrix(dist(A3, by_rows = T,upper = T,diag = T))^2)/num_features);
    result3=eigen(L3);
    sample3 <- sample_Y_from_eigendecomposition(k=16,result3);
    
    sampled_data = rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample1,],
                         simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample2,],
                         simulated_data[simulated_data$class==3,-ncol(simulated_data)][sample3,]);
    
    test_prediction <- knn_madd(data_without_label,test_data_without_label,
                                simulated_data$class,5,as.matrix(sampled_data))})[3]
  M2<- table(test_data$class,test_prediction)  
  
  mis_prob2 <- 1-sum(diag(M2))/sum(M2)
  
  ### Now just trying with random sample instead of Determinantal Point Process
  T3 <- system.time({sample_random1 <- sampleWithoutReplacement((num_samples/3),16);
                    sample_random2 <- sampleWithoutReplacement((num_samples/3),16);
                    
                    sample_random3 <- sampleWithoutReplacement((num_samples/3),16);
                    
                    
                    sampled_data_random <- rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample_random1,],
                                                simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample_random2,],
                                                simulated_data[simulated_data$class==3,-ncol(simulated_data)][sample_random3,]);
                    
                    
    test_prediction_random <- knn_madd(data_without_label,
                                                                 test_data_without_label,
                                                                 simulated_data$class,5,
                                                                 as.matrix(sampled_data_random))})[3]
  M3<- table(test_data$class,test_prediction_random) 
  
  mis_prob3 <- 1-sum(diag(M3))/sum(M3)
  #library(class)
  #T4<- system.time(knn_prediction <- knn(data_without_label,
  #                                       test_data_without_label,simulated_data$class,k=5))[3]
  #M4 <- table(test_data$class,knn_prediction) 
  #mis_prob4 <- 1-sum(diag(M4))/sum(M4)
  ### Now trying with traditional MADD
  
  #T5<- system.time(traditional_prediction <- knn_madd(data_without_label,
  #                test_data_without_label,simulated_data$class,5,data_without_label))[3]
  #M5<-table(test_data$class,traditional_prediction) 
  #mis_prob5 <- 1-sum(diag(M5))/sum(M5)
  c(print(R),mis_prob1,mis_prob2,mis_prob3,#mis_prob4,mis_prob5,
    T1,T2,T3#,T4,T5
  )
})
result<-result[-1,]
mean_result <- apply(result,1,mean)
se_result <- apply(result,1,sd)
mean_result
se_result
