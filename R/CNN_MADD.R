setwd("C:/Users/Adrija Saha/Desktop/Final Sem Project")
set.seed(2203)
library(Rcpp)
sourceCpp("modified.cpp")
library(proxy)
library(pracma)
condensed_set <- function(train_data, train_labels, threshold, k) {
  # Calculate Euclidean distance matrix
  dist_matrix <- matrix(0,nrow = nrow(train_data),ncol = nrow(train_data))
  for (i in 1:nrow(train_data)){
    for (j in 1:i){
      dist_matrix[i,j] <- calc_MADD(train_data[i,],train_data[j,],
                                    as.matrix(train_data))
      dist_matrix[j,i] <- dist_matrix[i,j]
    }
  }
  # Initialize empty list to store subset indices
  subset_indices <- NULL
  dp <- rep(0,nrow(train_data))
  dn <-rep(0,nrow(train_data))
  # Iterate through each row of the train data
  for (i in 1:nrow(train_data)) {
    # Calculate distance of the current observation to all other observations
    distances <- dist_matrix[i, ]
    
    # Sort distances in ascending order and get indices of k-nearest neighbors
    nearest_indices <- order(distances)[2:(k+1)]
    ifelse(train_labels[i]==train_labels[nearest_indices],
    dp[nearest_indices] <- dp[nearest_indices]+1,
    dn[nearest_indices] <- dn[nearest_indices]+1)
    
  }
  for (j in 1:nrow(train_data)){
  # Check if conditions are met
  if (dp[j] + dn[j] > 0 &&
            dp[j] - dn[j] > threshold){
    subset_indices <- c(subset_indices,j)
  }
  }
  # Return subset of train data
  return(list(condensed_data = train_data[subset_indices, ],
              condensed_label = train_labels[subset_indices] ))
}

# Function to find an orthonormal basis for the subspace of V orthogonal to e
orthogonalize_to_e <- function(V, e) {
  status = T
  j = 1; v_0 = NULL
  while(status){
    if(abs(t(V[, j])%*%e) > 1e-10){
      status = F
      v_0 = V[, j]
    }
    j=j+ 1 
  }
  v_0=as.matrix(v_0)
  # Compute the new matrix W
  W <- apply(V, 2, function(v) {
    v <- as.matrix(v)
    cont= (t(v) %*% (as.matrix(e)) / (t(v_0) %*% as.matrix(e)))[1,1]
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
  e <- matrix(0, nrow = N + 1, ncol = k+1)
  e[,1] <- 1
  
  # Compute elementary symmetric polynomials using the given algorithm
  for (l in 2:(k+1)) {
    for (n in 2:(N+1)) {
      e[n , l ] = e[n-1, l ] + eigenvalues[n-1] * e[n-1, l-1]
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
    e_i= rep(0,nrow(V))
    e_i[i]=1
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
num_samples <- 150
num_features <- 100
num_classes <- 2

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


# Create mean vectors for each class
class_means <- matrix(c(0,0.50), nrow = num_classes,ncol = num_features)

# Generate data using multivariate normal distribution with varying mean
generate_data <- function(num_samples, class_means) {
  class_labels <- rep(1:num_classes, each = num_samples/num_classes)
  data <- matrix(0, nrow = num_samples, ncol = num_features)
  
  for (i in 1:num_samples) {
    class_label <- class_labels[i]
    mean_vector <- class_means[class_label, ]
    data[i, ] <- mvrnorm(1, mu = mean_vector, Sigma = diag(num_features))
  }
  
  return(data)
}
R = 100
result <- replicate(R,{
  # Generate data
  simulated_data <- generate_data(num_samples, class_means)
  
  # Add class_labels to the dataset
  simulated_data <- data.frame(simulated_data, class = 
                                 factor(rep(1:num_classes, each = num_samples/num_classes)))
  
  # Generating Test Data
  # Generate data
  test_data <- generate_data(500, class_means)
  
  # Add class_labels to the dataset
  test_data <- data.frame(test_data, class = 
                            factor(rep(1:num_classes, each = 500/num_classes)))
  
  
  # Choice of L-ensemble matrix
  data_without_label = as.matrix(simulated_data[,-ncol(simulated_data)])
  test_data_without_label = as.matrix(test_data[,-ncol(test_data)])
  
  
  # Custom KNN function using MADD distance metric
  T1 <- system.time({# Sample from 1st population
    A1 = as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)])
    L1 = as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)])%*%t(
      as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]))
    result1=eigen(L1)
    
    # Sample from 2nd population
    A2 = as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)])
    L2 = as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)])%*%t(
      as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]))
    result2=eigen(L2)
    
    result2=eigen(L2);
    sample2 <- sample_Y_from_eigendecomposition(k=32,result2);
    sampled_data = rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample1,],
                         simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample2,]);
    
    test_prediction <- knn_madd(data_without_label,test_data_without_label,
                                simulated_data$class,5,as.matrix(sampled_data))})[3]
  M1<- table(test_data$class,test_prediction)  
  
  mis_prob1 <- 1-sum(diag(M1))/sum(M1)
  T2 <- system.time({# Sample from 1st population
    A1 = as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]);
    
    L1 = exp(-(as.matrix(dist(A1, by_rows = T,upper = T,diag = T))^2)/num_features);
    result1=eigen(L1);
    sample1 <- sample_Y_from_eigendecomposition(k=32,result1);
    # Sample from 2nd population
    A2 = as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]);
    L2 = exp(-(as.matrix(dist(A2, by_rows = T,upper = T,diag = T))^2)/num_features);
    result2=eigen(L2);
    sample2 <- sample_Y_from_eigendecomposition(k=4,result2);
    sampled_data = rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample1,],
                         simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample2,]);
    
    test_prediction <- knn_madd(data_without_label,test_data_without_label,
                                simulated_data$class,5,as.matrix(sampled_data))})[3]
  M2<- table(test_data$class,test_prediction)  
  
  mis_prob2 <- 1-sum(diag(M2))/sum(M2)
  T3 <- system.time({# Sample from 1st population
    A1 = as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]);
    
    L1 = exp(-(as.matrix(dist(A1, by_rows = T,upper = T,diag = T))^2)/num_features);
    result1=eigen(L1);
    sample1 <- sample_Y_from_eigendecomposition(k=8,result1);
    # Sample from 2nd population
    A2 = as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]);
    L2 = exp(-(as.matrix(dist(A2, by_rows = T,upper = T,diag = T))^2)/num_features);
    result2=eigen(L2);
    sample2 <- sample_Y_from_eigendecomposition(k=8,result2);
    sampled_data = rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample1,],
                         simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample2,]);
    
    test_prediction <- knn_madd(data_without_label,test_data_without_label,
                                simulated_data$class,5,as.matrix(sampled_data))})[3]
  M3<- table(test_data$class,test_prediction)  
  
  mis_prob3 <- 1-sum(diag(M3))/sum(M3)
  
  T4 <- system.time({# Sample from 1st population
    A1 = as.matrix(simulated_data[simulated_data$class==1,-ncol(simulated_data)]);
    
    L1 = exp(-(as.matrix(dist(A1, by_rows = T,upper = T,diag = T))^2)/num_features);
    result1=eigen(L1);
    sample1 <- sample_Y_from_eigendecomposition(k=16,result1);
    # Sample from 2nd population
    A2 = as.matrix(simulated_data[simulated_data$class==2,-ncol(simulated_data)]);
    L2 = exp(-(as.matrix(dist(A2, by_rows = T,upper = T,diag = T))^2)/num_features);
    result2=eigen(L2);
    sample2 <- sample_Y_from_eigendecomposition(k=16,result2);
    sampled_data = rbind(simulated_data[simulated_data$class==1,-ncol(simulated_data)][sample1,],
                         simulated_data[simulated_data$class==2,-ncol(simulated_data)][sample2,]);
    
    test_prediction <- knn_madd(data_without_label,test_data_without_label,
                                simulated_data$class,5,as.matrix(sampled_data))})[3]
  M4<- table(test_data$class,test_prediction)  
  
  mis_prob4 <- 1-sum(diag(M4))/sum(M4)
  
  T5 <- system.time({
    condensed_result <- condensed_set(data_without_label,simulated_data$class,threshold = 1,k = 5)
    condensed_train <- condensed_result$condensed_data
    condensed_train_label <- condensed_result$condensed_label
    test_prediction <- knn_madd(condensed_train,test_data_without_label,
                                condensed_train_label,5,as.matrix(condensed_train))
  })[3]
  M5<- table(test_data$class,test_prediction)  
  
  mis_prob5 <- 1-sum(diag(M5))/sum(M5)
  
  T6 <- system.time({
    condensed_result <- condensed_set(data_without_label,simulated_data$class,threshold = 2,k = 5)
    condensed_train <- condensed_result$condensed_data
    condensed_train_label <- condensed_result$condensed_label
    test_prediction <- knn_madd(condensed_train,test_data_without_label,
                                condensed_train_label,5,as.matrix(condensed_train))
  })[3]
  M6<- table(test_data$class,test_prediction)  
  
  mis_prob6<- 1-sum(diag(M6))/sum(M6)
  
  T7 <- system.time({
    condensed_result <- condensed_set(data_without_label,simulated_data$class,threshold = 3,k = 5)
    condensed_train <- condensed_result$condensed_data
    condensed_train_label <- condensed_result$condensed_label
    test_prediction <- knn_madd(condensed_train,test_data_without_label,
                                condensed_train_label,5,as.matrix(condensed_train))
  })[3]
  M7<- table(test_data$class,test_prediction)  
  
  mis_prob7 <- 1-sum(diag(M7))/sum(M7)
  
  T8 <- system.time({
    condensed_result <- condensed_set(data_without_label,simulated_data$class,threshold = 4,k = 5)
    condensed_train <- condensed_result$condensed_data
    condensed_train_label <- condensed_result$condensed_label
    test_prediction <- knn_madd(condensed_train,test_data_without_label,
                                condensed_train_label,5,as.matrix(condensed_train))
  })[3]
  M8<- table(test_data$class,test_prediction)  
  
  mis_prob8 <- 1-sum(diag(M8))/sum(M8)  
  library(class)
  T9<- system.time(knn_prediction <- knn(data_without_label,
                                         test_data_without_label,simulated_data$class,k=5))[3]
  M9 <- table(test_data$class,knn_prediction) 
  mis_prob9 <- 1-sum(diag(M9))/sum(M9)
  ### Now trying with traditional MADD
  
  T10<- system.time(traditional_prediction <- knn_madd(data_without_label,
              test_data_without_label,simulated_data$class,5,data_without_label))[3]
  M10<-table(test_data$class,traditional_prediction) 
  mis_prob10 <- 1-sum(diag(M10))/sum(M10)
  c(print(R),mis_prob1,mis_prob2,
    mis_prob3,
    mis_prob4,mis_prob5,mis_prob6,
    mis_prob7, mis_prob8,mis_prob9,mis_prob10,
    T1 ,T2, T3,T4, T5, T6, T7, T8, T9, T10
  )
})
result<-result[-1,]
#rownames(result) <- c("mis_prob1",
#"mis_prob2",
#  #"mis_prob3",
# "mis_prob4",
#  "T1"#,"T2",
#"T3",
# "T4"
#)
mean_result = apply(result,1,mean)
se_result = apply(result,1,sd)
mean_result
se_result
