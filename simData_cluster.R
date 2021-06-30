# Clustered Dataset Simulation

## 1. Define Parameter, Generate PI/P
G = 10; K = 2; J = 10; N = 1000
set.seed(1080)
PI = matrix(runif(G*K, 0, 1), nrow = G, ncol = K)
PI[,K] = 1 - rowSums(PI[,1:(K-1),drop = F])
P = matrix(0.2, nrow = J, ncol = K)
P[,2] = 1 - P[,1]
## 2. Generate Dataset
### 1 Cluster  + 1 latent class + J causes
data = matrix(0, nrow = N, ncol = 2 + J ) 
data[,1] = sample(1:G, N, replace = T, prob = rep(1/G,G))
for (i in 1:N){
  cls_i = data[i, 1]
  # generate latent class
  data[i,2] = sample(1:K, 1, prob = PI[cls_i,,drop = F] ) 
  # generate causes
  for (j in 1:J){
    data[i,j+2] = rbinom(n=1, size = 1, prob = P[j, data[,2]])
  }
  # The outcome data to be generated 
}
colnames(data) = c("G","U",paste("A",1:J,sep=""))

