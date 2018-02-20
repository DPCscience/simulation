###http://r.789695.n4.nabble.com/how-to-generate-a-random-correlation-matrix-with-restrictions-td903519.html
#Method 2
R <- matrix(runif(16), ncol=4) 
R <- (R * lower.tri(R)) + t(R * lower.tri(R)) 
diag(R) <- 1 
eigen(R)$val 
library(Matrix)
Q <- nearPD(R, posd.tol=1.e-04)$mat 
eigen(Q)$val 
max(abs(Q - R))  # maximum discrepancy between R and Q 

#Method 1
R <- matrix(runif(36), ncol=6) 
RtR <- R %*% t(R) 
Q <- cov2cor(RtR) 

#https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
#Method 3, creates a 1 in the diagonal unlike the other methods
d = 6;    
k = 6;      
W = randn(d,k)
S <- W%*%t(W)+ + diag(rand(1,d))
S <- diag(1./sqrt(diag(S))) %*% S %*% diag(1./sqrt(diag(S)))
S
