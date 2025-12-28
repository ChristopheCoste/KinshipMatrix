### Inputs
F = matrix(c(0.23, 0.51, 0.84, 0, 0, 0, 0 ,0 ,0),nrow=3,byrow=T) #fertility component of projection model
S = matrix(c(0, 0, 0, .57 ,0 ,0,0, .54 ,.46),nrow=3,byrow=T)   #survival component of projection model
gmax=3  #number of generations investigated
# same-litter sisters:
#Z= ?   #same-litter-sister matrix, enter it here (and uncomment) if there is more than one class of newborn and reproduction is neither Bernoulli nor Poisson
# Particular reproductive cases
Fstar=F # Poisson reproduction (comment/uncomment if reproduction is not/is Poisson)
#Fstar=matrix(0, dim(F)[1],dim(F)[1]) #Bernoulli reproduction (comment/uncomment if reproduction is not/is Bernoulli)
#Fstar= ? #% if reproduction is neither Poisson nor Bernoulli and there is sonly one newborn class

### eigen properties of the population projection matrix
M=F+S #projection matrix
ev <- eigen(M);lam <- max(Re(ev$values));
w_raw <- abs(Re(ev$vector))[, which.max(Re(ev$values))]
w <- w_raw/sum(w_raw)

### same-litter sisters
Dw <- diag(w)
Zdagger=Fstar%*%Dw%*%t(F)
#Z=lam*Z%*%Dw  #uncomment if NOT in one of the particular cases

###  Generation structured matrices/tensors
s=dim(F)[1]
myDiag <- function(x, vec, k) {  x[row(x) == col(x) - k] <- vec 
  x } #turns a vector (vec) into a zero matrix matrix but with vec on the kth diagonal
J1<-myDiag(matrix(0, gmax, gmax),  rep(1, gmax-1), -1) 
J2=matrix(0, gmax, gmax);J2[2,2]=1;
J3=matrix(0, gmax, gmax);J3[1,1]=1;
J=J3%x%diag(s)
Jy <- diag(gmax)%x%matrix(1, nrow = 1, ncol = s)
bigF=J1%x%F
bigS=diag(gmax)%x%S
bigM=bigF+bigS
bigZdagger=J2%x%Zdagger
bigDw <- diag(gmax)%x% Dw

### kinship matrix
vec <- function(x) {  y <- as.vector(x)  
return(y) }

vecKw=solve(lam*diag(s*s*gmax*gmax)-bigM%x%bigM)%*%((bigS%x% bigF   + bigF%x% bigS)%*%vec(J*bigDw) + vec(bigZdagger) )
Kw <- matrix(vecKw, nrow = gmax * s, ncol = gmax * s)
K <- matrix(vecKw, nrow = gmax * s, ncol = gmax * s) %*% solve(bigDw)
Ku <- matrix((Jy%x%Jy)%*%vecKw, nrow = gmax, ncol = gmax)

round(K, digits = 2)
round(Ku, digits = 2)
