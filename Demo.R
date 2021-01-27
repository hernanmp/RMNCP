rm(list = ls())

library(RMNCP)



library(MASS)
library(ks)
library(kdensity)

#RMNCP::MNWBS_full
#Lenght of time horizon
T= 300

#Location of change points
v =  c(floor(T/3),2*floor(T/3))
p=20

#Matrix for data
y =  matrix(0,T,p)

#Means of the data
mu0 =   rep(0,p)
mu1 = rep(0,p)
mu1[1:floor(p/2)] =  1

#Covariance matrices of the data
Sigma0 =  diag(p)
Sigma1 = diag(p)


#Generate data

for(t in 1:T)
{
  if(t <v[1] ||  t > v[2])
  {
    y[t,] = mvrnorm(n = 1, mu0, Sigma0)
  }
  
  if(t >=v[1] &&  t < v[2])
  {
    y[t,] = mvrnorm(n = 1, mu1,Sigma1)
  }
}## close for generate data


#Random intervals
M =   50
alpha =  sample.int(size =M  , n = T,replace = TRUE)
beta =   sample.int(size =M  , n = T,replace = TRUE)#alpha + floor((T- alpha)*runif(M))
#
for(j in 1:M)
{
  aux =  alpha[j]
  aux2 =  beta[j]
  #
  alpha[j] = min(aux,aux2)
  beta[j] = max(aux,aux2)
}
#Bandwith
K_max = 30
h = 5*(K_max*log(T)/T)^{1/p}  

#Run method
S =MNWBS_full(y,y,alpha,beta,h)
S

