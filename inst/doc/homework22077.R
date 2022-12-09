## -----------------------------------------------------------------------------
plot(iris)

## -----------------------------------------------------------------------------
data(iris)
head(iris,15)

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
x <- 2*(1-u)^(-1/2)
hist(x, prob = TRUE,breaks = 1000,xlim = c(2,8),main = expression(f(x)==8*x^(-3)),border = "red") #density histogram of sample
y <- seq(2, 8, .001)
lines(y, 8*y^(-3))
#density curve f(x)

## -----------------------------------------------------------------------------
n <- 1e4;j<-k<-0;y <- numeric(n)
while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1) #random variate from g(.)
  if ((16/9)*x^2*(1-x) > u) {
    #we accept x
    k <- k + 1
    y[k] <- x
  }
}
hist(y, prob = TRUE,breaks = 50,border = "red",main = expression(f(x)==12*x^2*(1-x)),)
z <- seq(0, 1, .01)
lines(z, 12*z^2*(1-z))
#density curve f(x)

## -----------------------------------------------------------------------------
n <- 1e3; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda) # the length of lambda = n
hist(x,breaks = 400,xlim = c(0,5),borde = "red")

## -----------------------------------------------------------------------------
n <- 1e4; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda) # the length of lambda = n
hist(x,prob = TRUE,breaks = 400,xlim = c(0,5),borde = "red",main = expression(f(x)==64*(2+x)^(-5)))
z <- seq(0, 5, .01)
lines(z, 64*(2+z)^(-5))

## -----------------------------------------------------------------------------
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}
test<-sample(1:1e4)
ans <- 0               
for (i in 1:100) {
  w<-system.time(quick_sort(test))[1]
  ans <- ans + w    
}
a1=ans/100
a1

test<-sample(1:2e4)
ans <- 0               
for (i in 1:100) {
  w<-system.time(quick_sort(test))[1]
  ans <- ans + w    
}
a2=ans/100
a2

test<-sample(1:4e4)
ans <- 0               
for (i in 1:100) {
  w<-system.time(quick_sort(test))[1]
  ans <- ans + w    
}
a3=ans/100
a3

test<-sample(1:6e4)
ans <- 0               
for (i in 1:100) {
  w<-system.time(quick_sort(test))[1]
  ans <- ans + w    
}
a4=ans/100
a4
test<-sample(1:8e4)
ans <- 0               
for (i in 1:100) {
  w<-system.time(quick_sort(test))[1]
  ans <- ans + w    
}
a5=ans/100
a5
an<-c(a1,a2,a3,a4,a5)
tn<-c(1e4*log(1e4),2e4*log(2e4),4e4*log(1e4),6e4*log(6e4),8e4*log(8e4))
model<-lm(an~tn)
plot(tn,an)
abline(model,)

## -----------------------------------------------------------------------------
m <- 1e4; x <- runif(m, min=0, max=1)
theta.hat1 <- mean(exp(x))
theta.hat1
theta.hat2 <- mean((exp(x)+exp(1-x))/2)
theta.hat2
var1<-mean(exp(2*x)-(theta.hat1)^2)
var1
var2<-mean((exp(x)+exp(1-x))^2/4-(theta.hat2)^2)
var2
r<-(var1-var2)/var1
r


## -----------------------------------------------------------------------------
n <- rgamma(1e4,3/2,1/2)+1# the random value of f1
m <- sqrt(2)*n^(1/2)*exp(-n/2)/(n-1)^(1/2)/exp(-(n-1)/2)*gamma(3/2)/sqrt(2*pi)# the random value of g/f1
p <- mean(m)# the mean of samples of g/f1
q <- (m-p)^2
r <- mean(q)# the variance produced by f1
r

## -----------------------------------------------------------------------------
n <- runif(1e4)
n0 <- -2*log(1-n)+1# the random value of f2
m <- exp(-(1/2))/sqrt(2*pi)*n0^(1/2)# the random value of g/f2
p <- mean(m)# the mean of samples of g/f2
q <- (m-p)^2
r <- mean(q)# the variance produced by f2
r

## -----------------------------------------------------------------------------
n <-runif(2000)
n01 <- -log((5-(1-exp(-1))*(1-1)-(1-exp(-1))*n)/5)
n02 <- -log((5-(1-exp(-1))*(2-1)-(1-exp(-1))*n)/5)
n03 <- -log((5-(1-exp(-1))*(3-1)-(1-exp(-1))*n)/5)
n04 <- -log((5-(1-exp(-1))*(4-1)-(1-exp(-1))*n)/5)
n05 <- -log((5-(1-exp(-1))*(5-1)-(1-exp(-1))*n)/5)
m1 <- (1-exp(-1))/(1+n01^2)/5
m2 <- (1-exp(-1))/(1+n02^2)/5
m3 <- (1-exp(-1))/(1+n03^2)/5
m4 <- (1-exp(-1))/(1+n04^2)/5
m5 <- (1-exp(-1))/(1+n05^2)/5
p <- mean(m1)+mean(m2)+mean(m3)+mean(m4)+mean(m5)
r <- mean((m1-mean(m1))^2)+mean((m2-mean(m2))^2)+mean((m3-mean(m3))^2)+mean((m4-mean(m4))^2)+mean((m5-mean(m5))^2)#the variance got by Stratified Importance Sampling
r

## -----------------------------------------------------------------------------
m <- 10000
g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u <- runif(m)
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat <- mean(fg)
r <- mean((fg-theta.hat)^2)#the variance got by Importance Sampling
r

## -----------------------------------------------------------------------------
set.seed(22077)
n<-1000
alpha<-0.05
muhat<-numeric(2)
Z<-numeric(1000)#the times of mu in the confidence interval
for(i in 1:1000){
  x<-rlnorm(n)
  muhat[1]<-mean(log(x))-sqrt(var(log(x)))*qt(1-alpha/2,df=n-1)/sqrt(n)
  muhat[2]<-mean(log(x))+sqrt(var(log(x)))*qt(1-alpha/2,df=n-1)/sqrt(n)
  if(muhat[1]<0 && muhat[2]>0){
    Z[i]=1
  }
}
mean(Z)

## -----------------------------------------------------------------------------
count5test <- function(x,y){
  X <- x-mean(x)
  Y <- y-mean(y)
  outx <- (sum(X>max(Y))+ sum(X<min(Y)))
  outy <- (sum(Y>max(X))+ sum(Y<min(X)))
  return(as.integer(max(c(outx,outy))>5))
}#Count Five test

## -----------------------------------------------------------------------------
alpha<-0.055
Ftest <- function(x,y,n){
  varx<-var(x)
  vary<-var(y)
  if(var(x)/var(y)>=qf(alpha/2,n-1,n-1) && var(x)/var(y)<=qf(1-alpha/2,n-1,n-1)){
    return(0)
  }
  else return(1)
}# F test with alpha=0.055 and n1=n2=n

## -----------------------------------------------------------------------------
set.seed(22077)
m<-10000
n<-20
sigma1<-1
sigma2<-1.5
powercount5<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  count5test(x,y)
}))
# the power of Count Five test
powerF<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  Ftest(x,y,n)
}))
# the power of F test
print(powercount5)
print(powerF)

## -----------------------------------------------------------------------------
set.seed(22077)
m<-10000
n<-100
sigma1<-1
sigma2<-1.5
powercount5<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  count5test(x,y)
}))
# the power of Count Five test
powerF<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  Ftest(x,y,n)
}))
# the power of F test
print(powercount5)
print(powerF)

## -----------------------------------------------------------------------------
set.seed(22077)
m<-10000
n<-1000
sigma1<-1
sigma2<-1.5
powercount5<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  count5test(x,y)
}))
# the power of Count Five test
powerF<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  Ftest(x,y,n)
}))
# the power of F test
print(powercount5)
print(powerF)

## -----------------------------------------------------------------------------
y <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
set.seed(22077)
s <- function(y,i){
  mean(y[i])
}
library(boot)
replicate <- boot(data=y,statistic=s,R=1e4)
# use bootstrap to estimate MLE
round(c(theta=replicate$t0,bias=mean(replicate$t)-replicate$t0,se=sd(replicate$t)),4)

## -----------------------------------------------------------------------------
print(boot.ci(replicate,type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
n <- 20;m=1000
sum1 <- numeric(m);sum2 <- numeric(m);sum3 <- numeric(m);sum4 <- numeric(m);
set.seed(22091)
r <- function(x,i){
  mean(x[i])
}
library(boot)
for(i in 1:m){
  x <-rnorm(n)
  obj <- boot(data=x,statistic=r,R=1e4)
  cl <- boot.ci(obj,type=c("norm","basic","perc","bca")) # Compute 95% bootstrap confidence intervals
  if(cl$norm[2]<0 && cl$norm[3]>0){
    sum1[i] <-1
  }
  if(cl$basic[4]<0 && cl$basic[5]>0){
    sum2[i] <-1
  }
  if(cl$perc[4]<0 && cl$perc[5]>0){
    sum3[i] <-1
  }
  if(cl$bca[4]<0 && cl$bca[5]>0){
    sum4[i] <-1
  }
}
round(c(norm=mean(sum1),basic=mean(sum2),perc=mean(sum3),bca=mean(sum4)),3)


## -----------------------------------------------------------------------------
rm=list(c())
#write the function
pca=function(x,i){
  val=eigen(cov(x[i,]))$values
  return(val[1]/sum(val))
}

library(bootstrap)
n=nrow(scor)
theta.hat <- pca(scor,1:n)
theta.jack <- numeric(n)
for(i in 1:n){
theta.jack[i] <- pca(scor,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),3)

## -----------------------------------------------------------------------------
rm=list(c())
library(DAAG)
attach(ironslag)
n=length(magnetic)
e1=e2=e3=e4=c()

for (k in 1:(n-1)) {
  for (l in (k+1):n) {
    y=magnetic[-c(k,l)]
    x=chemical[-c(k,l)]
    
    J1=lm(y~x)
    yhat11=J1$coef[1]+J1$coef[2]*chemical[k]
    yhat12=J1$coef[1]+J1$coef[2]*chemical[l]
    e1=c(e1,magnetic[k]-yhat11,magnetic[l]-yhat12)
    
    J2=lm(y~x+I(x^2))
    yhat21=J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
    yhat22=J2$coef[1]+J2$coef[2]*chemical[l]+J2$coef[3]*chemical[l]^2
    e2=c(e2,magnetic[k]-yhat21,magnetic[l]-yhat22)
    
    J3=lm(log(y)~x)
    yhat31=exp(J3$coef[1]+J3$coef[2]*chemical[k])
    yhat32=exp(J3$coef[1]+J3$coef[2]*chemical[l])
    e3=c(e3,magnetic[k]-yhat31,magnetic[l]-yhat32)
    
    J4=lm(log(y)~log(x))
    yhat41=exp(J4$coef[1]+J4$coef[2]*log(chemical[k]))
    yhat42=exp(J4$coef[1]+J4$coef[2]*log(chemical[l]))
    e4=c(e4,magnetic[k]-yhat41,magnetic[l]-yhat42)
  }
}
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))
detach(ironslag)

## -----------------------------------------------------------------------------
rm=list(c())
set.seed(11111)
x=rnorm(12)
y=rnorm(12)

R=999
z=c(x,y)
K=1:24
D=numeric(R)
options(warn=-1)
D0=cor(x,y,method = "spearman")
for(i in 1:R){
  k=sample(K,size=12,replace = FALSE)
  x1=z[k]
  y1=z[-k]
  D[i]=cor(x1,y1,method = "spearman")
}
p=mean(c(D0,D)>=D0)
options(warn=0)
p
cor.test(x,y)

## -----------------------------------------------------------------------------
set.seed(22077)
# use the function below to calculate the density of  Laplace distribution
f <- function(x){
  return(exp(-abs(x))/2)
}
# use the function below to generate chain with given parameters:sigma,n,X_0,N
mh <- function(sigma,N,x0){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for(i in 2:N){
    xt <- x[i-1]
    y <- rnorm(1,xt,sigma)
    num <- f(y)*dnorm(xt,y,sigma)
    den <- f(xt)*dnorm(y,xt,sigma)
    if(u[i] <= num/den) x[i] <- y else{
      x[i] <- xt 
      k <- k+1 # y is rejected
    }
  }
  return(list(x=x,k=k))
}

# generate chains with four differeint sigma

N <- 15000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
mh1 <- mh(sigma[1],N,x0)
mh2 <- mh(sigma[2],N,x0)
mh3 <- mh(sigma[3],N,x0)
mh4 <- mh(sigma[4],N,x0)
# calculate the acceptance rates
print(c((N-mh1$k)/N,(N-mh2$k)/N,(N-mh3$k)/N,(N-mh4$k)/N))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i_th row of X
  pis <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means)
  psi.w <- apply(psi,1,"var")
  W <- mean(psi.w)
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W            # G-R statistic
  return(r.hat)
}

k <- 4
b <- 1000
sigma <- 0.5

# choose overdispersed initial values
x0 <- c(-10,-5,5,10)

# generate the chains
X <- matrix(0,nrow=k,ncol=N)
for(i in 1:k){
  X[i,] <- mh(sigma,N,x0[i])$x
}

# compute diagnostic statistics
psi <- t(apply(X,1,cumsum))
for(i in 1:nrow(psi))
  psi[i,] <- psi[i,]/(1:ncol(psi))
print(Gelman.Rubin(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0,N)
for(j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N],type="l",xlab="",ylab="R")
abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
# initialize constants and parameters
N <- 5000          # length of chain
burn <- 1000       # burn-in length
X <- matrix(0,N,2) # the chain, a bivariate sample

rho <- 0.9         # correlation
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

##### generate the chain #####

X[1, ] <- c(mu1,mu2)  # initialize
for(i in 2:N){
  x2 <- X[i-1,2]
  m1 <- mu1 + rho*(x2-mu2)*sigma1/sigma2
  X[i,1] <- rnorm(1,m1,s1)
  x1 <- X[i,1]
  m2 <- mu2 + rho*(x1-mu1)*sigma2/sigma1
  X[i,2] <- rnorm(1,m2,s2)
}

# Remove the first 1000 observations and draw the graph generated by the sample
b <- burn + 1
x <- X[b:N, ]

plot(x,main="",cex=.5,xlab=bquote(X[1]),ylab=bquote(X[2]),ylim=range(X[,2]))

lm <- lm(x[,1]~x[,2])
summary(lm)
qqnorm(x)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i_th row of X
  pis <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means)
  psi.w <- apply(psi,1,"var")
  W <- mean(psi.w)
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W            # G-R statistic
  return(r.hat)
}

k <- 4
b <- 1000
sigma <- 0.6

# choose overdispersed initial values
x0 <- c(-10,-5,5,10)

# generate the chains
X <- matrix(0,nrow=k,ncol=N)
for(i in 1:k){
  X[i,] <- mh(sigma,N,x0[i])$x
}

# compute diagnostic statistics
psi <- t(apply(X,1,cumsum))
for(i in 1:nrow(psi))
  psi[i,] <- psi[i,]/(1:ncol(psi))
print(Gelman.Rubin(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0,N)
for(j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N],type="l",xlab="",ylab="R")
abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
set.seed(22077)
root <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1); x2 <- rexp(N,1); x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(alpha+b1*x1+b2*x2+b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
    }
  solution <- uniroot(g,c(0,10))
  solution$root
}
N <- 1e6; b1 <- 0; b2 <- 1; b3  <- -1; f0 <-c(.1,.01,.001,.0001)
alpha <- numeric(4)
for(i in 1:4){
  alpha[i] <- root(N,b1,b2,b3,f0[i])
}
alpha
plot(f0,alpha,xlim=c(0,0.1),ylim=c(0,10))

## -----------------------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
x <- numeric(10)
obj <- function(lambda,u,v){
  for(i in 1: 10){
    x[i] <- (v[i]-u[i])/(exp(lambda*(v[i]-u[i])))-u[i]
  }
  sum(x)
}
u <- uniroot(obj,lower=-10,upper=10,u= u,v= v)
lambda.hat <- u$root
lambda.hat

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
beauty<-data.frame(x=0:4,y=6:10)
lapply(beauty,scale01)

## -----------------------------------------------------------------------------
beauty<-data.frame(x=c(0,0.5,0.7,0.9,0.5),y=c(6,7,8,9,10),z=c('a','b','c','d','e'))
a = c(0)
  for (i in 1:ncol(beauty)) {
    if (class(beauty[,i])!='numeric') a=c(a,i) 
  }
beauty1 = beauty[,-a]
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
lapply(beauty1,scale01)

## -----------------------------------------------------------------------------
beauty<-data.frame(x=c(0,0.5,0.7,0.9,0.5),y=c(6,7,8,9,10))
vapply(beauty,sd,c(1))

## -----------------------------------------------------------------------------
beauty<-data.frame(x=c(0,0.5,0.7,0.9,0.5),y=c(6,7,8,9,10),z=c('a','b','c','d','e'))
a <-vapply(beauty,is.numeric,c(1))
b <-c(0)
  for (i in 1:length(a)) {
    if (class(a[i])!='numeric') b=c(b,i) 
  }
beauty1 = beauty[,-b]
vapply(beauty1,sd,c(1))

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp <- '../Rcpp/'
# Can create source file in Rstudio
sourceCpp('D:/meanC.cpp')
library(microbenchmark)
x <- runif(1e4); mean2 <- function(x)sum(x)/length(x)
ts <- microbenchmark(meanR=mean(x),meanR2=mean2(x),
meancpp=meanC(x))
summary(ts)[,c(1,3,5,6)]


## -----------------------------------------------------------------------------
N <- 5000
burn <- 1000 # burn-in length
X <- matrix(0,N,2)

rho <-0.9 # correlation
mu1 <- mu2 <-0
sigma1 <- sigma2 <-1


# generate the chain
gibbs_r <- function(N,burn,X,rho,mu1,mu2,sigma1,sigma2){
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X[1,] <- c(mu1,mu2)
  for(i in 2:N){
    x2 <- X[i-1,2]
    m1 <- mu1+rho*(x2-mu2)*sigma1/sigma2
    X[i,1] <- rnorm(1,m1,s1)
    x1 <- X[i,1]
    m2 <- mu2+rho*(x1-mu1)*sigma2/sigma1
    X[i,2] <- rnorm(1,m2,s2)
  }
  
  b <- burn +1 
  x <- X[b:N,]
}

x1 <- gibbs_r(N,burn,X,rho,mu1,mu2,sigma1,sigma2)

## -----------------------------------------------------------------------------
N <- 5000
burn <- 1000 # burn-in length
X <- matrix(0,N,2)

rho <-0.9 # correlation
mu1 <- mu2 <-0
sigma1 <- sigma2 <-1
library(Rcpp)
sourceCpp('D:/gibbs_cpp.cpp')
x2 <- gibbs_cpp(N,mu1,mu2,sigma1,sigma2,rho)

## -----------------------------------------------------------------------------
qqplot(x1,x2)

## -----------------------------------------------------------------------------
library(microbenchmark)
microbenchmark(gibbs_r(N,burn,X,rho,mu1,mu2,sigma1,sigma2),gibbs_cpp(N,mu1,mu2,sigma1,sigma2,rho))

