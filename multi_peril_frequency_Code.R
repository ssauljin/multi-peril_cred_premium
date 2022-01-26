#### Load data files and creating training and test sets ####
library(MASS)
library(copula)
library(glmnet)

# Rawdata can be accessed from the following URL: "https://sites.google.com/a/wisc.edu/jed-frees"

load("data.RData")
load("dataout.RData")

train <-   data[,c(1:2,9:14,21:22,25,34,37,26,29)] # Only IM, CN, and PN claims are used
rm(data)
head(train)

train$exposureIM <- 1*(train$CoverageIM>0)
train$exposureCN <- 1*(train$CoverageCN>0)
train$exposurePN <- 1*(train$CoveragePN>0)

temp <- as.data.frame(train[,c(1,16:18)])

temp <- aggregate(temp[,2:4], by = list(temp$PolicyNum), FUN = sum)
colnames(temp)[1] <- "PolicyNum"

temp$loweIM     <- (temp$exposureIM < temp$exposureCN/3 + temp$exposurePN/3)*1
temp$loweCN     <- (temp$exposureCN < temp$exposureIM/3 + temp$exposurePN/3)*1
temp$lowePN     <- (temp$exposurePN < temp$exposureIM/3 + temp$exposureCN/3)*1

temp <- temp[,c(1,5:7)]

test <- dataout[,c(1:2,9:14,21:22,25,34,37,26,29)]
rm(dataout)
head(test)

test$exposureIM <- 1*(test$CoverageIM>0)
test$exposureCN <- 1*(test$CoverageCN>0)
test$exposurePN <- 1*(test$CoveragePN>0)

test <- merge(x = test, y = temp, by = "PolicyNum", all.x=TRUE)
test$loweIM[is.na(test$loweIM)] <- 0
test$loweCN[is.na(test$loweCN)] <- 0
test$lowePN[is.na(test$lowePN)] <- 0

rm(temp)

#### Preliminary analysis ####

# Measures of Association
attach(train)
cor.test(FreqCN, FreqIM, method="kendall") 
cor.test(FreqCN, FreqIM, method="spearman") 
cor.test(FreqCN, FreqIM, method="pearson") 

cor.test(FreqPN, FreqIM, method="kendall") 
cor.test(FreqPN, FreqIM, method="spearman") 
cor.test(FreqPN, FreqIM, method="pearson") 

cor.test(FreqPN, FreqCN, method="kendall") 
cor.test(FreqPN, FreqCN, method="spearman") 
cor.test(FreqPN, FreqCN, method="pearson") 


# Frequency table
FCN <- train$FreqCN
FIM <- train$FreqIM
FPN <- train$FreqPN

FCN[FCN>=5] <- 5
FIM[FIM>=5] <- 5
FPN[FPN>=5] <- 5

table(FCN, FIM)
table(FPN, FIM)
table(FCN, FPN)

glm.poisIM  <- glm(train$FreqIM~.,data=train[,c(3,5:10)  ],
                   family="poisson", weights=train$exposureIM)             # Poisson model for IM
glm.poisCN  <- glm(train$FreqCN~.,data=train[,c(3,5:8,12)],
                   family="poisson", weights=train$exposureCN)             # Poisson model for CN
glm.poisPN  <- glm(train$FreqPN~.,data=train[,c(3,5:8,14)],
                   family="poisson", weights=train$exposurePN)             # Poisson model for PN


elnet.cv.poisIM <- cv.glmnet(as.matrix(train[,c(3,5:10)  ]), train$FreqIM, type.measure="deviance", alpha=.5,
                             family="poisson", weights=train$exposureIM)

elnet.cv.poisCN <- cv.glmnet(as.matrix(train[,c(3,5:8,12)  ]), train$FreqCN, type.measure="deviance", alpha=.5,
                             family="poisson", weights=train$exposureCN)

elnet.cv.poisPN <- cv.glmnet(as.matrix(train[,c(3,5:8,14)  ]), train$FreqPN, type.measure="deviance", alpha=.5,
                             family="poisson", weights=train$exposurePN)

par(mfrow = c(1, 3), mar=c(5.1, 4.1, 5.1, 2.1))
plot(elnet.cv.poisIM, main="Variable selection - FreqIM")
plot(elnet.cv.poisCN, main="Variable selection - FreqCN")
plot(elnet.cv.poisPN, main="Variable selection - FreqPN")

Rprox <- as.data.frame(cbind(train$Year, train$PolicyNum, 
                             train$FreqCN, train$FreqIM, train$FreqPN, 
                             fitted(glm.poisCN)*train$exposureCN,
                             fitted(glm.poisIM)*train$exposureIM,
                             fitted(glm.poisPN)*train$exposurePN ))
colnames(Rprox) <- c("Year", "PolicyNum", "FreqCN", "FreqIM", "FreqPN", "nuCN", "nuIM", "nuPN")

dsumCN   <- 0
dsumIM   <- 0
dsumPN   <- 0
odsum    <- 0
odsum1   <- 0
odsum2   <- 0
odsum3   <- 0
nsumCN   <- 0
nsumIM   <- 0
nsumPN   <- 0
onsum    <- 0
onsum1   <- 0
onsum2   <- 0
onsum3   <- 0
n        <- 0
polnums  <- unique(Rprox$PolicyNum)
M <- length(polnums)

for (i in 1:M) {
  incn <- sum(Rprox$PolicyNum == polnums[i])
  vec  <- c(Rprox$FreqCN[(n+1):(n+incn)]-Rprox$nuCN[(n+1):(n+incn)],
            Rprox$FreqIM[(n+1):(n+incn)]-Rprox$nuIM[(n+1):(n+incn)],
            Rprox$FreqPN[(n+1):(n+incn)]-Rprox$nuPN[(n+1):(n+incn)])
  nec  <- c(Rprox$nuCN[(n+1):(n+incn)], Rprox$nuIM[(n+1):(n+incn)], Rprox$nuPN[(n+1):(n+incn)])
  
  vec1 <- Rprox$FreqIM[(n+1):(n+incn)]-Rprox$nuIM[(n+1):(n+incn)]
  vec2 <- Rprox$FreqCN[(n+1):(n+incn)]-Rprox$nuCN[(n+1):(n+incn)]
  vec3 <- Rprox$FreqPN[(n+1):(n+incn)]-Rprox$nuPN[(n+1):(n+incn)]
  
  nec1 <- Rprox$nuIM[(n+1):(n+incn)]
  nec2 <- Rprox$nuCN[(n+1):(n+incn)]
  nec3 <- Rprox$nuPN[(n+1):(n+incn)]
  
  mat     <- outer(vec , vec )
  ma1     <- outer(vec1, vec1)
  ma2     <- outer(vec2, vec2)
  ma3     <- outer(vec3, vec3)
  nat     <- outer(nec , nec )
  na1     <- outer(nec1, nec1)
  na2     <- outer(nec2, nec2)
  na3     <- outer(nec3, nec3)
  
  
  dsumCN  <- dsumCN  + sum(diag(mat)[        1 :   incn ])
  dsumIM  <- dsumIM  + sum(diag(mat)[(  incn+1):(2*incn)])
  dsumPN  <- dsumPN  + sum(diag(mat)[(2*incn+1):(3*incn)])
  odsum   <- odsum   - sum(diag(mat)) + sum(mat)
  odsum1  <- odsum1  - sum(diag(ma1)) + sum(ma1)
  odsum2  <- odsum2  - sum(diag(ma2)) + sum(ma2)
  odsum3  <- odsum3  - sum(diag(ma3)) + sum(ma3)
  
  nsumCN  <- nsumCN  + sum(diag(nat)[        1 :   incn ])
  nsumIM  <- nsumIM  + sum(diag(nat)[(  incn+1):(2*incn)])
  nsumPN  <- nsumPN  + sum(diag(nat)[(2*incn+1):(3*incn)])
  onsum   <- onsum   - sum(diag(nat)) + sum(nat)
  onsum1  <- onsum1  - sum(diag(na1)) + sum(na1)
  onsum2  <- onsum2  - sum(diag(na2)) + sum(na2)
  onsum3  <- onsum3  - sum(diag(na3)) + sum(na3)
  
  n <- n + incn }

r   <- onsum /odsum
r1  <- onsum1/odsum1
r2  <- onsum2/odsum2
r3  <- onsum3/odsum3
wIM <- sum(Rprox$nuIM)/(dsumIM -nsumIM/r)
wCN <- sum(Rprox$nuCN)/(dsumCN -nsumCN/r)
wPN <- sum(Rprox$nuPN)/(dsumPN -nsumPN/r)

# Test for overdispersion

phihat_IM <- sum(residuals(glm.poisIM, type="pearson")^2)/df.residual(glm.poisIM) # phi_hat of IM
phihat_CN <- sum(residuals(glm.poisCN, type="pearson")^2)/df.residual(glm.poisCN) # phi_hat of CN
phihat_PN <- sum(residuals(glm.poisPN, type="pearson")^2)/df.residual(glm.poisCN) # phi_hat of PN

resIM <- ((fitted(glm.poisIM)*train$exposureIM-train$FreqIM)^2 - train$FreqIM)/ fitted(glm.poisIM)/train$exposureIM
mean(resIM, na.rm=TRUE)                                         # phi_tilde of IM
mean(resIM, na.rm=TRUE)/sd(resIM, na.rm=TRUE)*sqrt(length(na.omit(resIM)))           # test statistic for phi_tilde of IM

resCN <- ((fitted(glm.poisCN)*train$exposureCN-train$FreqCN)^2 - train$FreqCN)/ fitted(glm.poisCN)/train$exposureCN
mean(resCN, na.rm=TRUE)                                         # phi_tilde of CN
mean(resCN, na.rm=TRUE)/sd(resCN, na.rm=TRUE)*sqrt(length(na.omit(resCN)))           # test statistic for phi_tilde of CN

resPN <- ((fitted(glm.poisPN)-train$FreqPN)^2 - train$FreqPN)/ fitted(glm.poisPN)
mean(resPN, na.rm=TRUE)                                         # phi_tilde of PN
mean(resPN, na.rm=TRUE)/sd(resPN, na.rm=TRUE)*sqrt(length(na.omit(resPN)))           # test statistic for phi_tilde of PN

library(AER) # Quicker way to test the significance of phi_tilde
  dispersiontest(glm.poisIM)  
  dispersiontest(glm.poisCN)  
  dispersiontest(glm.poisPN)  

#### Model calibration ####

glm.nbIM <- glm.nb(train$FreqIM~.,data=train[,c(3,5:10)  ],
                   control=glm.control(maxit=50), weights=train$exposureIM)   # NB model for IM
glm.nbCN <- glm.nb(train$FreqCN~.,data=train[,c(3,5:8,12)],
                   control=glm.control(maxit=50), weights=train$exposureCN)   # NB model for CN
glm.nbPN <- glm.nb(train$FreqPN~.,data=train[,c(3,5:8,14)],
                   control=glm.control(maxit=50), weights=train$exposurePN)   # NB model for PN

require(pscl)
zip.IM <- zeroinfl(train$FreqIM~.,data=train[,c(3,5:10)  ],
                   dist = "poisson", weights=train$exposureIM)             # ZIP model for IM
zip.CN <- zeroinfl(train$FreqCN~.,data=train[,c(3,5:8,12)],
                   dist = "poisson", weights=train$exposureCN)             # ZIP model for CN
zip.PN <- zeroinfl(train$FreqPN~.,data=train[,c(3,5:8,14)],
                   dist = "poisson", weights=train$exposurePN)             # ZIP model for PN

# optimization of Shi and Valdez (2014) for multivariate frequency with Frank copula
copEst <- function(x1, x2, x3, n1, n2, n3, init.alpha1, init.alpha2, init.alpha3, init.phi) { 
  x1 <- cbind(rep(1,nrow(x1)),x1)
  colnames(x1)[1] <- "intercept"
  x2 <- cbind(rep(1,nrow(x2)),x2)
  colnames(x2)[1] <- "intercept"
  x3 <- cbind(rep(1,nrow(x3)),x3)
  colnames(x3)[1] <- "intercept"
  
  
  # marginal likelihood of multivariate NB distribution
  "negll.frankcop" <- function(parm) {
    e1 <- ncol(x1);
    e2 <- ncol(x2);
    e3 <- ncol(x3);
  
    cop <- frankCopula(parm[e1+e2+e3+1], dim=3)
    lambdaIM <- exp(as.matrix(x1) %*% parm[       1 :       e1 ] )*train$exposureIM
    lambdaCN <- exp(as.matrix(x2) %*% parm[(   e1+1):(   e1+e2)] )*train$exposureCN
    lambdaPN <- exp(as.matrix(x3) %*% parm[(e1+e2+1):(e1+e2+e3)] )*train$exposurePN
    
    c  <- ppois(cbind(n1  , n2   , n3  ), cbind(lambdaIM,lambdaCN,lambdaPN))
    c1 <- ppois(cbind(n1-1, n2-1 , n3-1), cbind(lambdaIM,lambdaCN,lambdaPN))
    likfrank <- pCopula(cbind( c[,1], c[,2], c[,3]),cop)-pCopula(cbind( c[,1], c[,2],c1[,3]),cop)+
                pCopula(cbind(c1[,1],c1[,2], c[,3]),cop)-pCopula(cbind(c1[,1],c1[,2],c1[,3]),cop)+
                pCopula(cbind( c[,1],c1[,2],c1[,3]),cop)-pCopula(cbind( c[,1],c1[,2], c[,3]),cop)+
                pCopula(cbind(c1[,1], c[,2],c1[,3]),cop)-pCopula(cbind(c1[,1], c[,2], c[,3]),cop)
    return(-sum(log(pmax(likfrank,10^-10))))
  } 
  
  init.est <- as.vector(c(init.alpha1,init.alpha2, init.alpha3, init.phi))
  # init.est <- as.vector(c(init.alpha,init.w))
  
  fit.frankcop <- optim(init.est, negll.frankcop, method="BFGS", hessian=TRUE)
  parm.hat <- fit.frankcop$par
  loglik.frankcop <- -fit.frankcop$value
  inv.frankcop.Hess <- solve(fit.frankcop$hessian);
  parm.se <- sqrt(diag(inv.frankcop.Hess));
  # put together the model with the est, se, t, pval, AIC, BIC
  dfe <- length(n1+n2+n3-length(parm.hat));
  t_ratio<-parm.hat/parm.se;
  #test if diff. from 1 t_ratio[1:3]<-(parm.hat[1:3]-1)/parm.se[1:3];
  pval <- pf(t_ratio*t_ratio,df1=1,df2=dfe,lower.tail=F);
  ttable <- cbind(parm.hat,parm.se,t_ratio,pval) 
  ttable <- round(ttable,digits=4)
  
  rownames(ttable)<- c(colnames(x1),colnames(x2),colnames(x3),"phi")
  # rownames(ttable)<- c(colnames(x),"weight")
  colnames(ttable)<- c("estimate", "std error", "t-val","Pr>|t|");
  
  #  AIC<- 2*negll.frankcop(parm.hat) + 2*length(parm.hat);
  #  BIC<- 2*negll.frankcop(parm.hat) + log(length(n))*length(parm.hat);
  #  loglik <- negll.frankcop(parm.hat)
  return(list(ttable=ttable, coef=parm.hat));
}

system.time(
  copM <- copEst(x1=train[,c(3,5:10)  ] ,n1=train$FreqIM, init.alpha1=coefficients(glm.poisIM), 
                 x2=train[,c(3,5:8,12)] ,n2=train$FreqCN, init.alpha2=coefficients(glm.poisCN),
                 x3=train[,c(3,5:8,14)] ,n3=train$FreqPN, init.alpha3=coefficients(glm.poisPN),
                 init.phi=1.0)) 


# Estimation of alpha via marginal frequency likelihood for proposed model (Multi-Peril)
NBEst2 <- function(x1, x2, n1, n2, id, init.alpha1, init.alpha2, w1=1, w2=1, r) { # optimization of proposed marginal likelihood for frequency
  x1 <- cbind(rep(1,nrow(x1)),x1)
  colnames(x1)[1] <- "intercept"
  x2 <- cbind(rep(1,nrow(x2)),x2)
  colnames(x2)[1] <- "intercept"
  
  # marginal likelihood of multivariate NB distribution
  "negll.NB" <- function(parm) {
    e1 <- ncol(x1);
    e2 <- ncol(x2);
    reg_eqn1 <- as.matrix(x1) %*% parm[1:e1]
    reg_eqn2 <- as.matrix(x2) %*% parm[(e1+1):(e1+e2)]
    data <- cbind(id,w1*exp(reg_eqn1)*train$exposureIM,w1*n1,
                     w2*exp(reg_eqn2)*train$exposureCN,w2*n2);
    #   data <- cbind(id,exp(reg_eqn),n);
    
    colnames(data)[2] <- "sv1";
    colnames(data)[3] <- "n1";
    colnames(data)[4] <- "sv2";
    colnames(data)[5] <- "n2";
    
    temp1 = sum(w1*n1*reg_eqn1+w2*n2*reg_eqn2-log(gamma(w1*n1+w2*n2+1)))+length(unique(id))*(r*log(r)-log(gamma(r)));
    temp2 = -sum((as.matrix(aggregate(n1~id,data,sum))[,2]+as.matrix(aggregate(n2~id,data,sum))[,2]+r)*log(
      as.matrix(aggregate(sv1~id,data,sum))[,2]+as.matrix(aggregate(sv2~id,data,sum))[,2]+r)
    )+sum(log(gamma(as.matrix(aggregate(n1~id,data,sum))[,2]+as.matrix(aggregate(n2~id,data,sum))[,2]+r)));
    result = -temp1-temp2
    return(result)
  } 
  init.est <- as.vector(c(init.alpha1,init.alpha2))
  # init.est <- as.vector(c(init.alpha,init.w))
  
  fit.NB <- optim(init.est, negll.NB, NULL)
  parm.hat <- fit.NB$par
  loglik.NB <- -fit.NB$value
  
  # next estimate the standard errors.
  library(nlme)
  negll.NB.Hess <- fdHess(parm.hat, negll.NB);
  inv.NB.Hess <- solve(negll.NB.Hess$Hessian);
  parm.se <- sqrt(diag(inv.NB.Hess));
  # put together the model with the est, se, t, pval, AIC, BIC
  dfe <- length(n1+n2-length(parm.hat));
  t_ratio<-parm.hat/parm.se;
  #test if diff. from 1 t_ratio[1:3]<-(parm.hat[1:3]-1)/parm.se[1:3];
  pval <- pf(t_ratio*t_ratio,df1=1,df2=dfe,lower.tail=F);
  ttable <- cbind(parm.hat,parm.se,t_ratio,pval) 
  ttable <- round(ttable,digits=4)
  
  rownames(ttable)<- c(colnames(x1),colnames(x2))
  # rownames(ttable)<- c(colnames(x),"weight")
  colnames(ttable)<- c("estimate", "std error", "t-val","Pr>|t|");
  
  #  AIC<- 2*negll.NB(parm.hat) + 2*length(parm.hat);
  #  BIC<- 2*negll.NB(parm.hat) + log(length(n))*length(parm.hat);
  #  loglik <- negll.NB(parm.hat)
  return(list(ttable=ttable, coef=parm.hat));
}

set.seed(100)
system.time(
  NBm_Joint <- NBEst2(x1=train[,c(3,5:10)     ] ,n1=train$FreqIM, init.alpha1=coefficients(glm.poisIM), w1=wIM,
                      x2=train[,c(3,5:8,12)] ,n2=train$FreqCN, init.alpha2=coefficients(glm.poisCN), w2=wCN,
                      id=train$PolicyNum, r=r)) 

# Estimation of alpha via marginal frequency likelihood for proposed model (Multi-Peril)
NBEst3 <- function(x1, x2, x3, n1, n2, n3, id, init.alpha1, init.alpha2, init.alpha3, w1=1, w2=1, w3=1, r) { # optimization of proposed marginal likelihood for frequency
  x1 <- cbind(rep(1,nrow(x1)),x1)
  colnames(x1)[1] <- "intercept"
  x2 <- cbind(rep(1,nrow(x2)),x2)
  colnames(x2)[1] <- "intercept"
  x3 <- cbind(rep(1,nrow(x3)),x3)
  colnames(x3)[1] <- "intercept"
  
  
  # marginal likelihood of multivariate NB distribution
  "negll.NB" <- function(parm) {
    e1 <- ncol(x1);e2 <- ncol(x2);e3 <- ncol(x3);

    reg_eqn1 <- as.matrix(x1) %*% parm[1:e1]
    reg_eqn2 <- as.matrix(x2) %*% parm[(e1+1):(e1+e2)]
    reg_eqn3 <- as.matrix(x3) %*% parm[(e1+e2+1) :(e1+e2+e3) ]
    
    data <- cbind(id, w1*exp(reg_eqn1)*train$exposureIM, w1*n1, 
                      w2*exp(reg_eqn2)*train$exposureCN, w2*n2,
                      w3*exp(reg_eqn3)*train$exposurePN, w3*n3);
    
    colnames(data)[2] <- "sv1";
    colnames(data)[3] <- "n1";
    colnames(data)[4] <- "sv2";
    colnames(data)[5] <- "n2";
    colnames(data)[6] <- "sv3";
    colnames(data)[7] <- "n3";
    
    temp1 = sum(w1*n1*reg_eqn1+w2*n2*reg_eqn2+w3*n3*reg_eqn3-log(gamma(w1*n1+w2*n2+w3*n3+1)))+
      length(unique(id))*(r*log(r)-log(gamma(r)));
    temp2 =   -sum((as.matrix(aggregate(n1 ~id,data,sum))[,2]+
                    as.matrix(aggregate(n2 ~id,data,sum))[,2]+ 
                    as.matrix(aggregate(n3 ~id,data,sum))[,2]+ r)*log(
                    as.matrix(aggregate(sv1~id,data,sum))[,2]+
                    as.matrix(aggregate(sv2~id,data,sum))[,2]+ 
                    as.matrix(aggregate(sv3~id,data,sum))[,2]+ r))+ sum( 
          log(gamma(as.matrix(aggregate(n1 ~id,data,sum))[,2]+
                    as.matrix(aggregate(n2 ~id,data,sum))[,2]+
                    as.matrix(aggregate(n3 ~id,data,sum))[,2]+ r)));
    result = -temp1-temp2
    return(result)
  } 
  init.est <- as.vector(c(init.alpha1, init.alpha2, init.alpha3))
  #  init.est <- as.vector(c(init.alpha1,init.alpha2, init.alpha3))
  
  fit.NB <- optim(init.est, negll.NB, NULL)
  parm.hat <- fit.NB$par
  loglik.NB <- -fit.NB$value
  
  # next estimate the standard errors.
  library(nlme)
  negll.NB.Hess <- fdHess(parm.hat, negll.NB);
  inv.NB.Hess <- solve(negll.NB.Hess$Hessian);
  parm.se <- sqrt(diag(inv.NB.Hess));
  # put together the model with the est, se, t, pval, AIC, BIC
  dfe <- length(n3+n1+n2-length(parm.hat));
  t_ratio<-parm.hat/parm.se;
  #test if diff. from 1 t_ratio[1:3]<-(parm.hat[1:3]-1)/parm.se[1:3];
  pval <- pf(t_ratio*t_ratio,df1=1,df2=dfe,lower.tail=F);
  ttable <- cbind(parm.hat,parm.se,t_ratio,pval) 
  ttable <- round(ttable,digits=4)
  
  rownames(ttable)<- c(colnames(x1), colnames(x2), colnames(x3))
  # rownames(ttable)<- c(colnames(x),"weight")
  colnames(ttable)<- c("estimate", "std error", "t-val","Pr>|t|");
  
  #  AIC<- 2*negll.NB(parm.hat) + 2*length(parm.hat);
  #  BIC<- 2*negll.NB(parm.hat) + log(length(n))*length(parm.hat);
  #  loglik <- negll.NB(parm.hat)
  return(list(ttable=ttable, coef=parm.hat));
}

set.seed(99)
U <- runif(7)/10-0.05

NBm_Jointt <- NBEst3(x1=train[,c(3,5:10)  ] ,n1=train$FreqIM, init.alpha1=NBm_Joint$coef[1 :8]      , w1=wIM,
                     x2=train[,c(3,5:8,12)] ,n2=train$FreqCN, init.alpha2=NBm_Joint$coef[9:15]      , w2=wCN,
                     x3=train[,c(3,5:8,14)] ,n3=train$FreqPN, init.alpha3=coefficients(glm.poisPN)+U, w3=wPN,
                     id=train$PolicyNum, r=r) 

options(scipen = 999)
GLMtable <- function(object) {
  coef.beta <- coef(object)
  vc <- object$vcov
  if (is.null(vc)) {vc <- vcov(object)}
  s.err <- sqrt(diag(vc))    
  err.beta <- s.err
  test.value <- coef.beta / err.beta
  dn <- c("Estimate", "s.e.")             
  pvalue <- round(2 * pt(-abs(test.value), object$df.residual),9)
  coef.table <- cbind(coef.beta, err.beta, pvalue)  
  dn2 <- "Pr(>|t|)"
  dimnames(coef.table) <- list(names(coef.beta), c(dn, dn2))
  return(coef.table) }

library(knitr)
library(kableExtra)
options(knitr.table.format = "latex")



# Summary of regression coefficients for IM claim from four models
IMtable <- round(cbind(GLMtable(glm.poisIM)[,-2], GLMtable(zip.IM)[1:8,-2], copM$ttable[1:8,c(1,4)], 
                       NBm_Jointt$ttable[1:8,c(1,4)]),digits=4)
colnames(IMtable) <- c("Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value")

kable(IMtable, booktabs = T, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{=========}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")

# Summary of regression coefficients for CN claim from four models
CNtable <- round(cbind(GLMtable(glm.poisCN)[,-2], GLMtable(zip.CN)[1:7,-2], copM$ttable[9:15,c(1,4)]+cbind(rep(0,7), c(rep(0,6),0.0001)),
                       NBm_Jointt$ttable[9:15,c(1,4)]),digits=4)
colnames(CNtable) <- c("Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value")

kable(CNtable, booktabs = T, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{=========}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")

# Summary of regression coefficients for CN claim from four models
PNtable <- round(cbind(GLMtable(glm.poisPN)[,-2], GLMtable(zip.PN)[1:7,-2], copM$ttable[16:22,c(1,4)]+cbind(rep(0,7), c(rep(0,6),0.0001)),
                       NBm_Jointt$ttable[16:22,c(1,4)]),digits=4)
colnames(PNtable) <- c("Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value")

kable(PNtable, booktabs = T, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{=========}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")



#### Out-of-sample validation ####
id <- train$PolicyNum
x1 <- train[,c(3,5:10)     ]
x1 <- cbind(rep(1,nrow(x1)),x1)
n1 <- train$FreqIM
#w1 <- 1/phihat_IM
w1 <- wIM  

x2 <- train[,c(3,5:8,12)]
x2 <- cbind(rep(1,nrow(x2)),x2)
n2 <- train$FreqCN
#w2 <- 1/phihat_CN
w2 <- wCN

x3 <- train[,c(3,5:8,14)]
x3 <- cbind(rep(1,nrow(x3)),x3)
n3 <- train$FreqPN
#w3 <- 1/phihat_CN
w3 <- wPN


# calculation of bonus-malus factor based on the proposed model for frequency
Nreg_eqn1 <- as.matrix(x1) %*% as.matrix(NBm_Jointt$coef[1 :8 ])
Nreg_eqn2 <- as.matrix(x2) %*% as.matrix(NBm_Jointt$coef[9 :15])
Nreg_eqn3 <- as.matrix(x3) %*% as.matrix(NBm_Jointt$coef[16:22])
Ndata <- cbind(id,exp(Nreg_eqn1)*train$exposureIM,n1,
                  exp(Nreg_eqn2)*train$exposureCN,n2,
                  exp(Nreg_eqn3)*train$exposurePN,n3);
colnames(Ndata)[2] <- "nv1";
colnames(Ndata)[4] <- "nv2";
colnames(Ndata)[6] <- "nv3";

Npost     <- aggregate(n1 ~id,Ndata,sum)      # aggregate actual   IM claim counts for years
Npost$nv1 <- aggregate(nv1~id,Ndata,sum)[,2]  # aggregate expected IM claim counts for years
Npost$n2  <- aggregate(n2 ~id,Ndata,sum)[,2]  # aggregate actual   CN claim counts for years
Npost$nv2 <- aggregate(nv2~id,Ndata,sum)[,2]  # aggregate expected CN claim counts for years
Npost$n3  <- aggregate(n3 ~id,Ndata,sum)[,2]  # aggregate actual   PN claim counts for years
Npost$nv3 <- aggregate(nv3~id,Ndata,sum)[,2]  # aggregate expected PN claim counts for years

Npost$nweight <- (w1*Npost$n1 +w2*Npost$n2  +w3*Npost$n3+ r) / (w1*Npost$nv1 + w2*Npost$nv2 + w3*Npost$nv3+ r)
Npost$nweigh1 <- (   Npost$n1                           +r1) / (   Npost$nv1                              +r1)
Npost$nweigh2 <- (                Npost$n2              +r2) / (                  Npost$nv2               +r2)
Npost$nweigh3 <- (                              Npost$n3+r3) / (                                 Npost$nv3+r3)

# bonus-malus factor for a policyholder on frequency
colnames(Npost)[1] <- "PolicyNum"

# attach the bonus-malus factor for each policyholder on the test set
Ptest <- merge(x = test, y = Npost, by = "PolicyNum", all.x = TRUE)
Ptest$nweight[is.na(Ptest$nweight)] <- 1
Ptest$nweigh1[is.na(Ptest$nweigh1)] <- 1
Ptest$nweigh2[is.na(Ptest$nweigh2)] <- 1
Ptest$nweigh3[is.na(Ptest$nweigh3)] <- 1


xt1 <- test[,c(3,5:10)     ]
xt1 <- cbind(rep(1,nrow(xt1)),xt1)
xt2 <- test[,c(3,5:8,12)]
xt2 <- cbind(rep(1,nrow(xt2)),xt2)
xt3 <- test[,c(3,5:8,14)]
xt3 <- cbind(rep(1,nrow(xt3)),xt3)

n1_poispred <- exp(as.matrix(xt1) %*% coefficients(glm.poisIM)) *test$exposureIM # IM frequency premium with Poisson model
n1_coppred <- exp(as.matrix(xt1) %*% as.matrix(copM$coef[1:8])) *test$exposureIM # IM frequency premium with copula model
n1_zippred   <- predict(zip.IM, xt1) *test$exposureIM                            # IM frequency premium with ZIP model
n1_propred <- exp(as.matrix(xt1) %*% as.matrix(NBm_Jointt$coef[1:8])
                  )*Ptest$nweight *test$exposureIM                               # IM frequency premium with Proposed model
n1_unipred <- exp(as.matrix(xt1) %*% coefficients(glm.poisIM) 
                  )*Ptest$nweigh1 *test$exposureIM                               # IM frequency premium with Uni-cred model

n2_poispred <- exp(as.matrix(xt2) %*% coefficients(glm.poisCN)) *test$exposureCN # CN frequency premium with Poisson model
n2_coppred <- exp(as.matrix(xt2) %*% as.matrix(copM$coef[9:15]))*test$exposureCN    # CN frequency premium with copula model
n2_zippred   <- predict(zip.CN, xt2) *test$exposureCN                           # CN frequency premium with ZIP model
n2_propred <- exp(as.matrix(xt2) %*% as.matrix(NBm_Jointt$coef[9:15])
                  )*Ptest$nweight *test$exposureCN                              # CN frequency premium with Proposed model
n2_unipred <- exp(as.matrix(xt2) %*% coefficients(glm.poisCN)
                  )*Ptest$nweigh2 *test$exposureCN                              # CN frequency premium with Uni-cred model

n3_poispred <- exp(as.matrix(xt3) %*% coefficients(glm.poisPN))  *test$exposurePN # PN frequency premium with Poisson model
n3_coppred <- exp(as.matrix(xt3) %*% as.matrix(copM$coef[16:22]))*test$exposurePN    # PN frequency premium with copula model
n3_zippred   <- predict(zip.PN, xt3) *test$exposurePN                           # PN frequency premium with ZIP model
n3_propred <- exp(as.matrix(xt3) %*% as.matrix(NBm_Jointt$coef[16:22])
                  )*Ptest$nweight    *test$exposurePN                           # PN frequency premium with Proposed model
n3_unipred <- exp(as.matrix(xt3) %*% coefficients(glm.poisPN)
                  )*Ptest$nweigh3    *test$exposurePN                           # PN frequency premium with Uni-cred model

# n_poispred  <- c(n1_poispred, n2_poispred, n3_poispred)
# n_coppred   <- c(n1_coppred , n2_coppred , n3_coppred )
# n_zippred   <- c(n1_zippred , n2_zippred , n3_zippred )
# n_propred   <- c(n1_propred , n2_propred , n3_propred )
# n_actual    <- c(test$FreqIM, test$FreqCN, test$FreqPN)
n1_actual   <- test$FreqIM
n2_actual   <- test$FreqCN
n3_actual   <- test$FreqPN

# RMSE (root mean squared errors) of four models
# RMSE_pois      <- sqrt(mean((n_poispred - n_actual)^2))
# RMSE_cop       <- sqrt(mean((n_coppred  - n_actual)^2))
# RMSE_zip       <- sqrt(mean((n_zippred  - n_actual)^2))
# RMSE_proposed  <- sqrt(mean((n_propred  - n_actual)^2))
# 
# RMSE_pois
# RMSE_cop
# RMSE_zip
# RMSE_proposed

RMSE_pois1     <- sqrt(mean((n1_poispred - n1_actual)^2))
RMSE_cop1      <- sqrt(mean((n1_coppred  - n1_actual)^2))
RMSE_zip1      <- sqrt(mean((n1_zippred  - n1_actual)^2))
RMSE_proposed1 <- sqrt(mean((n1_propred  - n1_actual)^2))
RMSE_unicred1  <- sqrt(mean((n1_unipred  - n1_actual)^2))

RMSE_pois1
RMSE_cop1
RMSE_zip1
RMSE_proposed1
RMSE_unicred1

RMSE_pois2     <- sqrt(mean((n2_poispred - n2_actual)^2))
RMSE_cop2      <- sqrt(mean((n2_coppred  - n2_actual)^2))
RMSE_zip2      <- sqrt(mean((n2_zippred  - n2_actual)^2))
RMSE_proposed2 <- sqrt(mean((n2_propred  - n2_actual)^2))
RMSE_unicred2  <- sqrt(mean((n2_unipred  - n2_actual)^2))

RMSE_pois2
RMSE_cop2
RMSE_zip2
RMSE_proposed2
RMSE_unicred2

RMSE_pois3     <- sqrt(mean((n3_poispred - n3_actual)^2))
RMSE_cop3      <- sqrt(mean((n3_coppred  - n3_actual)^2))
RMSE_zip3      <- sqrt(mean((n3_zippred  - n3_actual)^2))
RMSE_proposed3 <- sqrt(mean((n3_propred  - n3_actual)^2))
RMSE_unicred3  <- sqrt(mean((n3_unipred  - n3_actual)^2))

RMSE_pois3
RMSE_cop3
RMSE_zip3
RMSE_proposed3
RMSE_unicred3


# MAE (mean absoulte errors) for four models
# MAE_pois    <- mean(abs(n_poispred - n_actual))
# MAE_cop     <- mean(abs(n_coppred - n_actual))
# MAE_zip      <- mean(abs(n_zippred - n_actual))
# MAE_proposed <- mean(abs(n_propred - n_actual))
# 
# MAE_pois
# MAE_cop
# MAE_zip
# MAE_proposed

MAE_pois1     <- mean(abs(n1_poispred - n1_actual))
MAE_cop1      <- mean(abs(n1_coppred - n1_actual))
MAE_zip1      <- mean(abs(n1_zippred - n1_actual))
MAE_proposed1 <- mean(abs(n1_propred - n1_actual))
MAE_unicred1  <- mean(abs(n1_unipred - n1_actual))


MAE_pois1
MAE_cop1
MAE_zip1
MAE_proposed1
MAE_unicred1 

MAE_pois2     <- mean(abs(n2_poispred - n2_actual))
MAE_cop2      <- mean(abs(n2_coppred - n2_actual))
MAE_zip2      <- mean(abs(n2_zippred - n2_actual))
MAE_proposed2 <- mean(abs(n2_propred - n2_actual))
MAE_unicred2  <- mean(abs(n2_unipred - n2_actual))

MAE_pois2
MAE_cop2
MAE_zip2
MAE_proposed2
MAE_unicred2

MAE_pois3     <- mean(abs(n3_poispred - n3_actual))
MAE_cop3      <- mean(abs(n3_coppred - n3_actual))
MAE_zip3      <- mean(abs(n3_zippred - n3_actual))
MAE_proposed3 <- mean(abs(n3_propred - n3_actual))
MAE_unicred3  <- mean(abs(n3_unipred - n3_actual))

MAE_pois3
MAE_cop3
MAE_zip3
MAE_proposed3
MAE_unicred3

dev <- function(pred, actual) {
  pred   <- pred   + 1e-23
  actual <- actual + 1e-23
  result <- 2*sum(actual*log(actual/pred)-actual+pred)
  return(result) }

DEV_pois1     <- dev(n1_poispred,n1_actual)
DEV_cop1      <- dev(n1_coppred, n1_actual)
DEV_zip1      <- dev(n1_zippred, n1_actual)
DEV_proposed1 <- dev(n1_propred, n1_actual)
DEV_unicred1  <- dev(n1_unipred, n1_actual)

DEV_pois1
DEV_cop1
DEV_zip1
DEV_proposed1
DEV_unicred1

DEV_pois2     <- dev(n2_poispred,n2_actual)
DEV_cop2      <- dev(n2_coppred, n2_actual)
DEV_zip2      <- dev(n2_zippred, n2_actual)
DEV_proposed2 <- dev(n2_propred, n2_actual)
DEV_unicred2  <- dev(n2_unipred, n2_actual)

DEV_pois2
DEV_cop2
DEV_zip2
DEV_proposed2
DEV_unicred2

DEV_pois3     <- dev(n3_poispred,n3_actual)
DEV_cop3      <- dev(n3_coppred, n3_actual)
DEV_zip3      <- dev(n3_zippred, n3_actual)
DEV_proposed3 <- dev(n3_propred, n3_actual)
DEV_unicred3  <- dev(n3_unipred, n3_actual)

DEV_pois3
DEV_cop3
DEV_zip3
DEV_proposed3
DEV_unicred3
# 
# RMSE_pois1l     <- sqrt(sum(test$loweIM*(n1_poispred - n1_actual)^2)/sum(test$loweIM))
# RMSE_cop1l      <- sqrt(sum(test$loweIM*(n1_coppred  - n1_actual)^2)/sum(test$loweIM))
# RMSE_zip1l      <- sqrt(sum(test$loweIM*(n1_zippred  - n1_actual)^2)/sum(test$loweIM))
# RMSE_proposed1l <- sqrt(sum(test$loweIM*(n1_propred  - n1_actual)^2)/sum(test$loweIM))
# RMSE_unicred1l  <- sqrt(sum(test$loweIM*(n1_unipred  - n1_actual)^2)/sum(test$loweIM))
# 
# RMSE_pois1l
# RMSE_cop1l
# RMSE_zip1l
# RMSE_proposed1l
# RMSE_unicred1l
# 
# RMSE_pois2l     <- sqrt(sum(test$loweCN*(n2_poispred - n2_actual)^2)/sum(test$loweCN))
# RMSE_cop2l      <- sqrt(sum(test$loweCN*(n2_coppred  - n2_actual)^2)/sum(test$loweCN))
# RMSE_zip2l      <- sqrt(sum(test$loweCN*(n2_zippred  - n2_actual)^2)/sum(test$loweCN))
# RMSE_proposed2l <- sqrt(sum(test$loweCN*(n2_propred  - n2_actual)^2)/sum(test$loweCN))
# RMSE_unicred2l  <- sqrt(sum(test$loweCN*(n2_unipred  - n2_actual)^2)/sum(test$loweCN))
# 
# RMSE_pois2l
# RMSE_cop2l
# RMSE_zip2l
# RMSE_proposed2l
# RMSE_unicred2l
# 
# RMSE_pois3l     <- sqrt(sum(test$lowePN*(n3_poispred - n3_actual)^2)/sum(test$lowePN))
# RMSE_cop3l      <- sqrt(sum(test$lowePN*(n3_coppred  - n3_actual)^2)/sum(test$lowePN))
# RMSE_zip3l      <- sqrt(sum(test$lowePN*(n3_zippred  - n3_actual)^2)/sum(test$lowePN))
# RMSE_proposed3l <- sqrt(sum(test$lowePN*(n3_propred  - n3_actual)^2)/sum(test$lowePN))
# RMSE_unicred3l  <- sqrt(sum(test$lowePN*(n3_unipred  - n3_actual)^2)/sum(test$lowePN))
# 
# RMSE_pois3l
# RMSE_cop3l
# RMSE_zip3l
# RMSE_proposed3l
# RMSE_unicred3l
# 
# MAE_pois1l     <- mean(test$loweIM*abs(n1_poispred - n1_actual))
# MAE_cop1l      <- mean(test$loweIM*abs(n1_coppred - n1_actual))
# MAE_zip1l      <- mean(test$loweIM*abs(n1_zippred - n1_actual))
# MAE_proposed1l <- mean(test$loweIM*abs(n1_propred - n1_actual))
# MAE_unicred1l  <- mean(test$loweIM*abs(n1_unipred - n1_actual))
# 
# MAE_pois1l
# MAE_cop1l
# MAE_zip1l
# MAE_proposed1l
# MAE_unicred1l
# 
# MAE_pois2l     <- mean(test$loweCN*abs(n2_poispred - n2_actual))
# MAE_cop2l      <- mean(test$loweCN*abs(n2_coppred - n2_actual))
# MAE_zip2l      <- mean(test$loweCN*abs(n2_zippred - n2_actual))
# MAE_proposed2l <- mean(test$loweCN*abs(n2_propred - n2_actual))
# MAE_unicred2l  <- mean(test$loweCN*abs(n2_unipred - n2_actual))
# 
# MAE_pois2l
# MAE_cop2l
# MAE_zip2l
# MAE_proposed2l
# MAE_unicred2l
# 
# 
# MAE_pois3l     <- mean(test$lowePN*abs(n3_poispred - n3_actual))
# MAE_cop3l      <- mean(test$lowePN*abs(n3_coppred - n3_actual))
# MAE_zip3l      <- mean(test$lowePN*abs(n3_zippred - n3_actual))
# MAE_proposed3l <- mean(test$lowePN*abs(n3_propred - n3_actual))
# MAE_unicred3l  <- mean(test$lowePN*abs(n3_unipred - n3_actual))
# 
# MAE_pois3l
# MAE_cop3l
# MAE_zip3l
# MAE_proposed3l
# MAE_unicred3l


predtable1 <- rbind(
c(RMSE_pois1, RMSE_unicred1, RMSE_cop1, RMSE_zip1, RMSE_proposed1),
c(MAE_pois1, MAE_unicred1, MAE_cop1, MAE_zip1, MAE_proposed1),
c(DEV_pois1, DEV_unicred1, DEV_cop1, DEV_zip1, DEV_proposed1))
rownames(predtable1) <- c("RMSE", "MAE", "DEV")

kable(predtable1, booktabs = T, digits=4, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{======}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")


predtable2 <- rbind(
  c(RMSE_pois2, RMSE_unicred2, RMSE_cop2, RMSE_zip2, RMSE_proposed2),
  c(MAE_pois2, MAE_unicred2, MAE_cop2, MAE_zip2, MAE_proposed2),
  c(DEV_pois2, DEV_unicred2, DEV_cop2, DEV_zip2, DEV_proposed2))
rownames(predtable2) <- c("RMSE", "MAE", "DEV")

kable(predtable2, booktabs = T, digits=4, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{======}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")



predtable3 <- rbind(
  c(RMSE_pois3, RMSE_unicred3, RMSE_cop3, RMSE_zip3, RMSE_proposed3),
  c(MAE_pois3, MAE_unicred3, MAE_cop3, MAE_zip3, MAE_proposed3),
  c(DEV_pois3, DEV_unicred3, DEV_cop3, DEV_zip3, DEV_proposed3))
rownames(predtable3) <- c("RMSE", "MAE", "DEV")

kable(predtable3, booktabs = T, digits=4, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{======}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")


Npost$n <- Npost$n1 + Npost$n2 + Npost$n3
Npost$n[Npost$n> 2 ]<- "> 2"
Npost$n[Npost$n==0] <- "= 0"
Npost$n[Npost$n==1] <- "= 1"
Npost$n[Npost$n==2] <- "= 2"
#Npost$n[Npost$n> 2] <- "> 2"

Npost$n <- factor(Npost$n, levels = c("= 0","= 1", "= 2", "> 2"))
library(ggplot2)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



ggplot(Npost, aes(x=n, y=nweight, color=n)) +  geom_boxplot() +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  labs(x="Number of claims observed", y="Credibility factors",
       title = "Credibility factors for frequency\nunder the proposed model") +
  theme(legend.position = c(.10, .82), plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="#0072B2") +
  scale_y_continuous(labels = scales::percent) + guides(fill=guide_legend(ncol=2))

#, axis.text.x=element_blank()
