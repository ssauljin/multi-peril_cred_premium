#### Load data files and creating training and test sets ####
library(MASS)
library(copula)

# load rawdata from URL
load(url("https://sites.google.com/a/wisc.edu/jed-frees/home/documents/data.RData"))
load(url("https://sites.google.com/a/wisc.edu/jed-frees/home/documents/dataout.RData"))

train <- data[,c(1:2,9:14,21:23,25,34:35,37)] # Only IM and CN claims are used
rm(data)
head(train)

test <- dataout[,c(1:2,9:14,21:23,25,34:35,37)]
rm(dataout)
head(test)

#### Preliminary analysis ####

# Measures of Association
attach(train)
cor.test(FreqCN, FreqIM, method="kendall") 
cor.test(FreqCN, FreqIM, method="spearman") 
cor.test(FreqCN, FreqIM, method="pearson") 

# Frequency table
FCN <- train$FreqCN
FIM <- train$FreqIM

FCN[FCN>=5] <- 5
FIM[FIM>=5] <- 5

table(FCN, FIM)

glm.poisIM  <- glm(train$FreqIM~.,data=train[-c(1,2,4,12:15)], family="poisson")               # Poisson model for IM
glm.poisCN  <- glm(train$FreqCN~.,data=train[-c(1,2,4,9:12,15)], family="poisson")             # Poisson model for CN

Rprox <- as.data.frame(cbind(train$Year, train$PolicyNum, train$FreqCN,
                             train$FreqIM, fitted(glm.poisCN), fitted(glm.poisIM)))
colnames(Rprox) <- c("Year", "PolicyNum", "FreqCN", "FreqIM", "nuCN", "nuIM")

dsumCN   <- 0
odsum    <- 0
dsumIM   <- 0
nsumCN   <- 0
onsum    <- 0
nsumIM   <- 0
n        <- 0
polnums  <- unique(Rprox$PolicyNum)
M <- length(polnums)

for (i in 1:M) {
  incn <- sum(Rprox$PolicyNum == polnums[i])
  vec  <- c(Rprox$FreqCN[(n+1):(n+incn)]-Rprox$nuCN[(n+1):(n+incn)],
            Rprox$FreqIM[(n+1):(n+incn)]-Rprox$nuIM[(n+1):(n+incn)])
  nec  <- c(Rprox$nuCN[(n+1):(n+incn)], Rprox$nuIM[(n+1):(n+incn)])
  
  mat     <- outer(vec, vec)
  nat     <- outer(nec, nec)
  
  dsumCN  <- dsumCN  + sum(diag(mat)[1:incn])
  dsumIM  <- dsumIM  + sum(diag(mat)) - sum(diag(mat)[1:incn])
  odsum   <- odsum   - sum(diag(mat)) + sum(mat)
  
  nsumCN  <- nsumCN  + sum(diag(nat)[1:incn])
  nsumIM  <- nsumIM  + sum(diag(nat)) - sum(diag(nat)[1:incn])
  onsum   <- onsum   - sum(diag(nat)) + sum(nat)
  
  n <- n + incn }

r   <- onsum/odsum
wIM <- sum(Rprox$nuIM)/(dsumIM -nsumIM/r)
wCN <- sum(Rprox$nuCN)/(dsumCN -nsumCN/r)

# Test for overdispersion

phihat_IM <- sum(residuals(glm.poisIM, type="pearson")^2)/df.residual(glm.poisIM) # phi_hat of IM
phihat_CN <- sum(residuals(glm.poisCN, type="pearson")^2)/df.residual(glm.poisCN) # phi_hat of CN

resIM <- ((fitted(glm.poisIM)-train$FreqIM)^2 - train$FreqIM)/ fitted(glm.poisIM)
mean(resIM)                                         # phi_tilde of IM
mean(resIM)/sd(resIM)*sqrt(length(resIM))           # test statistic for phi_tilde of IM

resCN <- ((fitted(glm.poisCN)-train$FreqCN)^2 - train$FreqCN)/ fitted(glm.poisCN)
mean(resCN)                                         # phi_tilde of CN
mean(resCN)/sd(resCN)*sqrt(length(resCN))           # test statistic for phi_tilde of CN

library(AER) # Quicker way to test the significance of phi_tilde
dispersiontest(glm.poisIM)  
dispersiontest(glm.poisCN)  

#### Model calibration ####

glm.nbIM <- glm.nb(train$FreqIM~.,data=train[-c(1,2,4,12:15)],control=glm.control(maxit=50))   # NB model for IM
glm.nbCN <- glm.nb(train$FreqCN~.,data=train[-c(1,2,4,9:12,15)],control=glm.control(maxit=50)) # NB model for CN

require(pscl)
zip.IM <- zeroinfl(train$FreqIM~.,data=train[-c(1,2,4,12:15)], dist = "poisson")               # ZIP model for IM
zip.CN <- zeroinfl(train$FreqCN~.,data=train[-c(1,2,4,9:12,15)], dist = "poisson")             # ZIP model for CN

# optimization of Shi and Valdez (2014) for multivariate frequency with Frank copula
copEst <- function(x1, x2, n1, n2, init.alpha1, init.alpha2, init.phi) { 
  x1 <- cbind(rep(1,nrow(x1)),x1)
  colnames(x1)[1] <- "intercept"
  x2 <- cbind(rep(1,nrow(x2)),x2)
  colnames(x2)[1] <- "intercept"
  
  # marginal likelihood of multivariate NB distribution
  "negll.frankcop" <- function(parm) {
    e1 <- ncol(x1);
    e2 <- ncol(x2);
    cop <- frankCopula(parm[18])
    lambdaIM <- exp(as.matrix(x1) %*% parm[1:e1] )
    lambdaCN <- exp(as.matrix(x2) %*% parm[(e1+1):(e1+e2)] )
    c  <- ppois(cbind(n1  , n2  ), cbind(lambdaIM,lambdaCN))
    c1 <- ppois(cbind(n1-1, n2-1), cbind(lambdaIM,lambdaCN))
    likfrank <- pCopula(c,cop)+pCopula(c1,cop)-pCopula(cbind(c[,1],c1[,2]),cop)-pCopula(cbind(c[,2],c1[,1]),cop)
    return(-sum(log(pmax(likfrank,10^-10))))
  } 
  
  init.est <- as.vector(c(init.alpha1,init.alpha2, init.phi))
  # init.est <- as.vector(c(init.alpha,init.w))
  
  fit.frankcop <- optim(init.est, negll.frankcop, method="BFGS", hessian=TRUE)
  parm.hat <- fit.frankcop$par
  loglik.frankcop <- -fit.frankcop$value
  inv.frankcop.Hess <- solve(fit.frankcop$hessian);
  parm.se <- sqrt(diag(inv.frankcop.Hess));
  # put together the model with the est, se, t, pval, AIC, BIC
  dfe <- length(n1+n2-length(parm.hat));
  t_ratio<-parm.hat/parm.se;
  #test if diff. from 1 t_ratio[1:3]<-(parm.hat[1:3]-1)/parm.se[1:3];
  pval <- pf(t_ratio*t_ratio,df1=1,df2=dfe,lower.tail=F);
  ttable <- cbind(parm.hat,parm.se,t_ratio,pval) 
  ttable <- round(ttable,digits=4)
  
  rownames(ttable)<- c(colnames(x1),colnames(x2),"phi")
  # rownames(ttable)<- c(colnames(x),"weight")
  colnames(ttable)<- c("estimate", "std error", "t-val","Pr>|t|");
  
  #  AIC<- 2*negll.frankcop(parm.hat) + 2*length(parm.hat);
  #  BIC<- 2*negll.frankcop(parm.hat) + log(length(n))*length(parm.hat);
  #  loglik <- negll.frankcop(parm.hat)
  return(list(ttable=ttable, coef=parm.hat));
}

system.time(
  copM <- copEst(x1=train[-c(1,2,4,12:15)]  ,n1=train$FreqIM,init.alpha1=coefficients(glm.poisIM), 
                 x2=train[-c(1,2,4,9:12,15)],n2=train$FreqCN,
                 init.alpha2=coefficients(glm.poisCN), init.phi=1.07)) 

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
    data <- cbind(id,w1*exp(reg_eqn1),w1*n1,w2*exp(reg_eqn2),w2*n2);
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
  NBm_Joint <- NBEst2(x1=train[-c(1,2,4,12:15)],n1=train$FreqIM,init.alpha1=coefficients(glm.poisIM), w1=wIM,
                     x2=train[-c(1,2,4,9:12,15)],id=train$PolicyNum, n2=train$FreqCN,
                     init.alpha2=coefficients(glm.poisCN), w2=wCN, r=r)) 

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
IMtable <- round(cbind(GLMtable(glm.poisIM)[,-2], GLMtable(zip.IM)[1:9,-2], copM$ttable[1:9,c(1,4)], 
                       NBm_Joint$ttable[1:9,c(1,4)]),digits=4)
colnames(IMtable) <- c("Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value")

kable(IMtable, booktabs = T, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{=========}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")

# Summary of regression coefficients for CN claim from four models
CNtable <- round(cbind(GLMtable(glm.poisCN)[,-2], GLMtable(zip.CN)[1:8,-2], copM$ttable[10:17,c(1,4)]+cbind(rep(0,8), c(rep(0,7),0.0001)),
                       NBm_Joint$ttable[-(1:9),c(1,4)]),digits=4)
colnames(CNtable) <- c("Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value", "Estimate", "p-value")

kable(CNtable, booktabs = T, linesep = c("","","","","","","","","","\\hline"),
      bottomrule="\\hhline{=========}",escape=FALSE) %>% kable_styling(latex_options = "hold_position")

#### Out-of-sample validation ####

x1 <- train[-c(1,2,4,12:15)]
x1 <- cbind(rep(1,nrow(x1)),x1)
n1 <- train$FreqIM
#w1 <- 1/phihat_IM
w1 <- wIM  
x2 <- train[-c(1,2,4,9:12,15)]
x2 <- cbind(rep(1,nrow(x2)),x2)
id <- train$PolicyNum
n2 <- train$FreqCN
#w2 <- 1/phihat_CN
w2 <- wCN


# calculation of bonus-malus factor based on the proposed model for frequency
Nreg_eqn1 <- as.matrix(x1) %*% as.matrix(NBm_Joint$coef[1:9])
Nreg_eqn2 <- as.matrix(x2) %*% as.matrix(NBm_Joint$coef[10:17])
Ndata <- cbind(id,exp(Nreg_eqn1),n1, exp(Nreg_eqn2),n2);
colnames(Ndata)[2] <- "nv1";
colnames(Ndata)[4] <- "nv2";


Npost     <- aggregate(n1~id,Ndata,sum)      # aggregate actual   IM claim counts for years
Npost$nv1 <- aggregate(nv1~id,Ndata,sum)[,2] # aggregate expected IM claim counts for years
Npost$n2  <- aggregate(n2~id,Ndata,sum)[,2]  # aggregate actual   CN claim counts for years
Npost$nv2 <- aggregate(nv2~id,Ndata,sum)[,2] # aggregate expected CN claim counts for years

Npost$nweight <- (w1*Npost$n1 +w2*Npost$n2 + r) / (w1*Npost$nv1 + w2*Npost$nv2 + r)
# bonus-malus factor for a policyholder on frequency
colnames(Npost)[1] <- "PolicyNum"

# attach the bonus-malus factor for each policyholder on the test set
Ptest <- merge(x = test, y = Npost, by = "PolicyNum", all.x = TRUE)
Ptest$nweight[is.na(Ptest$nweight)] <- 1

xt1 <- test[-c(1,2,4,12:15)]
xt1 <- cbind(rep(1,nrow(xt1)),xt1)
xt2 <- test[-c(1,2,4,9:12,15)]
xt2 <- cbind(rep(1,nrow(xt2)),xt2)

n1_poispred <- exp(as.matrix(xt1) %*% coefficients(glm.poisIM)) # IM frequency premium with Poisson model
n1_coppred <- exp(as.matrix(xt1) %*% as.matrix(copM$coef[1:9]))     # IM frequency premium with copula model
n1_zippred   <- predict(zip.IM, xt1)                            # IM frequency premium with ZIP model
n1_propred <- exp(as.matrix(xt1) %*% as.matrix(NBm_Joint$coef[1:9])
                  )*Ptest$nweight                               # IM frequency premium with Proposed model

n2_poispred <- exp(as.matrix(xt2) %*% coefficients(glm.poisCN)) # CN frequency premium with Poisson model
n2_coppred <- exp(as.matrix(xt2) %*% as.matrix(copM$coef[10:17]))    # CN frequency premium with copula model
n2_zippred   <- predict(zip.CN, xt2)                            # CN frequency premium with ZIP model
n2_propred <- exp(as.matrix(xt2) %*% as.matrix(NBm_Joint$coef[10:17])
                  )*Ptest$nweight                               # CN frequency premium with Proposed model

n_poispred  <- c(n1_poispred, n2_poispred)
n_coppred  <- c(n1_coppred, n2_coppred)
n_zippred    <- c(n1_zippred, n2_zippred)
n_propred  <- c(n1_propred, n2_propred)
n_actual <- c(test$FreqIM, test$FreqCN)


# RMSE (root mean squared errors) of four models
RMSE_pois    <- sqrt(mean((n_poispred - n_actual)^2))
RMSE_cop       <- sqrt(mean((n_coppred - n_actual)^2))
RMSE_zip      <- sqrt(mean((n_zippred   - n_actual)^2))
RMSE_proposed <- sqrt(mean((n_propred - n_actual)^2))

RMSE_pois
RMSE_cop
RMSE_zip
RMSE_proposed


# MAE (mean absoulte errors) for four models
MAE_pois    <- mean(abs(n_poispred - n_actual))
MAE_cop       <- mean(abs(n_coppred - n_actual))
MAE_zip      <- mean(abs(n_zippred - n_actual))
MAE_proposed <- mean(abs(n_propred - n_actual))

MAE_pois
MAE_cop
MAE_zip
MAE_proposed

Npost$n <- Npost$n1 + Npost$n2
Npost$n[Npost$n> 2 ]<- "> 2"
Npost$n[Npost$n==0] <- "= 0"
Npost$n[Npost$n==1] <- "= 1"
Npost$n[Npost$n==2] <- "= 2"
#Npost$n[Npost$n> 2] <- "> 2"

Npost$n <- factor(Npost$n, levels = c("= 0","= 1", "= 2", "> 2"))
library(ggplot2)

ggplot(Npost, aes(x=n, y=nweight, color=n)) +  geom_boxplot() +
  labs(x="Number of claims observed", y="Posterior rating factors",
       title = "Posterior rating factors for frequency\nunder the proposed model") +
  theme(legend.position = c(.10, .82), plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), axis.text.x=element_blank()) +
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="cornflowerblue") +
  scale_y_continuous(labels = scales::percent) + guides(fill=guide_legend(ncol=2))
