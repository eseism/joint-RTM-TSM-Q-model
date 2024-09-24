
# Joint RTM-TSM attenuation model

This Rscript can be used to produce the to reproduce the main results in the journal article ''Joint Bayesian Hierarchical Model for Lg Wave Attenuation and Station Site Response Using Both RTM and TSM paths in The Middle East".

The model is a joint Bayesian shared parameter framework that simulataneoulsy models Reverse Two Station (RTM) paths and Two Station  (TSM) paths while accounting for station site responses. 


## Authors

- Saikat Nandy [@SaikatNandy](https://github.com/SaikatNandy)


## Installation


```bash
#install.packages("readr",dependencies=TRUE) #install devtools from CRAN
#library(devtools) #load devtools
libraries=c("boot","rstan","dplyr","Metrics","rstudioapi","ggplot2","hrbrthemes","Matrix","maps","gridExtra")
install.packages(libraries,dependencies=TRUE)
sapply(libraries,require,character=TRUE)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)

```

    
## Data
The following data files are used to reproduce the main results in the journal article ''Joint Bayesian Hierarchical Model for Lg Wave Attenuation and Station Site Response Using Both RTM and TSM paths in The Middle East".
Three datafiles are used an input data one for the RTM paths, one for TSM paths, and the third for station site responses. 

Column information of the RTM dataset: \
1: Event Longitude \
2: Event Latitude\
3: Station Longitude\
4: Station Latitude\
5: distance between Event 1 and Station 1 (See Figure 3 in manuscript for reference)\
6: distance between Event 1 and Station 2 (See Figure 3 in manuscript for reference)\
7: distance between Event 2 and Station 1 (See Figure 3 in manuscript for reference)\
8: distance between Event 2 and Station 2 (See Figure 3 in manuscript for reference)\
9: Total distance travelled by each ray path\
10: Spectral amplitude Ratio\
11-1897: Subdistance Matrix

```bash
RTM <- read_table("ME_RTM.out", col_names = FALSE)
subdis = sapply(X = 1:1887, FUN = function(x) paste("subdis", x, sep = "")) 
RTM_col_names = c("stlo1", "stla1", "stlo2", "stla2", "d11", "d12", "d21", "d22", "L", "AR", subdis)
names(RTM) = RTM_col_names
```

Column information of the TSM dataset: \
1: Event Longitude \
2: Event Latitude\
3: Station Longitude\
4: Station Latitude\
5: distance between Event and Station 1 (See Figure 2 in manuscript for reference)\
6: distance between Event and Station 2 (See Figure 2 in manuscript for reference)\
7: Total distance travelled by each ray path\
8: Spectral amplitude Ratio\
9-1895: Subdistance Matrix

```bash
TSM <- read_table("ME_TSM.out", col_names = FALSE)
subdis = sapply(X = 1:1887, FUN = function(x) paste("subdis", x, sep = "")) 
TSM_col_names = c("stlo1", "stla1", "stlo2", "stla2","d1", "d2", "L", "AT", subdis)
names(TSM) = TSM_col_names
```

Column information of the Site Response dataset: \
1: Event Longitude \
2: Event Latitude\
3: Station Longitude\
4: Station Latitude\
5: Relative distance between stations (See Figure 3 in manuscript for reference)\
6: Relative site response on the log scale\
7-1893: Subdistance Matrix

```bash
SSmeshFile <- read_table("ME_SS.out", col_names = FALSE)
subdis = sapply(X = 1:1887, FUN = function(x) paste("subdis", x, sep = "")) 
SS_col_names = c("stlo1", "stla1", "stlo2", "stla2", "L", "Y_S", subdis)
names(SSmeshFile) = SS_col_names
```

Incidence Matrices with entries for station pairs corresponding to station site responses
```bash
Phi_T=read_table("phiT.txt",col_names = F)
Phi_S=read_table("phiS.txt",col_names = F)
```

## Usage/Examples

Prepare the response vectors and sub distance matrices


```javascript
m=0.5  // Geometric spreading term
v=3.5  // group velocity
f=1    // frequency

Y_R = 0.5 * log(RTM$AR * (RTM$d12*RTM$d21/RTM$d11/RTM$d22)^m)     // log amplitude ratio for RTM paths
X_R = as.matrix(-pi/v*f*RTM[,11:1897])   // subdistance matrix for RTM paths
Y_T = log(TSM$AT *(TSM$d12/TSM$d11)^m)  // log amplitude ratio for TSM paths
X_T = as.matrix(-pi/v*f*TSM[,9:1895])   // subdistance matrix for TSM paths 

Y_S=SSmeshFile$Y_S   // log relative site response
```

Filter the subdistance matrix, to discard cells with no ray paths passing through it. We keep cells with at least one RTM or TSM ray paths passing through it. 

```bash
X = rbind(X_R, X_T)
columns.Sum = colSums(X)
columns.Sum = (columns.Sum != 0)
X_R = X_R[,columns.Sum]
X_T = X_T[,columns.Sum]
```

Draw the MCMC samples
```bash
model_dat <- list(n_R = nrow(X_R),
                   n_T = nrow(X_T),
                   n_S = length(Y_S),
                   p = ncol(X_R),
                   r = ncol(Phi_T),
                   Y_R = Y_R,
                   Y_T = Y_T,
                   X_R = as.matrix(X_R),
                   X_T = as.matrix(X_T),
                   Y_S = Y_S,
                   Phi_T = as.matrix(Phi_T),
                   Phi_S = as.matrix(Phi_S))

fit <- stan(file = "joint_model.stan", chains = 1, data = model_dat, iter=8000, warmup=5000)
```

Compute the posterior estimates for Q
```bash
pars = rstan::extract(fit1)

est=summary(fit,pars="beta",probs=c(0.025,0.50,0.75,0.975))$summary
beta.est=est[,1]
beta.se=est[,3]
QQ=1/pars$beta
Q.hat=apply(QQ,2,mean)
Q.hat.sd=apply(QQ,2,sd)
Q.hat.lower=apply(QQ,2,quantile,probs=c(0.025))
Q.hat.upper=apply(QQ,2,quantile,probs=c(0.975))
Q.hat.75=apply(QQ,2,quantile,probs=c(0.75))
```

Visualize the Lg Q estimates for the study area
```bash
meshFile <- read_table("mesh.out", col_names = FALSE)
meshFile=as.matrix(stations[columns.Sum,-1])

dat_joint=data.frame(X1=meshFile[,1],X2=meshFile[,2],Q.est=Q.hat,beta.est=beta.est,beta.sd=est[,3],beta.025=est[,4],beta0.975=est[,7],Q.025=Q.hat.lower,Q.975=Q.hat.upper,Q.sd=Q.hat.sd,Q.075=Q.hat.75,beta.075=est[,6])

q.est.plot.joint =   ggplot() + geom_raster(data=dat_joint, aes(x=X1, y = X2, fill = Q.est))+scale_fill_viridis_c()+ggtitle("Estimated Q, using joint RTM-TSM")+labs(x="Longitude",y="Latitude")+theme(legend.title=element_text(hjust=0.5),legend.position = c(.85, .1),legend.justification = c("left", "bottom"), legend.box.just = "center", legend.margin = margin(4, 4, 4, 4))
q.sd.plot.joint =    ggplot() + geom_raster(data=dat_joint, aes(x=X1, y = X2, fill = Q.sd))+scale_fill_viridis_c()+ggtitle("Standard deviation in Q estimates, joint RTM-TSM")+labs(x="Longitude",y="Latitude")+theme(legend.title=element_text(hjust=0.5),legend.position = c(.85, .1),legend.justification = c("left", "bottom"), legend.box.just = "center", legend.margin = margin(4, 4, 4, 4))
q.lower.plot.joint = ggplot() + geom_raster(data=dat_joint, aes(x=X1, y = X2, fill = Q.025))+scale_fill_viridis_c()+ggtitle("Lower 95% CI for Q estimates, joint RTM-TSM")+labs(x="Longitude",y="Latitude")+theme(legend.title=element_text(hjust=0.5),legend.position = c(.85, .1),legend.justification = c("left", "bottom"), legend.box.just = "center", legend.margin = margin(4, 4, 4, 4))
q.upper.plot.joint = ggplot() + geom_raster(data=dat_joint, aes(x=X1, y = X2, fill = Q.975))+scale_fill_viridis_c()+ggtitle("Upper 95% CI for Q estimates, joint RTM-TSM")+labs(x="Longitude",y="Latitude")+theme(legend.title=element_text(hjust=0.5),legend.position = c(.85, .1),legend.justification = c("left", "bottom"), legend.box.just = "center", legend.margin = margin(4, 4, 4, 4))



grid.arrange(q.est.plot.joint,q.sd.plot.joint,q.lower.plot.joint,q.upper.plot.joint,ncol=3)

```

Figures in the manuscript were generated using GMT.
