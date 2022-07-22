
##############################
## Written by Janie Coulombe
## For manuscript entitled
## 'Estimating individualized treatment rules in longitudinal studies with covariate-driven observation times'
## Co-authors Erica EM Moodie, Susan M Shortreed, Christel Renoux
## July 2022
##############################



rm(list=ls(all=TRUE))


 
## Packages to download

library(survival)
library(splines)


## FUNCTION SIMUL, can be used to simulate the datasets and to estimate the blip function (i.e., its coefficients)
## using the 6 different estimators



simul<- function(sampsize=500, nbsimul=1 ,print ){
    ## Change according to sample size (sampsize) and nb of simulations (nsbimul) desired
    ## print=TRUE shows the time needed to perform only the analysis part (fit the regression models, etc.) for DW1 estimator


## Maximum follow-up time TAU
TAU<-1
 
## Parameters outcome model
beta4<-0.4   # age
beta5<-0.05  # sex
beta6<- -0.6 # height
bint<- 0.5   # Q*treatment interaction
bint2<- -1   # Age*treatment interaction
betaX<- -2

## For error term in outcome model
sigmaepsilond<- 0.1
meanphid=0; 
sigmaphid=0.2; 

## Parameters related to mediator Z
mu1d=4; 
sigma2_1d=2; 
mu2d=2; 
sigma2_2d=1; 
betaZ<-  2.5

## Treatment model parameters
beta0=0.5 
beta1=0.55  
beta2=-0.2  
beta3=-1

# For defining the 3 confounders "age, male, height"
agemean=1;  
agesd=1; 	
p_male=0.55; 
heightmean=0; 
heightsd=1;  

## Outcome observation model (gamma parameters)
gamma1<- 0.3
gamma2<- -0.6
gamma3<- -0.4
gamma4<- -0.3
  
nbvisits<-c() ## to count visits

## Loop, ran "nbsimul" times 

for( S in 1: nbsimul){


      print(S)

	## create empty matrix with time-dependent variables to keep
	mat<- matrix(NA, nrow=sampsize*TAU*100, ncol= 11) 
	colnames(mat)<- c('ID','time','X','Y','age','sex','height','Qvar','Fother','visit','Z')
	count<- 0

	for( i in 1:sampsize){

            ## Define confounders
      	age<-rnorm(mean=agemean,sd=agesd,n=1) 
      	sex<-rbinom(p=p_male,n=1,size=1)
            height<-rnorm(mean=heightmean,sd=heightsd,n=1)

	      ## Treatment depending on confounders
	      p_treated <- exp(beta0+beta1*age +beta2*sex+beta3* height )/(1+exp(beta0+beta1*age +beta2*sex+beta3* height )) 
		visitcount<-0

 		for(t in 1:(TAU*100)){ 

      			count<- count+1

				mat[count,5] <- age
				mat[count, 6]<- sex
				mat[count, 7]<- height
				mat[count ,3] <- rbinom(p=p_treated,n=1,size=1) ## X: the exposure 

				## Mediator
				Z1a<-rnorm(mean=mu1d,sd=sqrt(sigma2_1d),n=1) ## if X=0
      			Z1b<-rnorm(mean=mu2d,sd=sqrt(sigma2_2d),n=1) ## if X=1
            		mat[count, 11]<- ifelse(mat[count,3]==1,Z1b,Z1a)
				center<-ifelse(mat[count,3]==1,mu2d,mu1d)

				## Simulate tailoring variable Q
				mat[count, 8]<- rnorm(mean=0.5, sd=0.5, n=1) 
		
				## Simulate outcome
				alpha0<- sqrt(t/100)
				phi<-rnorm(mean=meanphid,sd=sigmaphid,n=1)
            		epsilon<-rnorm(mean=phi,sd=sigmaepsilond,n=1)

				mat[count, 4]<- alpha0 + betaX*mat[count ,3] + beta4*mat[count ,5] + beta5*mat[count ,6] + beta6*mat[count ,7] +
			               	betaZ*(mat[count,11]-center) + bint*(mat[count ,8])*mat[count ,3] + bint2* mat[count ,3]*mat[count ,5] + epsilon
		
				# Add ID into the matrix, as well as time
				mat[count, 1]<- i
				mat[count, 2]<- t/100
			
				## Fother variable (a pure predictor) - this one was not used in the following simulations... could be removed
				mat[count, 9] <- rbinom(p= 0.5, n=1, size=1)  

				ratei<-  exp(  gamma1*mat[count,3]  + gamma2*mat[count,11] + gamma3* mat[count, 6] + gamma4*mat[count,7] )*0.1
             		ratei[ratei>1]<-1 # remove rates higher than 1, should act as a probability of jump over tiny time increment
            		mat[count, 10]<-  rbinom(size=1,n=1,prob= ratei)

	 			ifelse( mat[count, 10]==1 , visitcount<- visitcount + 1, visitcount<- visitcount)  ## to count the number of visits
						} ## End of the loop for the follow-up of one individual

			nbvisits<- append(nbvisits, visitcount)
					} ## End of the loop for creating the dataset (loop ran once for each individual)
 
      start_time2 <- Sys.time()
 
	## Keep only the visits
	data<- data.frame(mat)
	datvisit<- data[data$visit==1,]
	data[data$visit==0,]$Y<-NA ## Put outcome as missing if no visit

	## Visit intensity model (need counting process format)
	data$t1<-data$time-0.01
    	data$t2<-data$time

  	gamma<-coxph(Surv(data$t1, data$t2, data$visit )~ data$X +data$Z+ data$sex + data$height )$coef 
	## 2 other models, for the wrongly specified visit models:
   	gamma2b<-coxph(Surv(data$t1, data$t2, data$visit )~ data$X + data$Z  )$coef 
 	gamma3b<-coxph(Surv(data$t1, data$t2, data$visit )~ data$X + data$sex )$coef 
	
      data$rho_i<-    exp(gamma[1]*data$X+ gamma[2]*data$Z+gamma[3]*data$sex + gamma[4]*data$height)    ## right model
      data$rho_i2<-   exp(gamma2b[1]*data$X+ gamma2b[2]*data$Z )    ## wrong model but no expected bias
      data$rho_i3<-   exp(gamma3b[1]*data$X+ gamma3b[2]*data$sex )    ## wrong model and expected bias
 
 	## Compute the propensity scores and two types of IPTW (one correct iptw, one wrong iptw2) 
	
	ps2<-predict(glm(data$X ~ I(data$age^2) + data$sex + I(sin(data$height^2)), family='binomial'),type='response')
	data$iptw2<- 1/ifelse(data$X==1, ps2, (1-ps2))	

	ps<-predict(glm(data$X ~ data$age + data$sex + data$height, family='binomial'),type='response')
	data$iptw<- 1/ifelse(data$X==1, ps, (1-ps))	

 
	############################################
	## Compute the different blip functions   ##
      ############################################

      ## DOUBLY WEIGHTED DW1 ##
      (coef_DW1_x <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw*1/rho_i, data=data)$coef[5] )
	(coef_DW1_age <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw*1/rho_i, data=data)$coef[10]) 
	(coef_DW1_Q <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw*1/rho_i, data=data)$coef[11] )

      end_time2 <- Sys.time()
      if(print=='TRUE'){print(end_time2-start_time2) }

      ## DOUBLY WEIGHTED DW2 ##
	(coef_DW2_x <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X  + height, weight= iptw*1/rho_i2, data=data)$coef[5] )
	(coef_DW2_age <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + height, weight= iptw*1/rho_i2, data=data)$coef[9]) 
	(coef_DW2_Q <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X  + height, weight= iptw*1/rho_i2, data=data)$coef[10] )
   
      ## DOUBLY WEIGHTED DW3 ##	 
	(coef_DW3_x <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw2*1/rho_i2, data=data)$coef[5] )
	(coef_DW3_age <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw2*1/rho_i2, data=data)$coef[10]) 
	(coef_DW3_Q <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw2*1/rho_i2, data=data)$coef[11] )
 
      ## DOUBLY WEIGHTED DW4 ##
	(coef_DW4_x <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw*1/rho_i3, data=data)$coef[5] )
	(coef_DW4_age <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw*1/rho_i3, data=data)$coef[10]) 
	(coef_DW4_Q <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw*1/rho_i3, data=data)$coef[11] )
 
	## OLS ESTIMATOR ##
	(coef_OLS_x <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height  , data=data)$coef[5] )
	(coef_OLS_age <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height  , data=data)$coef[10]) 
	(coef_OLS_Q <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height  , data=data)$coef[11] )

      ## IPT ESTIMATOR ##
	(coef_IPT_x <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw , data=data)$coef[5] )
	(coef_IPT_age <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw , data=data)$coef[10]) 
	(coef_IPT_Q <-lm(Y ~ bs( time, degree=3)+ X +  age*X + age + Qvar + Qvar*X + sex + height, weight= iptw , data=data)$coef[11] )

      ## add a <write> function to print them somewhere at each round of the loop... 
 } ## End loop for one simulation

 
} ## End function 


#########################################################################################

## Running times for sample sizes of 250, 500, 50000 
 ## for 1 simulation, including the dataset creation:

for(sizes in c(250, 500, 50000 )){
 start_time <- Sys.time()
 simul(sampsize=sizes, nbsimul=1, print='FALSE')
 end_time <- Sys.time()
 print(end_time-start_time)
 }

[1] 1
Time difference of 2.041787 secs
[1] 1
Time difference of 3.788557 secs
[1] 1
Time difference of 11.60191 mins

#########################################################################################

## Running times for sample size of 250, 500, 50000 
 ## for the data analysis part only, for 1 simulation, 
  ## i.e., only the time needed for fitting the regression models for DW1 estimator

for(sizes in c(250, 500, 50000)){
 simul(sampsize=sizes, nbsimul=1, print='TRUE')
 }

[1] 1
Time difference of 0.5056951 secs
[1] 1
Time difference of 1.0131 secs
[1] 1
Time difference of 4.120133 mins

##########################################################################################


