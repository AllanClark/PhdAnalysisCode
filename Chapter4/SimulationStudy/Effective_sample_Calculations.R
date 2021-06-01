
require(coda)

#The output for one of the data sets will be of dimension 500*10000 by 12
#extract a 10 000 by 12 matrix and then perform the calculations

effective_Sample_calcs<-function(samples_SS, times_SS=NULL){
     #The function calculates the Effective sample size as well as effective sampling rate
     #for the jags samples, drum method as well as Polya-Gamma method

     effectiveSS<-as.mcmc(samples_SS)

     x<-effectiveSize(effectiveSS) #effective samples sizes
     xt<-x #rescaled effective sample sizes

     xt[1:4]<-x[1:4]/times_SS[1]
     xt[5:8]<-x[5:8]/times_SS[2]
     xt[9:12]<-x[9:12]/times_SS[3]

     alpha0<-c(x[1], x[5], x[9])
     alpha1<-c(x[2], x[6], x[10])
     beta0<-c(x[3], x[7], x[11])
     beta1<-c(x[4], x[8], x[12])

     alpha0t<-c(xt[1], xt[5], xt[9])
     alpha1t<-c(xt[2], xt[6], xt[10])
     beta0t<-c(xt[3], xt[7], xt[11])
     beta1t<-c(xt[4], xt[8], xt[12])

     eff_all_coefs<-c(alpha0, alpha1, beta0, beta1)
     efft_all_coefs<-c(alpha0t, alpha1t, beta0t, beta1t)

     list( ESS=eff_all_coefs, RESS=efft_all_coefs, 
           med_ESS=c(median(x[1:4]),median(x[5:8]),median(x[9:12])),
           med_RESS=c(median(xt[1:4]),median(xt[5:8]),median(xt[9:12])),
           min_ESS=c(min(x[1:4]),min(x[5:8]),min(x[9:12])),
           max_ESS=c(max(x[1:4]),max(x[5:8]),max(x[9:12])),
           min_RESS=c(min(xt[1:4]),min(xt[5:8]),min(xt[9:12])),
           max_RESS=c(max(xt[1:4]),max(xt[5:8]),max(xt[9:12]))
          )
}

#effective_Sample_calcs(samples_SS=xs, xtimes)

#acfplot_samples=function(colindex, lm=100){
#     plot(acf(xs[,colindex[1]], lag.max=lm, plot=F)$acf, type="l", ylim=c(-1,1), ylab="Sample ACF", xlab="Lag")
#     lines(acf(xs[,colindex[2]], plot=F,lag.max=lm)$acf, col="blue")
#     lines(acf(xs[,colindex[3]], plot=F, lag.max=lm)$acf, col="red")
#     abline(h=0, lty=2)
#     abline(v=0, lty=2)
#}

#par(mfrow=c(2,2), mar=c(4,4,0,1))

#acfplot_samples(c(1,5,9), 500)
#acfplot_samples(c(2,6,10), 500)
#acfplot_samples(c(3,7,11), 500)
#acfplot_samples(c(4,8,12), 500)

#applying to all sets
#----------------------------------------------------------------------------------------------


require(coda)

#The output for one of the data sets will be of dimension 500*10000 by 12
#extract a 10 000 by 12 matrix and then perform the calculations

effective_Sample_calcs<-function(samples_SS, times_SS=NULL){
     #The function calculates the Effective sample size as well as effective sampling rate
     #for the jags samples, drum method as well as Polya-Gamma method

     effectiveSS<-as.mcmc(samples_SS)

     x<-effectiveSize(effectiveSS) #effective samples sizes
     xt<-x #rescaled effective sample sizes

     xt[1:4]<-x[1:4]/times_SS[1]
     xt[5:8]<-x[5:8]/times_SS[2]
     xt[9:12]<-x[9:12]/times_SS[3]

     alpha0<-c(x[1], x[5], x[9])
     alpha1<-c(x[2], x[6], x[10])
     beta0<-c(x[3], x[7], x[11])
     beta1<-c(x[4], x[8], x[12])

     alpha0t<-c(xt[1], xt[5], xt[9])
     alpha1t<-c(xt[2], xt[6], xt[10])
     beta0t<-c(xt[3], xt[7], xt[11])
     beta1t<-c(xt[4], xt[8], xt[12])

     eff_all_coefs<-c(alpha0, alpha1, beta0, beta1)
     efft_all_coefs<-c(alpha0t, alpha1t, beta0t, beta1t)

     list( ESS=eff_all_coefs, RESS=efft_all_coefs, 
           med_ESS=c(median(x[1:4]),median(x[5:8]),median(x[9:12])),
           med_RESS=c(median(xt[1:4]),median(xt[5:8]),median(xt[9:12])),
           min_ESS=c(min(x[1:4]),min(x[5:8]),min(x[9:12])),
           max_ESS=c(max(x[1:4]),max(x[5:8]),max(x[9:12])),
           min_RESS=c(min(xt[1:4]),min(xt[5:8]),min(xt[9:12])),
           max_RESS=c(max(xt[1:4]),max(xt[5:8]),max(xt[9:12]))
          )
}


#Run for the first set of data (p=0.5)
Data_sets<-c("Logit_simulations_n50_J2_case1","Logit_simulations_n50_J5_case1",
             "Logit_simulations_n100_J2_case1","Logit_simulations_n100_J5_case1",
             "Logit_simulations_n50_J2_case2","Logit_simulations_n50_J5_case2",
             "Logit_simulations_n100_J2_case2","Logit_simulations_n100_J5_case2"
)


setwd("C:\\Users\\User\\Documents\\Research\\logit_occ\\OfficeCluster\\Simulations")


All_Results <- list(NULL)

for (idata in 1:8){

	filename<-paste(Data_sets[idata],".RData", sep="")
      load(filename)

	if (idata==1) { output_used = output1 }
	if (idata==2) { output_used = output2 }
	if (idata==3) { output_used = output3 }
	if (idata==4) { output_used = output4 }
	if (idata==5) { output_used = output5 }
	if (idata==6) { output_used = output6 }
	if (idata==7) { output_used = output7 }
	if (idata==8) { output_used = output8 }
	#if (idata==1) { output_used = output1 }

	R1<-matrix(0, nrow=500, ncol=12)
	R2<-R1
	R3<-R1[,1:3]
	R4<-R3
	R5<-R3
	R6<-R3
	R7<-R3
	R8<-R3

	for (i in 1:500){
		lo<-1+10000*(i-1)
		hi<-lo+10000-1
		xs<-output_used[ lo:hi, 1:12]
		xtimes<-output_used[ lo, 13:15]

		es<-effective_Sample_calcs(samples_SS=xs, xtimes)

		R1[i,]<-es$ESS
		R2[i,]<-es$RESS

		R3[i,]<-es$med_ESS
		R4[i,]<-es$med_RESS

	      R5[i,]<-es$min_ESS
		R6[i,]<-es$min_RESS

	      R7[i,]<-es$max_ESS
		R8[i,]<-es$max_RESS
	}

	Results<-cbind(R1, R2, R3, R4, R5, R6, R7, R8)

	All_Results[[idata]]<-Results
	save.image("Logit_simulations_All_Results.RData")

	print(idata)
}


rm( list=ls()[-1])
save.image("Logit_simulations_All_Results_condense.RData")

#par(mfrow=c(2,2), mar=c(4,4,0,1))
#for (i in 1:4){boxplot(cbind(R1[,(1:3)+3*(i-1)], R2[,(1:3)+3*(i-1)]), outline=F, col=c(rep("red",3),rep("yellow",3)))}
#par(mfrow=c(1,2), mar=c(4,4,0,1))
#boxplot(R3[,1:3])
#boxplot(R4[,1:3])




#Run for the second set of data (p=0.7)
#------------------------------------------------------

Data_sets<-c("Logit_simulations_n50_J2_case5","Logit_simulations_n50_J5_case5",
             "Logit_simulations_n100_J2_case5","Logit_simulations_n100_J5_case5",
             "Logit_simulations_n50_J2_case6","Logit_simulations_n50_J5_case6",
             "Logit_simulations_n100_J2_case6","Logit_simulations_n100_J5_case6"
)


setwd("C:\\Users\\User\\Documents\\Research\\logit_occ\\OfficeCluster\\Simulations")


All_Results <- list(NULL)

for (idata in 1:8){

	filename<-paste(Data_sets[idata],".RData", sep="")
      load(filename)

	if (idata==1) { output_used = output13 }
	if (idata==2) { output_used = output14 }
	if (idata==3) { output_used = output15 }
	if (idata==4) { output_used = output16 }
	if (idata==5) { output_used = output17 }
	if (idata==6) { output_used = output18 }
	if (idata==7) { output_used = output19 }
	if (idata==8) { output_used = output20 }

	R1<-matrix(0, nrow=500, ncol=12)
	R2<-R1
	R3<-R1[,1:3]
	R4<-R3
	R5<-R3
	R6<-R3
	R7<-R3
	R8<-R3

	for (i in 1:500){
		lo<-1+10000*(i-1)
		hi<-lo+10000-1
		xs<-output_used[ lo:hi, 1:12]
		xtimes<-output_used[ lo, 13:15]

		es<-effective_Sample_calcs(samples_SS=xs, xtimes)

		R1[i,]<-es$ESS
		R2[i,]<-es$RESS

		R3[i,]<-es$med_ESS
		R4[i,]<-es$med_RESS

	      R5[i,]<-es$min_ESS
		R6[i,]<-es$min_RESS

	      R7[i,]<-es$max_ESS
		R8[i,]<-es$max_RESS
	}

	Results<-cbind(R1, R2, R3, R4, R5, R6, R7, R8)

	All_Results[[idata]]<-Results
	save.image("Logit_simulations_All_Results2.RData")

	print(idata)
}


rm( list=ls()[-1])
save.image("Logit_simulations_All_Results_condense2.RData")



