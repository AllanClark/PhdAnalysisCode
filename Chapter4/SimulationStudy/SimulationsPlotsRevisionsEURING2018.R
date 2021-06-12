#setwd("~/Documents/AllanClark/Stats/Research/Clark_AE/PhD/Paper3/EURINGmanuscript/Analysis/SimulationStudy")
  
median.stats<-function(working){
  W<-working
  m1<-apply(W[,1:4], 1, function(x){median(x, na.rm=T)})
  m2<-apply(W[,5:8], 1, function(x){median(x, na.rm=T)})
  m3<-apply(W[,9:12], 1, function(x){median(x, na.rm=T)})
  m4<-apply(W[,13:16], 1, function(x){median(x, na.rm=T)})
  
  return( cbind(m1, m2, m3, m4) )
}

Data_sets<-c("Logit_simulations_n50_J3_case1","Logit_simulations_n50_J5_case1",
             "Logit_simulations_n100_J3_case1","Logit_simulations_n100_J5_case1",
             "Logit_simulations_n50_J3_case2","Logit_simulations_n50_J5_case2",
             "Logit_simulations_n100_J3_case2","Logit_simulations_n100_J5_case2",
             "Logit_simulations_n50_J3_case5","Logit_simulations_n50_J5_case5",
             "Logit_simulations_n100_J3_case5","Logit_simulations_n100_J5_case5",
             "Logit_simulations_n50_J3_case6","Logit_simulations_n50_J5_case6",
             "Logit_simulations_n100_J3_case6","Logit_simulations_n100_J5_case6",
             "Logit_simulations_n500_J5_case3", "Logit_simulations_n500_J10_case3", 
             "Logit_simulations_n500_J5_case4", "Logit_simulations_n500_J10_case4"
)

Data_sets<-paste(Data_sets,"_EURING_revision", sep="")

#----------------------------------------------------------------------------------
#Produce a plot of the ESS, ESR and Timing stats for the different methods
#consider n=50 n=100 and case 1
#----------------------------------------------------------------------------------

#pdf("ESS_ESR_n50_n_100_J3_case1.pdf") #1,3
#pdf("ESS_ESR_n50_n_100_J3_case2.pdf") #5,7
pdf("ESS_ESR_n50_n_100_J5_case1.pdf") #2,4
#pdf("ESS_ESR_n50_n_100_J5_case2.pdf") #6,8
#pdf("ESS_ESR_n50_n_100_J3_case6.pdf") #13,15
#pdf("ESS_ESR_n50_n_100_J5_case6.pdf") #14,16
l1<-layout(matrix( 1:6, nrow=2, byrow=T), widths=c( 3,3,2))
#layout.show(l1)

icounter<-0
for (idata in c(17,18)){
  icounter<-icounter+1
  
  if (icounter<2){par(mar=c(0,0,1,0), oma=c(6,4,0,1), new=F)
  }else{
    par(mar=c(0,0,0,0), new=F)
  }
  filename<-paste(Data_sets[idata],".RData", sep="")
  load(filename)
  
  if (idata==1) { output = output1 }
  if (idata==2) { output = output2 }
  if (idata==3) { output = output3 }
  if (idata==4) { output = output4 }
  if (idata==5) { output = output5 }
  if (idata==6) { output = output6 }
  if (idata==7) { output = output7 }
  if (idata==8) { output = output8 }
  
  if (idata==9) { output = output13 }
  if (idata==10) { output = output14 }
  if (idata==11) { output = output15 }
  if (idata==12) { output = output16 }
  if (idata==13) { output = output17 }
  if (idata==14) { output = output18 }
  if (idata==15) { output = output19 }
  if (idata==16) { output = output20 }
  
  if (idata==17) { output = output10 }
  if (idata==18) { output = output11 }
  
  outputESS<-output
  outputESR<-output
  labs<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  #ESS
  i1<-c(1,5,9,13)
  outputESS<-outputESS[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESS)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  boxplot(outputESS, col=c(2:5), outline=F, 
          names=NULL, cex.axis=.65, horizontal=T, ylim=c(0,12000),
          xaxt="n", yaxt="n", medlty=rep(1,4), medlwd=.5)
  
  if (icounter>1){
    axis(1, at=c(0, 2500, 5000, 7500, 10000),labels= c(0, 2500, 5000, 7500, 10000), las=2, cex=.8)
    
    mtext(expression(bold("Expected Sample size")), outer=TRUE, line=4, side=1, 
               col="black", adj=0.13, font=2, cex=.8, lwd=5)
    
    mtext(expression(bold("Expected Sampling rate")), outer=TRUE, line=4, side=1, 
          col="black", adj=0.58, font=2, cex=.8, lwd=2)
    
    mtext(expression(bold("Efficiency")), outer=TRUE, line=4, side=1, 
          col="black", adj=0.9, font=2, cex=.8)
  }
  abline(h=4.5, lty=3)
  abline(h=8.5, lty=3)
  abline(h=12.5, lty=3)
  
  mtext(expression("(n=50," ~ psi ~ "=0.3, J=5)"), outer=TRUE, line=2.5, side=2, 
        col="black", adj=0.75, font=2, cex=.7)
  
  mtext(expression(~ alpha[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.58, font=2, cex=.85)
  mtext(expression(~ alpha[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.68, font=2, cex=.85)
  mtext(expression(~ beta[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.79, font=2, cex=.85)
  mtext(expression(~ beta[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.9, font=2, cex=.85)
  
  #ESR
  if (icounter<2){par(mar=c(0,0,1,1.25), new=F)
  }else{
    par(mar=c(0,0,0,1.25), new=F)
    }
  outputESR[,1:4]<-outputESR[,1:4]/outputESR[,17]
  outputESR[,5:8]<-outputESR[,5:8]/outputESR[,18]
  outputESR[,9:12]<-outputESR[,9:12]/outputESR[,19]
  outputESR[,13:16]<-outputESR[,13:16]/outputESR[,20]
  outputESR<-outputESR[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESR)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)

  boxplot(outputESR, col=c(2:5), outline=F, 
          names=NULL, cex.axis=.65, horizontal=T, ylim=c(0,6000),
          xaxt="n", yaxt="n", medlty=rep(1,4), medlwd=.5)
  if (icounter>1){axis(1, at=c(0, 1500, 3000, 4500, 6000),labels= c(0, 1500, 3000, 4500, 6000), las=2, cex=.8)}
  abline(h=4.5, lty=3)
  abline(h=8.5, lty=3)
  abline(h=12.5, lty=3)
  
  mtext(expression("(n=100," ~ psi ~ "=0.3, J=5)"), outer=TRUE, line=2.5, side=2, 
        col="black", adj=0.2, font=2, cex=.7)
  
  mtext(expression(~ alpha[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.05, font=2, cex=.85)
  mtext(expression(~ alpha[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.18, font=2, cex=.85)
  mtext(expression(~ beta[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.29, font=2, cex=.85)
  mtext(expression(~ beta[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.42, font=2, cex=.85)
  
  #Timer
  boxplot(output[,c(17:20)]/output[,c(19)], names=c("JAGS", "dRUM", "PG", "RStan"), outline=F,
           cex.axis=.65, horizontal=T, ylim=c(0,30),
          xaxt="n", yaxt="n", col=c(2:5), medlty=rep(1,4), medlwd=.5)
  
  axis(4, at=1:4,labels= c("JAGS", "dRUM", "PG", "Stan")
         , las=3, cex=.5)
  
  for(i in 2:5){
    axis(side=4, at=i-1, col.axis=i, labels= c(expression(bold("JAGS")), expression(bold("dRUM")), expression(bold("PG")), expression(bold("Stan")))[i-1] , las=3)
    
  }
  
  if (icounter>1){axis(1, at=c(0, 15, 30),labels= c(0, 15, 30), las=3, cex=.8)}

  rm(output)
}
dev.off()

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
#Produce a plot of the ESS, ESR and Timing stats for the different methods
#for the big data set!
#----------------------------------------------------------------------------------

pdf("Logit_simulations_n500_J5_case3.pdf") #17
  l2<-layout(matrix( 1:3, nrow=1, byrow=T), widths=c( 3,3,2))
  idata<-17
  icounter<-1
  
  if (icounter<2){par(mar=c(0,0,1,0), oma=c(6,4,0,1), new=F)
  }else{
    par(mar=c(0,0,0,0), new=F)
  }
  filename<-paste(Data_sets[idata],".RData", sep="")
  load(filename)
  
  if (idata==17) { output = output9 }
  
  outputESS<-output
  outputESR<-output
  labs<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  #ESS
  i1<-c(1,5,9,13)
  outputESS<-outputESS[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESS)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  boxplot(outputESS, col=c(2:5), outline=F, 
          names=NULL, cex.axis=.65, horizontal=T, ylim=c(0,12000),
          xaxt="n", yaxt="n", medlty=rep(1,4), medlwd=.5)
  
    axis(1, at=c(0, 3000, 6000, 9000, 12000),labels= c(0, 3000, 6000, 9000, 12000), las=2, cex=1.5)
    
    mtext(expression(bold("Expected Sample size")), outer=TRUE, line=4, side=1, 
          col="black", adj=0.13, font=2, cex=.8, lwd=5)
    
    mtext(expression(bold("Expected Sampling rate")), outer=TRUE, line=4, side=1, 
          col="black", adj=0.58, font=2, cex=.8, lwd=2)
    
    mtext(expression(bold("Run-time\nEfficiency")), outer=TRUE, line=4, side=1, 
          col="black", adj=0.9, font=2, cex=.8)
    
  abline(h=4.5, lty=3)
  abline(h=8.5, lty=3)
  abline(h=12.5, lty=3)
  
  mtext(expression("(n=500," ~ psi ~ "=0.3, J=10)"), outer=TRUE, line=2.5, side=2, 
        col="black", adj=0.5, font=2, cex=1)
  
  mtext(expression(~ alpha[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.1, font=2, cex=.85)
  mtext(expression(~ alpha[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.36, font=2, cex=.85)
  mtext(expression(~ beta[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.61, font=2, cex=.85)
  mtext(expression(~ beta[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.87, font=2, cex=.85)
  
  #ESR
  if (icounter<2){par(mar=c(0,0,1,1.25), new=F)
  }else{
    par(mar=c(0,0,0,1.25), new=F)
  }
  outputESR[,1:4]<-outputESR[,1:4]/outputESR[,17]
  outputESR[,5:8]<-outputESR[,5:8]/outputESR[,18]
  outputESR[,9:12]<-outputESR[,9:12]/outputESR[,19]
  outputESR[,13:16]<-outputESR[,13:16]/outputESR[,20]
  outputESR<-outputESR[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESR)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  boxplot(outputESR, col=c(2:5), outline=F, 
          names=NULL, cex.axis=.65, horizontal=T, ylim=c(0,500),
          xaxt="n", yaxt="n", medlty=rep(1,4), medlwd=.5)
  axis(1, at=c(0, 125, 250, 375, 500),labels= c(0, 125, 250, 375, 500), las=2, cex=1.5)
  abline(h=4.5, lty=3)
  abline(h=8.5, lty=3)
  abline(h=12.5, lty=3)
  
  #Timer
  boxplot((output[,c(17:20)])/(output[,19]), names=c("JAGS", "dRUM", "PG", "RStan"), outline=F,
          cex.axis=1.5, horizontal=T, 
           yaxt="n", col=c(2:5),
          medcol=1, medlty=rep(1,4), medlwd=.5)
  
  axis(4, at=1:4,labels= 
         c("JAGS", "dRUM", "PG", "Stan"), las=3, cex=1.5)
  
  for(i in 2:5){
    axis(side=4, at=i-1, col.axis=i, labels= c(expression(bold("JAGS")), expression(bold("dRUM")), expression(bold("PG")), expression(bold("Stan")))[i-1] , las=3)
    
  }
  rm(output)
dev.off()

#idata =18
pdf("Logit_simulations_n500_J10_case3.pdf") #18
  l2<-layout(matrix( 1:3, nrow=1, byrow=T), widths=c( 3,3,2))
  
  idata<-18
  icounter<-1
  
  if (icounter<2){par(mar=c(0,0,1,0), oma=c(6,4,0,1), new=F)
  }else{
    par(mar=c(0,0,0,0), new=F)
  }
  filename<-paste(Data_sets[idata],".RData", sep="")
  load(filename)
  
  if (idata==18) { output = output10 }
  
  outputESS<-output
  outputESR<-output
  labs<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  #ESS
  i1<-c(1,5,9,13)
  outputESS<-outputESS[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESS)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  boxplot(outputESS, col=c(2:5), outline=F, 
          names=NULL, cex.axis=.65, horizontal=T, ylim=c(0,13000),
          xaxt="n", yaxt="n", medlty=rep(1,4), medlwd=.5)
  
  axis(1, at=seq(0,13000, length=5),labels= seq(0,13000, length=5), las=2, cex=1.5)
  
  mtext(expression(bold("Expected Sample size")), outer=TRUE, line=4, side=1, 
        col="black", adj=0.13, font=2, cex=.8, lwd=5)
  
  mtext(expression(bold("Expected Sampling rate")), outer=TRUE, line=4, side=1, 
        col="black", adj=0.58, font=2, cex=.8, lwd=2)
  
  mtext(expression(bold("Run-time\nEfficiency")), outer=TRUE, line=4, side=1, 
        col="black", adj=0.9, font=2, cex=.8)
  
  abline(h=4.5, lty=3)
  abline(h=8.5, lty=3)
  abline(h=12.5, lty=3)
  
  mtext(expression("(n=500," ~ psi ~ "=0.5, J=5)"), outer=TRUE, line=2.5, side=2, 
        col="black", adj=0.5, font=2, cex=1)
  
  mtext(expression(~ alpha[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.1, font=2, cex=.85)
  mtext(expression(~ alpha[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.36, font=2, cex=.85)
  mtext(expression(~ beta[0] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.61, font=2, cex=.85)
  mtext(expression(~ beta[1] ), outer=TRUE, line=1, side=2, 
        col="black", adj=0.87, font=2, cex=.85)
  
  #ESR
  if (icounter<2){par(mar=c(0,0,1,1.25), new=F)
  }else{
    par(mar=c(0,0,0,1.25), new=F)
  }
  outputESR[,1:4]<-outputESR[,1:4]/outputESR[,17]
  outputESR[,5:8]<-outputESR[,5:8]/outputESR[,18]
  outputESR[,9:12]<-outputESR[,9:12]/outputESR[,19]
  outputESR[,13:16]<-outputESR[,13:16]/outputESR[,20]
  outputESR<-outputESR[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESR)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  boxplot(outputESR, col=c(2:5), outline=F, 
          names=NULL, cex.axis=.65, horizontal=T, ylim=c(0,600),
          xaxt="n", yaxt="n", medlty=rep(1,4), medlwd=.5)
  axis(1, at=seq(0,600, length=5),labels= seq(0,600, length=5), las=2, cex=1.5)
  abline(h=4.5, lty=3)
  abline(h=8.5, lty=3)
  abline(h=12.5, lty=3)
  
  #Timer
  boxplot((output[,c(17:20)])/(output[,19]), names=c("JAGS", "dRUM", "PG", "RStan"), outline=F,
          cex.axis=1.5, horizontal=T, 
          yaxt="n", col=c(2:5),
          medcol=1, medlty=rep(1,4), medlwd=.5)
  
  axis(4, at=1:4,labels= 
         c("JAGS", "dRUM", "PG", "Stan"), las=3, cex=1.5)
  
  for(i in 2:5){
    axis(side=4, at=i-1, col.axis=i, labels= c(expression(bold("JAGS")), expression(bold("dRUM")), expression(bold("PG")), expression(bold("Stan")))[i-1] , las=3)
    
  }
  rm(output)
dev.off()
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#produce the plots that contains the median ESS and ESR
#detection probability approximately 0.5
pdf("ESSRESS1.pdf")
l3<-layout(matrix(c(1,2,17, 9,10,3,4,18, 11,12,5,6, 19, 13,14,7,8, 20,15,16), nrow=4, byrow=T), widths=c(1,1,.75,1,1))

par(mar=c(0,0,.5,0), oma=c(5,6,4,1))
par(tcl=-0.25)
par(mgp=c(2,.6,0))

#Calculation of the ESS
for (idata in 1:6){
  filename<-paste(Data_sets[idata],".RData", sep="")
  load(filename)
  
  if (idata==1) { output = output1 }
  if (idata==2) { output = output2 }
  if (idata==3) { output = output3 }
  if (idata==4) { output = output4 }
  if (idata==5) { output = output5 }
  if (idata==6) { output = output6 }
  #if (idata==7) { output = output7 }
  #if (idata==8) { output = output8 }
  
  working<-output
  
  medians<-median.stats(working)
  
  if ( (idata ==1) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white","white", "grey","white"), 
            ylim=c(0,10000), axes=F, medlwd=.5)
    
    axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)
  }
  if ( (idata ==2) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,10000), axes=F, medlwd=.5)
    
    #axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)    
  }
  if ( (idata ==3) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,10000), axes=F, medlwd=.5)
    
    axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)
  }
  if ( (idata ==4) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,10100), axes=F, medlwd=.5)
    
    #axis(2, at=c(0, 5200/2, 5200),labels= c(0, 5200/2, 5200), las=2)
  }
  if ( (idata ==5) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,12000), axes=F, medlwd=.5)
    
    axis(2, at=c(0, 6000, 12000),labels= c(0, 6000, 12000), las=2)
  }
  if ( (idata ==6) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,12000), axes=F, medlwd=.5)
    
    #axis(2, at=c(0, 5200/2, 5200),labels= c(0, 5200/2, 5200), las=2)
  }

  box(col="black")
}#end idata

#output7 used
filename<-paste(Data_sets[7],".RData", sep="")
load(filename)
output = output7

working<-output

medians<-median.stats(working)

boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,13000), axes=F, medlwd=.5)
  
axis(2, at=c(0, 6500, 13000),labels= c(0, 6500, 13000), las=2)
axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
box(col="black")
rm(output)


#output8 used
filename<-paste(Data_sets[8],".RData", sep="")
load(filename)
output = output8

working<-output

medians<-median.stats(working)

boxplot(medians, xaxt="n", 
        outline=F, col=c("white","white", "grey","white"), 
        ylim=c(0,13000), axes=F, medlwd=.5)
axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)

box(col="black")
rm(output)

mtext("Median ESS", outer=TRUE, line=2.5, side=2.5, col="black", adj=0.94, font=2, cex=.8)
mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.64, font=2, cex=.8)
mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.35, font=2, cex=.8)
mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.05, font=2, cex=.8)

mtext(expression("(n=50," ~ psi ~ "=0.3)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.935, font=2, cex=.7)

mtext(expression("(n=100," ~ psi ~ "=0.3)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.635, font=2, cex=.7)

mtext(expression("(n=50,"  ~ psi ~ "=0.5)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.345, font=2, cex=.7)

mtext(expression("(n=100," ~ psi ~ "=0.5)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.05, font=2, cex=.7)

mtext(expression("(J=3)"), outer=TRUE, line=1, side=3, col="black", adj=.1, font=2, cex=0.7)
mtext(expression("(J=5)"), outer=TRUE, line=1, side=3, col="black", adj=.3, font=2, cex=0.7)


#Calculation of the ESR

for (idata in 1:6){
  filename<-paste(Data_sets[idata],".RData", sep="")
  load(filename)
  
  if (idata==1) { output = output1 }
  if (idata==2) { output = output2 }
  if (idata==3) { output = output3 }
  if (idata==4) { output = output4 }
  if (idata==5) { output = output5 }
  if (idata==6) { output = output6 }
  if (idata==7) { output = output7 }
  if (idata==8) { output = output8 }
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  if ( (idata ==1) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white","white", "grey","white"), 
            ylim=c(0,4000), axes=F, medlwd=.5)
    
    axis(2, at=c(0, 2000, 4000),labels= c(0, 2000, 4000), las=2)
  }
  if ( (idata ==2) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,4000), axes=F, medlwd=.5)
    
    #axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)    
  }
  if ( (idata ==3) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,2100), axes=F, medlwd=.5)
    
    axis(2, at=c(0, 1000, 2000),labels= c(0, 1000, 2000), las=2)
  }
  if ( (idata ==4) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,2100), axes=F, medlwd=.5)
    
    #axis(2, at=c(0, 5200/2, 5200),labels= c(0, 5200/2, 5200), las=2)
  }
  if ( (idata ==5) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,4200), axes=F, medlwd=.5)
    
    axis(2, at=c(0, 2000, 4000),labels= c(0, 2000, 4000), las=2)
  }
  if ( (idata ==6) ){
    boxplot(medians, xaxt="n", 
            outline=F, col=c("white", "white", "grey","white"), 
            ylim=c(0,4200), axes=F, medlwd=.5)
    
    #axis(2, at=c(0, 5200/2, 5200),labels= c(0, 5200/2, 5200), las=2)
  }
  
  box(col="black")
}#end idata

  
  #output7 used
  filename<-paste(Data_sets[7],".RData", sep="")
  load(filename)
  output = output7
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,2300), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 1000, 2000),labels= c(0, 1000, 2000), las=2)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  #output8 used
  filename<-paste(Data_sets[8],".RData", sep="")
  load(filename)
  output = output8
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,2300), axes=F, medlwd=.5)
  
  box(col="black")
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  rm(output)
  
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.94, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.64, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.35, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.05, font=2, cex=.8)
  
  #for (i in 1:50){mtext(paste(2*i/100), outer=TRUE, line=5, side=2.5, col="black", adj=2*(i/100), font=2, cex=.3)}
  
  mtext(expression("(n=50," ~ psi ~ "=0.3)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.935, font=2, cex=.7)
  mtext(expression("(n=100," ~ psi ~ "=0.3)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.635, font=2, cex=.7)
  mtext(expression("(n=50,"  ~ psi ~ "=0.5)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.345, font=2, cex=.7)
  mtext(expression("(n=100," ~ psi ~ "=0.5)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.05, font=2, cex=.7)
  
  mtext(expression("(J=3)"), outer=TRUE, line=1, side=3, col="black", adj=.7, font=2, cex=0.7)
  mtext(expression("(J=5)"), outer=TRUE, line=1, side=3, col="black", adj=.91, font=2, cex=0.7)

dev.off()
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#produce the plots that contains the median ESS and ESR
#detection probability approximately 0.7
pdf("ESSRESS2.pdf")
  l3<-layout(matrix(c(1,2,17, 9,10,3,4,18, 11,12,5,6, 19, 13,14,7,8, 20,15,16), nrow=4, byrow=T), widths=c(1,1,.75,1,1))
  
  par(mar=c(0,0,.5,0), oma=c(5,6,4,1))
  par(tcl=-0.25)
  par(mgp=c(2,.6,0))
  
  #Calculation of the ESS
  for (idata in 9:14){
    filename<-paste(Data_sets[idata],".RData", sep="")
    load(filename)
    
    if (idata==9) { output = output13 }
    if (idata==10) { output = output14 }
    if (idata==11) { output = output15 }
    if (idata==12) { output = output16 }
    if (idata==13) { output = output17 }
    if (idata==14) { output = output18 }
    #if (idata==15) { output = output19 }
    #if (idata==16) { output = output20 }
    
    working<-output
    
    medians<-median.stats(working)
    
    if ( (idata ==9) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white","white", "grey","white"), 
              ylim=c(0,8500), axes=F, medlwd=.5)
      
      axis(2, at=c(0, 8500/2, 8500),labels= c(0, 8500/2, 8500), las=2)
    }
    if ( (idata ==10) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,8500), axes=F, medlwd=.5)
      
      #axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)    
    }
    if ( (idata ==11) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,8600), axes=F, medlwd=.5)
      
      axis(2, at=c(0, 8600/2, 8600),labels= c(0, 8600/2, 8600), las=2)
    }
    if ( (idata ==12) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,8600), axes=F, medlwd=.5)
    }
    if ( (idata ==13) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,10000), axes=F, medlwd=.5)
      
      axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)
    }
    if ( (idata ==14) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,10000), axes=F, medlwd=.5)
    }
    
    box(col="black")
    
  }#end idata
  
  #output19 used
  filename<-paste(Data_sets[15],".RData", sep="")
  load(filename)
  output = output19
  
  working<-output
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,10200), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 10200/2, 10200),labels= c(0, 10200/2, 10200), las=2)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  box(col="black")
  rm(output)
  
  
  #output20 used
  filename<-paste(Data_sets[16],".RData", sep="")
  load(filename)
  output = output20
  
  working<-output
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,10200), axes=F, medlwd=.5)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  mtext("Median ESS", outer=TRUE, line=2.5, side=2.5, col="black", adj=0.94, font=2, cex=.8)
  mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.64, font=2, cex=.8)
  mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.35, font=2, cex=.8)
  mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.05, font=2, cex=.8)
  
  mtext(expression("(n=50," ~ psi ~ "=0.3)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.935, font=2, cex=.7)
  
  mtext(expression("(n=100," ~ psi ~ "=0.3)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.635, font=2, cex=.7)
  
  mtext(expression("(n=50,"  ~ psi ~ "=0.5)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.345, font=2, cex=.7)
  
  mtext(expression("(n=100," ~ psi ~ "=0.5)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.05, font=2, cex=.7)
  
  mtext(expression("(J=3)"), outer=TRUE, line=1, side=3, col="black", adj=.1, font=2, cex=0.7)
  mtext(expression("(J=5)"), outer=TRUE, line=1, side=3, col="black", adj=.3, font=2, cex=0.7)
  
  for (idata in 9:14){
    filename<-paste(Data_sets[idata],".RData", sep="")
    load(filename)
    
    if (idata==9) { output = output13 }
    if (idata==10) { output = output14 }
    if (idata==11) { output = output15 }
    if (idata==12) { output = output16 }
    if (idata==13) { output = output17 }
    if (idata==14) { output = output18 }
    
    working<-output
    working[,1:4]<-working[,1:4]/working[,17]
    working[,5:8]<-working[,5:8]/working[,18]
    working[,9:12]<-working[,9:12]/working[,19]
    working[,13:16]<-working[,13:16]/working[,20]
    
    medians<-median.stats(working)
    
    if ( (idata ==9) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white","white", "grey","white"), 
              ylim=c(0,4000), axes=F, medlwd=.5)
      
      axis(2, at=c(0, 2000, 4000),labels= c(0, 2000, 4000), las=2)
    }
    if ( (idata ==10) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,4000), axes=F, medlwd=.5)
      
      #axis(2, at=c(0, 5000, 10000),labels= c(0, 5000, 10000), las=2)    
    }
    if ( (idata ==11) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,2100), axes=F, medlwd=.5)
      
      axis(2, at=c(0, 1000, 2000),labels= c(0, 1000, 2000), las=2)
    }
    if ( (idata ==12) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,2100), axes=F, medlwd=.5)
      
      #axis(2, at=c(0, 5200/2, 5200),labels= c(0, 5200/2, 5200), las=2)
    }
    if ( (idata ==13) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,4000), axes=F, medlwd=.5)
      
      axis(2, at=c(0, 2000, 4000),labels= c(0, 2000, 4000), las=2)
    }
    if ( (idata ==14) ){
      boxplot(medians, xaxt="n", 
              outline=F, col=c("white", "white", "grey","white"), 
              ylim=c(0,4000), axes=F, medlwd=.5)
      
      #axis(2, at=c(0, 5200/2, 5200),labels= c(0, 5200/2, 5200), las=2)
    }
    
    box(col="black")
    rm(output)
  }#end idata
  
  
  #output19 used
  filename<-paste(Data_sets[15],".RData", sep="")
  load(filename)
  output = output19
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,2300), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 1000, 2000),labels= c(0, 1000, 2000), las=2)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  #output20 used
  filename<-paste(Data_sets[16],".RData", sep="")
  load(filename)
  output = output20
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,2300), axes=F, medlwd=.5)
  
  box(col="black")
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  rm(output)
  
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.94, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.64, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.35, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.25, side=2, col="black", adj=0.05, font=2, cex=.8)
  
  #for (i in 1:50){mtext(paste(2*i/100), outer=TRUE, line=5, side=2.5, col="black", adj=2*(i/100), font=2, cex=.3)}
  
  mtext(expression("(n=50," ~ psi ~ "=0.3)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.935, font=2, cex=.7)
  mtext(expression("(n=100," ~ psi ~ "=0.3)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.635, font=2, cex=.7)
  mtext(expression("(n=50,"  ~ psi ~ "=0.5)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.345, font=2, cex=.7)
  mtext(expression("(n=100," ~ psi ~ "=0.5)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.05, font=2, cex=.7)
  
  mtext(expression("(J=3)"), outer=TRUE, line=1, side=3, col="black", adj=.7, font=2, cex=0.7)
  mtext(expression("(J=5)"), outer=TRUE, line=1, side=3, col="black", adj=.91, font=2, cex=0.7)

dev.off()

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#produce the plots that contains the median ESS and ESR
#for the larger data sets

pdf("ESSRESS3.pdf")
  l1<-layout(matrix(c(1,2, 9 ,5,6, 3, 4, 10, 7, 8), nrow=2, byrow=T), widths=c(1,1,.75,1,1))
  
  par(mar=c(0,0,.5,0), oma=c(5,6,4,1))
  par(tcl=-0.25)
  par(mgp=c(2,.6,0))
  
  #output9 used
  filename<-paste(Data_sets[17],".RData", sep="")
  load(filename)
  output = output9
  
  working<-output
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,10000), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 10000/2, 10200),labels= c(0, 10000/2, 10000), las=2)
  
  box(col="black")
  rm(output)
  
  #output10 used
  filename<-paste(Data_sets[18],".RData", sep="")
  load(filename)
  output = output10
  
  working<-output
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,10000), axes=F, medlwd=.5)
  #axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  #output11 used
  filename<-paste(Data_sets[19],".RData", sep="")
  load(filename)
  output = output11
  
  working<-output
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,15600), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 15000/2, 15600),labels= c(0, 15000/2, 15000), las=2)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  box(col="black")
  rm(output)
  
  #output12 used
  filename<-paste(Data_sets[20],".RData", sep="")
  load(filename)
  output = output12
  
  working<-output
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,15600), axes=F, medlwd=.5)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  mtext("Median ESS", outer=TRUE, line=2.5, side=2.5, col="black", adj=0.8, font=2, cex=.8)
  mtext("Median ESS", outer=TRUE, line=2.5, side=2, col="black", adj=0.2, font=2, cex=.8)
  
  mtext(expression("(n=500," ~ psi ~ "=0.3)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.8, font=2, cex=.7)
  mtext(expression("(n=500," ~ psi ~ "=0.5)"), outer=TRUE, line=4.5, side=2, col="black", adj=0.2, font=2, cex=.7)
  
  mtext(expression("(J=5)"), outer=TRUE, line=1, side=3, col="black", adj=.1, font=2, cex=0.7)
  mtext(expression("(J=10)"), outer=TRUE, line=1, side=3, col="black", adj=.3, font=2, cex=0.7)
  
  
  #===============================================================
  
  #Plot of the ESR histograms for the median values
  
  #output9 used
  filename<-paste(Data_sets[17],".RData", sep="")
  load(filename)
  output = output9
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,400), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 200, 400),labels= c(0, 200, 400), las=2)
  
  box(col="black")
  rm(output)
  
  #output10 used
  filename<-paste(Data_sets[18],".RData", sep="")
  load(filename)
  output = output10
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,400), axes=F, medlwd=.5)
  #axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  #output11 used
  filename<-paste(Data_sets[19],".RData", sep="")
  load(filename)
  output = output11
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,500), axes=F, medlwd=.5)
  
  axis(2, at=c(0, 250, 500),labels= c(0, 250, 500), las=2)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  box(col="black")
  rm(output)
  
  #output12 used
  filename<-paste(Data_sets[20],".RData", sep="")
  load(filename)
  output = output12
  
  working<-output
  working[,1:4]<-working[,1:4]/working[,17]
  working[,5:8]<-working[,5:8]/working[,18]
  working[,9:12]<-working[,9:12]/working[,19]
  working[,13:16]<-working[,13:16]/working[,20]
  
  medians<-median.stats(working)
  
  boxplot(medians, xaxt="n", 
          outline=F, col=c("white","white", "grey","white"), 
          ylim=c(0,500), axes=F, medlwd=.5)
  axis(1, at=c(1,2,3,4),labels= c("J","d","P", "S"), las=1)
  
  box(col="black")
  rm(output)
  
  mtext("Median ESR", outer=TRUE, line=-24.5, side=2, col="black", adj=0.8, font=2, cex=.8)
  mtext("Median ESR", outer=TRUE, line=-24.5, side=2, col="black", adj=0.2, font=2, cex=.8)
  
  mtext(expression("(n=500," ~ psi ~ "=0.3)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.8, font=2, cex=.7)
  mtext(expression("(n=500," ~ psi ~ "=0.5)"), outer=TRUE, line=-22.5, side=2, col="black", adj=0.2, font=2, cex=.7)
  
  mtext(expression("(J=5)"), outer=TRUE, line=1, side=3, col="black", adj=.7, font=2, cex=0.7)
  mtext(expression("(J=10)"), outer=TRUE, line=1, side=3, col="black", adj=.91, font=2, cex=0.7)

dev.off()
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------


#Tabulate the run-time efficiency of the methods
TimerEfficieny<- matrix(0, nrow=20, ncol=10)
AllTimers<-matrix(0, nrow=500, ncol=1)
medianESR<-matrix(0, nrow=20, ncol=16)

for (idata in 1:20){

  filename<-paste(Data_sets[idata],".RData", sep="")
  load(filename)
  
  if (idata==1) { output = output1 }
  if (idata==2) { output = output2 }
  if (idata==3) { output = output3 }
  if (idata==4) { output = output4 }
  if (idata==5) { output = output5 }
  if (idata==6) { output = output6 }
  if (idata==7) { output = output7 }
  if (idata==8) { output = output8 }
  
  if (idata==9) { output = output13 }
  if (idata==10) { output = output14 }
  if (idata==11) { output = output15 }
  if (idata==12) { output = output16 }
  if (idata==13) { output = output17 }
  if (idata==14) { output = output18 }
  if (idata==15) { output = output19 }
  if (idata==16) { output = output20 }
  
  if (idata==17) { output = output9 }
  if (idata==18) { output = output10 }
  if (idata==19) { output = output11 }
  if (idata==20) { output = output12 }

  #Timer<-na.omit((output[,c(17:20)])/output[,19])
  Timer<-(output[,c(17:20)])/output[,19]
  Timer<-Timer[,-3]
  
  AllTimers<-cbind(AllTimers, Timer)
  Timer<-na.omit(Timer)
  TimerEfficieny[idata,]<-c( apply(Timer,2,median), apply(Timer,2,mean), apply(Timer,2,sd), nrow(Timer))  
  
  #Obtain a comparison between the median effective sampling rate
  outputESS<-output
  outputESR<-output
  labs<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  #ESS
  i1<-c(1,5,9,13)
  outputESS<-outputESS[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESS)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  outputESR[,1:4]<-outputESR[,1:4]/outputESR[,17]
  outputESR[,5:8]<-outputESR[,5:8]/outputESR[,18]
  outputESR[,9:12]<-outputESR[,9:12]/outputESR[,19]
  outputESR[,13:16]<-outputESR[,13:16]/outputESR[,20]
  outputESR<-outputESR[,c(i1, i1+1, i1+2, i1+3)]
  colnames(outputESR)<-rep(c("JAGS", "dRUM", "PG", "RStan"),4)
  
  outputESRchange<-outputESR
  outputESRchange[,1:4]<-1/(outputESRchange[,1:4]/outputESRchange[,3])
  outputESRchange[,5:8]<-1/(outputESRchange[,5:8]/outputESRchange[,7])
  outputESRchange[,9:12]<-1/(outputESRchange[,9:12]/outputESRchange[,11])
  outputESRchange[,13:16]<-1/(outputESRchange[,13:16]/outputESRchange[,15])
  
  #outputESRchangem<-outputESRchange
  #outputESRchangem[,1:4]<-outputESRchangem[,c(1,)]
  
  cat("\n ", apply(na.omit(outputESRchange), 2, median))
  
  medianESR[idata,]<-apply(na.omit(outputESRchange), 2, median)
  rm(output)

}

matrix(apply(outputESRchange, 2, median), byrow=T, ncol=4)
matrix(apply(outputESRchange, 2, median), byrow=T, ncol=4)
matrix(apply(outputESRchange, 2, median), byrow=T, ncol=4)


colnames(TimerEfficieny)<-c("Jags-med", "drum-med","Stan-med",
                            "Jags-mean", "drum-mean","Stan-mean",
                            "Jags-sd", "drum-sd","Stan-sd",
                            "n")

TimerEfficieny<-TimerEfficieny[c( 
                                c(1,5,9,13),
                                c(2,6,10,14), 
                                c(3,7,11,15),
                                c(4,8,12,16),
                                c(17:20)), ]

AllTimers<-AllTimers[,-1]
AllTimers<-AllTimers[c( 
  c(1,5,9,13),
  c(2,6,10,14), 
  c(3,7,11,15),
  c(4,8,12,16),
  c(17:20)), ]

#pdf("Timer.pdf", paper='USr') 
#this picture was saved manually!
#first eps and then to pdf
  par(mfrow=c(1,1), mar=c(0,3,.5,0), oma=c(6,6,4,0.5), xpd=FALSE)
  boxplot(AllTimers, outline=F, 
          col=rep(c("red", "white", "grey"), 20), medlwd=.5, 
          names=NULL, cex.axis=.75, xaxt="n")
  abline(h=1, lty=2)
  abline(v=48.5, lty=2, lwd=2)
  abline(v=24.5, lty=2, lwd=2)
  abline(v=12.5, lty=2, lwd=1)
  abline(v=36.5, lty=2, lwd=1)
  par(xpd=NA)
  segments(1, 72, 24, 72, lwd=2)
  segments(25, 72, 48, 72, lwd=2)
  segments(49, 72, 60, 72, lwd=2)
  mtext("n=50", side=3, outer=TRUE, line=.25, col="black", adj=0.26, font=2, cex=.8, lwd=5)
  mtext("n=100", side=3, outer=TRUE, line=.25, col="black", adj=0.62, font=2, cex=.8, lwd=5)
  mtext("n=500", side=3, outer=TRUE, line=.25, col="black", adj=0.9, font=2, cex=.8, lwd=5)
  axis(2, at=35, labels="Timing-efficiency", line =2, tick=F)
  
  par(xpd=NA)
  for (i in 1:20){
    segments(1+3*(i-1), -2, 3+3*(i-1), -2, lwd=2, xpd=NA)
  }
  
  index<-seq(from=2, by=3, length=20)
  psivals<-c(rep(c(0.3, 0.5, 0.3, 0.5  ),4), 0.3, 0.3, 0.5, 0.5)
  for (i in 1:20){
    pv<-psivals[i]
    axis(1, at=index[i],labels=bquote(psi~"="~.(pv)), tick=F, line=-.95, cex.axis=.6)
  }
 
  for (i in 1:9){
    segments(1+6*(i-1), -5, 6+6*(i-1), -5, lwd=2, xpd=NA)
  }
  segments(49, -5, 60, -5, lwd=2, xpd=NA)
  
  index<-c(seq(from=3.5, by=6, length=8), 54.5)
  psivals<-c(rep(c(0.5, 0.7 ),4), 0.5)
  for (i in 1:9){
    pv<-psivals[i]
    axis(1, at=index[i],labels=bquote("p="~.(pv)), tick=F, line=0.25, cex.axis=.75)
  }

  rect(1, 66, 2, 68, col = "red", border="black")
  rect(1, 67, 2, 67, col = "red", border="black")
  text(x=3.5, y=67, "J", cex=.7)
  
  rect(1, 63, 2, 65, col = "white", border="black")
  rect(1, 64, 2, 64, col = "white", border="black")
  text(x=3.5, y=64, "d", cex=.7)
  
  rect(1, 60, 2, 62, col = "grey", border="black")
  rect(1, 61, 2, 61, col = "grey", border="black")
  text(x=3.5, y=61, "S", cex=.7)
#dev.off()

pdf("TimeEfficiency.pdf")
  par(mfrow=c(2,1),mar=c(0,1,.5,0), oma=c(6,6,4,1))
  
  plot(1:20, seq(0, 65, length=20), type="n", xaxt="n")
  axis(2, at=35, labels="Mean Timing-efficiency", line =2, tick=F)
  text(TimerEfficieny[,4], labels=rep(c("J3","J5"), 10), col=1, cex=.6)
  text(TimerEfficieny[,5], labels=rep(c("d3","d5"), 10), col=2, cex=.6)
  text(TimerEfficieny[,6], labels=rep(c("S3","S5"), 10), col="blue", cex=.6)
  abline(v=8.5, lty=2, lwd=2)
  abline(v=16.5, lty=2, lwd=2)
  abline(v=4.5, lty=2)
  abline(v=12.5, lty=2)
  segments(1, 70, 8, 70, lwd=2, xpd=TRUE)
  segments(9, 70, 16, 70, lwd=2, xpd=TRUE)
  segments(17, 70, 20, 70, lwd=2, xpd=TRUE)
  
  plot(1:20, seq(0, 8, length=20), type="n", xaxt="n")
  axis(2, at=4, labels="Std of Timing-efficiency", line =2, tick=F)
  text(TimerEfficieny[,7], labels=rep(c("J3","J5"), 10), col=1, cex=.6)
  text(TimerEfficieny[,8], labels=rep(c("d3","d5"), 10), col=2, cex=.6)
  text(TimerEfficieny[,9], labels=rep(c("S3","S5"), 10), col="blue", cex=.6)
  abline(v=8.5, lty=2, lwd=2)
  abline(v=16.5, lty=2, lwd=2)
  abline(v=4.5, lty=2)
  abline(v=12.5, lty=2)
  
  mtext("n=50", side=3, outer=TRUE, line=.25, col="black", adj=0.22, font=2, cex=.8, lwd=5)
  mtext("n=100", side=3, outer=TRUE, line=.25, col="black", adj=0.62, font=2, cex=.8, lwd=5)
  mtext("n=500", side=3, outer=TRUE, line=.25, col="black", adj=0.92, font=2, cex=.8, lwd=5)
  
  axis(1, at=1:20, 
       labels=c(rep(c("0.3","0.5", "0.3","0.5"),4), "0.3", "0.3", "0.5","0.5"), 
       tick=F, line=0.25, cex.axis=.7)
  
  mtext(expression(psi), side=1, outer=TRUE, line=1.35, col="black", adj=0.025, 
        font=2, las=1)
  
  par(xpd=NA)
  segments(1, -2, 2, -2, lwd=2, xpd=NA)
  segments(3, -2, 4, -2, lwd=2, xpd=NA)
  
  segments(5, -2, 6, -2, lwd=2, xpd=NA)
  segments(7, -2, 8, -2, lwd=2, xpd=NA)
  
  segments(9, -2, 10, -2, lwd=2, xpd=NA)
  segments(11, -2, 12, -2, lwd=2, xpd=NA)
  
  segments(13, -2, 14, -2, lwd=2, xpd=NA)
  segments(15, -2, 16, -2, lwd=2, xpd=NA)
  
  segments(17, -2, 20, -2, lwd=2, xpd=NA)
  
  mtext("p=0.5", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.073, font=2, cex=.65, lwd=5)
  
  mtext("p=0.7", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.173, font=2, cex=.65, lwd=5)
  
  mtext("p=0.5", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.271, font=2, cex=.65, lwd=5)
  
  mtext("p=0.7", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.37, font=2, cex=.65, lwd=5)
  
  mtext("p=0.5", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.468, font=2, cex=.65, lwd=5)
  
  mtext("p=0.7", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.568, font=2, cex=.65, lwd=5)
  
  mtext("p=0.5", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.665, font=2, cex=.65, lwd=5)
  
  mtext("p=0.7", side=1, outer=TRUE, line=2.2, 
        col="black", adj=0.765, font=2, cex=.65, lwd=5)
  
  mtext("p=0.5", side=1, outer=TRUE, line=2.2, 
        col="black", adj=.91, font=2, cex=.65, lwd=5)
dev.off()

