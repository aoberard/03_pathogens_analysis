#######################################################
#
#   Use Screening Analysis to Investigate Associations
#
#######################################################


SCN<-function(test,pathos){

########################################################
#     Libraries & Functions          
########################################################
########################################################

# pour le calcul d'enveloppe de confiance ‡ 95 %
library(boot)
# library(plotrix)
# library(lattice)
# Fonction utilisÈe :
#source("FctTestScreenENV.txt")

########################################################
#        Run the SCN              
########################################################
########################################################

  set.seed(100) 
  
  ############################
  # DonnÈes initiales :
  Xdat <- test
  Xtab<-test[,c(pathos)]
    # Struture de sauvegarde des rÈsultats :
    ResSCRENV <- FctTestScreenENV(Xtab)
    #
    # test significatif si > 0	
    ResSCRENV[[1]];	
    #
    # Pvalue			
    ResSCRENV[[2]];	
    #
    xNC <- length(ResSCRENV[[3]]);
    x <- (1:xNC);
    #     yinf <- ResSCRENV[[8]];
    #     ysup <- ResSCRENV[[7]];
    #     yobs <- ResSCRENV[[6]];
    #
    #     plot(x,yobs,type="n",ylim=range(c(0.0,yinf,ysup,yobs)),
    #          xlab="indice de combinaison",ylab="effectif",
    #          main="jeu de donnÈes BB/BR")
    #     lines(x,yinf,lty=1,col=4,lwd=2) #lower 95% 
    #     lines(x,ysup,lty=1,col=3,lwd=2) #upper 95%
    #     points(x,yobs,col=2) #observed
    valx <- x[(ResSCRENV[[4]]!=0)];
    #     for (i in 1:length(valx) ) {
    #          abline(v=valx[i],lty=3,col=1,lwd=1)
    #          }
    #
    # combinaison significatives des germes	 
    (ResSCRENV[[3]])[valx];	
    #
    # Pvalue des combinaisons significatives 
    (ResSCRENV[[5]])[valx];	
    #
    cbind(ResSCRENV[[3]],ResSCRENV[[4]],ResSCRENV[[5]]);
    # Combinaison - compris dans l'enveloppe (1) ou non (0) - P-value associÈe ‡ chaque combinaison
    cbind(ResSCRENV[[3]][valx],ResSCRENV[[6]][valx],ResSCRENV[[8]][valx],ResSCRENV[[7]][valx]);
    # Combinaison hors de l'enveloppe - obs - IC
    tt<-cbind(ResSCRENV[[3]],ResSCRENV[[6]],ResSCRENV[[8]],ResSCRENV[[7]],ResSCRENV[[4]],ResSCRENV[[5]]);
    # Combinaisons - obs - IC - significativitÈ - P-value
    tt;
    
########################################################
#       Define directionality of each association              
########################################################
########################################################

    # number of pathogens
    n <-dim(Xtab)[2]
    
    # number of combinations
    ncombos<-1/(0.5^(n))
    
direction<-rep(NA,ncombos)
tt<-as.data.frame(cbind(tt,direction))
for (qq in 1:ncombos){
  if(tt[qq,2]>ceiling((tt[qq,4]-tt[qq,3])/2)+tt[qq,3]){tt[qq,"direction"]<-as.character("Frequent")}else
    if(tt[qq,2]<floor((tt[qq,4]-tt[qq,3])/2)+tt[qq,3]){tt[qq,"direction"]<-"Rare"}else{tt[qq,"direction"]<-"Random"}
}

########################################################
#        Return Results              
########################################################
########################################################

names(tt)<-c(paste(names(Xtab),collapse="-"), "observed","lower bound", "upper bound","significant(Y/N)","pseudo-P value","direction")
return(tt)
}



