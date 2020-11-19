#############################################################################################
#############################################################################################
## CODE FOR REPONSE TO REVISION                                                            ##
## FOR THE PAPER: "TRANS-ACTING GENETIC VARIATION AFFECTS THE EXPRESSION OF ADJACENT GENES ##
## CODE WRITTEN BY: BY KRISNA VAN DYKE                                                     ##               
#############################################################################################
#############################################################################################

# README
# this code exists to address any concerns by our reviewers
# # it adds a new supplementary table and conducts some rudimentary analyses
# this code assumes the following
# # you have finished running the code provided in the main body
# # you can run this on a clean session assuming you did that, should run anyway though 
# all necessary libraries and funtions are included separately here 

# SET UP WORKING DIRECTORY
wd <-'SELECT_DIRECTORY'
setwd(wd)

# DECIDE IF YOU WANTS ANALYSES AND PLOTS DONE OR JUST MAKE OBJECTS
# you cannot make plots without doing analyses
# you cannot check sanity without making plots
# DECIDE IF YOU WANTS ANALYSES AND PLOTS DONE OR JUST MAKE OBJECTS
# you cannot make plots without doing analyses
# you cannot check sanity without making plots
do_analyses <- T
make_plots <- T
check_sanity <- T


# LOAD RELEVANT LIBRARIES
library(car) # for the Type III ANOVA
library(Hmisc) # various useful functions
library(magrittr) # for the ever useful pipe operator
library(MASS) # for negative binomial regression
library(preprocessCore) # used for efficient quantile normalization
library(Rcpp) # allows C++ code to be run in R
library(reshape2) # for melting the correlation matrix
library(rtracklayer) # for the import function and getting things as granges
library(tidyverse) # for some data processing

##################
# RCPP FUNCTIONS #
##################


cppFunction('int nCount(NumericVector v, IntegerVector key) {
  int m = key.size();
  int q = 0;
    
  for (int j = 0; j < m; j++){
  
    if (((v[key[j]] > 0) & (v[key[j]-1] > 0)) | ((v[key[j]] < 0) & (v[key[j]-1] < 0))){
      q++;
    }
  }
  
  return q;
}')


cppFunction('IntegerVector permute(int t, NumericVector hs, IntegerVector key) {
  srand (time(NULL));
  IntegerVector out(t);
  int m = key.size();
  int n = hs.size();

  for(int i = 0; i < t; i++) {
  
    // SHUFFLE THE VECTOR USING FISHER-YATES
    for (int j = 0; j < n; j++) {
      int p = rand() % n;
      float temp = hs[p];
      hs[p] = hs[j];
      hs[j] = temp;
      }
    
    // COUNT THE NUMBER OF NEIGHBORS      
    out[i] = 0;
    for (int k = 0; k < m; k++){
        if (((hs[key[k]] > 0) & (hs[key[k]-1] > 0)) | ((hs[key[k]] < 0) & (hs[key[k]-1] < 0))){
            out[i]++;
        }
    }
    
  }

  return out;
}')

###############
# R FUNCTIONS #
###############

# gets the adjacent-to-diagonal elements of a matrix
diagonal_grabber <- function(x){
  d <- c(rep(NA,ncol(x)))
  for (i in 2:length(d)){
    d[i] <- x[i,i-1]
  }
  return(d)
}

# retrieves neighbors and removes false neighbors
neighbor_grab <- function(x,key){
  diagonal_grabber(x)[-key]
}


# for distance relationships for each member of a vector
sqgen <- function(y){
  return(matrix(rep(y,length(y)),ncol=length(y),nrow=length(y))) 
}

# generates a distance matrix from a vector
distancer <- function(x){
  return(abs(sqgen(x) - t(sqgen(x))))
}

# shorthand which extracts and returns the upper triangle of a matrix
tri <- function(x){
  return(x[upper.tri(x,diag=F)])
}

sdize <- function(x){
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

#########################################
# LOAD OBJECTS PRODUCED BY  MAIN SCRIPT #
#########################################

# These all should have been generated already if you ran the main script 

load(paste0(wd,'/intermediary/b_mod.RData'))
load(paste0(wd,'/intermediary/metaNT.RData'))
load(paste0(wd,'/intermediary/chrkey.RData'))
load(paste0(wd,'/intermediary/bplus.RData'))
load(paste0(wd,'/intermediary/','jacTF.RData'))
load(paste0(wd,'/intermediary/baseWeiner.RData'))
load(paste0(wd,'/intermediary/deltaWeiner.RData'))
load(paste0(wd,'/intermediary/','ssWang.RData'))
load(paste0(wd,'/intermediary/','giSim.RData'))

# produce the pair orientation vector
pOrient <- rep(NA,nrow(metaNT))
for(i in 2:nrow(metaNT)){
  lead <- metaNT$strand[i] 
  lag <- metaNT$strand[i-1]
  if(lead < 0 & lag > 0){
    pOrient[i] <- 'c'
  }else if(lead > 0 & lag < 0){
    pOrient[i] <- 'd'
  } else {
    pOrient[i] <- 't'
  }
}

# prepare the basepair distance matrix
temp <- distancer(metaNT$chromosome) == 0
# make sure everything not on the same chromosome becomes NA
temp[!temp] <- NA
# get the values
intrachr_dist <- distancer(metaNT$start)*temp 
rm(temp)

#####################################################
# EXCESS BY WHICH REAL DOUBLETS EXCEED PERMUTATIONS #
#####################################################

# histogram for observed doublets minus mean doublets 
  # the red bar indicates a difference of 5
  # this is not to be included as a figure but is included in the respnse to our reviewers 
if(make_plots){
  temp <- read.csv(file=paste0(wd,'/tables/','Table_S1_01_doublets.csv'))
  
  excess <- temp$obsDoublets-temp$Mean
  
  # summary stats for those significantly in excess
  description <- excess[temp$pVal < 0.05] %>% summary
  
  # plot the excess versus the -log10(p)
  # # above the red line are significant entries at nominal 0.05
  # # the black vertical line is the median excess of those significantly in excess
  png(paste0(wd,'/plots/','revResponse_ExcessVsP.png'),width=800,height=600)
  plot(excess,-log10(temp$pVal),
       xlab = 'Excess of doublets over permutation mean',
       ylab= expression('-log'[10]*'(p-values)'))
  abline(h=-log10(0.05),lwd=3,col='red')
  abline(v=description[3])
  dev.off()
  
  # histogram of how many hotspots affect gene pairs
  png(paste0(wd,'/plots/','revResponse_HotspotsOnPairs.png'),width=800,height=600)
  hist(neighbor_grab(bplus,chrkey),
       breaks=20,
       main = 'Histogram of doublet counts per gene pair',
       xlab='Number of hotspots affecting gene pair')
  dev.off()
  
  # how many are in excess by greater than 5 or 20 (the reviewer bounds)
  sum(excess[temp$pVal < 0.05] > 5)/length((excess[temp$pVal < 0.05]))
  sum(excess[temp$pVal < 0.05] > 20)/length((excess[temp$pVal < 0.05]))
  
  # remove unnecessary objects
  rm(temp,excess,description)  
}

####################################################
# REAL DOUBLETS EXCEED PERMUTATIONS WITHOUT TOP 20 #
####################################################

# prepare indices for the creation of the permutation
numGenes <- 1:(nrow(b_mod)-1)
doublets <- diagonal_grabber(bplus)
top20 <- which(doublets %in% (sort(doublets,decreasing=T)[1:20] %>% unique()))
if(check_sanity){
 # since there are ties, we cut not just the top 20 pairs but the top 23
  sum(doublets %in% (sort(doublets,decreasing=T)[1:20] %>% unique()))
  # look at what these numbers are
  doublets[which(doublets %in% (sort(doublets,decreasing=T)[1:20] %>% unique()))] %>% sort()
}
n20 <- numGenes[!(numGenes %in% sort(unique(c(chrkey,top20))))]

# GENERATE PERMUTATIONS MINUS THE TOP 20 MOST AFFECTED GENE PAIRS (~18 min)
load(paste0(wd,'/intermediary/top20_permutations.RData'))
if(!exists('top20_permutations')){
  top20_permutations <- apply(b_mod,2,function(x)permute(1e5,x,n20))
  save(top20_permutations,file=paste0(wd,'/intermediary/top20_permutations.RData'))
}

if(do_analyses){
  
  # doublet counts for each hotspot after omission of the top 20 most affected gene pairs
  n20_nCounts <- apply(b_mod,2,function(x)nCount(x,n20))
  
  # p-values comparing real data doublets to permutation doublets
  # # any p-values appearing to be zero are actually p < number_permutations
  n20_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(top20_permutations),n20_nCounts)
  
  n20Perm <- as.data.frame(matrix(NA,nrow=1,ncol=6)) %>%
    `colnames<-`(c('Bonf_Count','Bonf_pVal','Nom_Count','Nom_pVal','Med_Count','Med_pVal')) %>%
    `rownames<-`(c('n20'))
  
  # binomial tests at different significance thresholds
  n20Perm['n20','Bonf_Count'] <- sum(n20_pVals < 0.05/length(n20_pVals))
  n20Perm['n20','Bonf_pVal'] <- binom.test(sum(n20_pVals < 0.05/length(n20_pVals)),
                                                 length(n20_pVals),0.05/length(n20_pVals),alternative='greater')$p.value
  n20Perm['n20','Nom_Count'] <- sum(n20_pVals < 0.05)
  n20Perm['n20','Nom_pVal'] <-binom.test(sum(n20_pVals < 0.05),length(n20_pVals),0.05,alternative='greater')$p.value
  n20Perm['n20','Med_Count'] <- sum(n20_pVals < 0.5)
  n20Perm['n20','Med_pVal'] <- binom.test(sum(n20_pVals < 0.5),length(n20_pVals),0.5,alternative='greater')$p.value
  
  if(make_plots){
    write.csv(n20Perm,file=paste0(wd,'/tables/','revResponse_n20Perm.csv'))
  }
  
  rm(n20_nCounts,n20_pVals,n20Perm)
  
}

rm(numGenes,doublets,top20,n20,top20_permutations)

#####################################
# THE MODEL WITH NORMALIZED WEIGHTS #
#####################################

# The final Type III ANOVA results are, as expected, unchanged 
if(do_analyses){
  
  # filter out the zeros you created with the additon of null rows 
  z <- neighbor_grab(bplus,chrkey) != 0
  
  # provide appropriate names to use in the model 
  Closeness <- log(1/neighbor_grab(intrachr_dist,chrkey))[z]
  doubletScore <- neighbor_grab(bplus,chrkey)[z]
  PairOrientation <- pOrient[-chrkey][z]
  TFInventory <- neighbor_grab(jacTF,chrkey)[z]
  BaselineChromatin <- baseWeiner[z]
  DeltaChromatin <- deltaWeiner[z]
  
  # functional similarity 
  geneInteractionSimilarity <- neighbor_grab(giSim,chrkey)[z]
  geneOntologySimilarity <- neighbor_grab(ssWang,chrkey)[z]
  
  # for the factorCorr matrix
  divergent <- PairOrientation == 'd'
  tandem <-  PairOrientation == 't'
  convergent <- PairOrientation == 'c'
  

  # get correlations between factors   
  factMat <- cbind(divergent, tandem, convergent, TFInventory,Closeness,BaselineChromatin,
                   DeltaChromatin,geneInteractionSimilarity,geneOntologySimilarity) 
  factCorr <- factMat %>% rcorr(type='spearman')
  
  # get the tables of these
  write.csv(as.data.frame(factCorr$r), file=paste0(wd,'/tables/','Table_S3_factCorrR.csv'))
  write.csv(as.data.frame(factCorr$P), file=paste0(wd,'/tables/','Table_S3_factCorrP.csv'))
  
  
  rm(factMat,factCorr,divergent,tandem,convergent)
  
  # provide exposure variables for both genes in each pair 
  exposure <- apply(b_mod,1,function(x)sum(x!=0))
  g1exposure <- c(NA,exposure[2:length(exposure)])[-chrkey][z]
  g2exposure <- c(NA,exposure[1:(length(exposure)-1)])[-chrkey][z]
  
  # apply(factMat[,3:ncol(factMat)],2,function(x)(x-mean(x,na.rm=T))/sd(x,na.rm=T)) %>% 
  
  fxNorm_model <- glm.nb(doubletScore ~ 
                          PairOrientation + 
                           sdize(TFInventory) +
                           sdize(Closeness) +
                           sdize(BaselineChromatin) +
                           sdize(DeltaChromatin) +
                           sdize(geneInteractionSimilarity) +
                           sdize(geneOntologySimilarity) +
                          offset(log(g1exposure)) + 
                          offset(log(g2exposure)),
                        maxit=1000)
  
  fxNorm_t3 <- fxNorm_model %>% 
    Anova(type='III') %>%
    `names<-`(c('LR_Chisq','Df','pVals'))
  
  
  # format it
  fxNorm_t3 <- fxNorm_t3 %>% add_column(Factor=rownames(fxNorm_t3)) %>%
    mutate(Significance=cut(pVals, breaks=c(-Inf,0.001,0.05,Inf), 
                            labels=c("p < 0.001","p < 0.05","p > 0.05")))
  
  write.csv(fxNorm_t3[-(which(names(fxNorm_t3) == 'Significance'))] %>% `rownames<-`(c(NULL)),
            file=paste0(wd,'/tables/','revResponse_FxNorm.csv'))
  
  
}

rm(list=ls())