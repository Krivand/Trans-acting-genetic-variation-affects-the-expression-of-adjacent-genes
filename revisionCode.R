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
library(GenomicRanges) # used for efficient accessing of basepair data
library(Hmisc) # various useful functions
library(magrittr) # for the ever useful pipe operator
library(MASS) # for negative binomial regression
library(preprocessCore) # used for efficient quantile normalization
library(Rcpp) # allows C++ code to be run in R
library(readxl) # read from MS excel
library(rtracklayer) # for the import function and getting things as granges
library(tidyverse) # for some data processing
library(viridis) # for manual coloring using viridis

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

###################
# C.U.T. ANALYSIS #
###################

if(do_analyses){
  
  # PREPARE THE CUTS AS A GRANGE
  cuts <- read_excel(paste0(wd,'/raw/',"12864_2016_2622_MOESM5_ESM.xlsx"), sheet = 1) %>% 
    as.data.frame() %>% rename(chromosome=Chromosome, end=Stop, start=Start) %>%
    mutate(chromosome = 
             recode_factor(chromosome,
                           chrI='1',chrII='2',chrIII='3',chrIV='4',
                           chrV='5',chrVI='6',chrVII='7',chrVIII='8',
                           chrIX = '9',chrX ='10',chrXI='11',chrXII='12',
                           chrXIII='13',chrXIV='14',chrXV='15',chrXVI='16')) %>% 
    mutate(chromosome=as.numeric(as.character(chromosome))) %>% makeGRangesFromDataFrame()
  
  # GET THE STRAND INVARIANT START AND END POSITIONS 
  temp <- metaNT[c('start','end','chromosome')]
  
  for(i in 1:nrow(temp)){
    if(temp$end[i] < temp$start[i]){
      flipper <- temp$start[i]
      temp$start[i] <- temp$end[i]
      temp$end[i] <- flipper
    }
  }
  
  rm(flipper)
  
  # PREPARE A GRANGE OF INTERGENIC SPACES
  # inelegant
  
  intergenic <- as.data.frame(matrix(NA,ncol=3,nrow=nrow(temp))) %>%
    `colnames<-`(c('start','end','chromosome'))
  
  for(i in 2:nrow(intergenic)){
    intergenic$start[i] <- temp$end[i-1]
    intergenic$end[i] <- temp$start[i]
    intergenic$chromosome[i] <- temp$chromosome[i]
  }
  
  intergenic <- intergenic[-chrkey,]
  
  # make sure overlapping genes don't get any overlaps but numbering is consistent
  for(i in 1:nrow(intergenic)){
    if(intergenic$end[i] < intergenic$start[i]){
      intergenic$start[i] <- 1
      intergenic$end[i] <- 1
      intergenic$chromosome[i] <- 1e10
    }
  }
  
  
  intergenic <- intergenic %>% makeGRangesFromDataFrame()
  
  # count the number of annotated CUTs in each intergenic region
  isect <- countOverlaps(intergenic, cuts)
  
  # examine how many CUTs are found per intergenic region 
  table(isect)
  
  # see if the presence of a cut changes pairing score
  wilcox.test(bplus[isect != 0],bplus[isect == 0]) %>% print()
  
  rm(intergenic,isect)
  
  # PREPARE A GRANGE OF BOTH GENES AND THEIR INTERGENIC SPACES
  # inelegant
  
  pairSpan <- as.data.frame(matrix(NA,ncol=3,nrow=nrow(temp))) %>%
    `colnames<-`(c('start','end','chromosome'))
  
  for(i in 2:nrow(pairSpan)){
    pairSpan$start[i] <- temp$start[i-1]
    pairSpan$end[i] <- temp$end[i]
    pairSpan$chromosome[i] <- temp$chromosome[i]
  }
  
  pairSpan <- pairSpan[-chrkey,]
  
  # make sure overlapping genes don't get any overlaps but numbering is consistent
  for(i in 1:nrow(pairSpan)){
    if(pairSpan$end[i] < pairSpan$start[i]){
      pairSpan$start[i] <- 1
      pairSpan$end[i] <- 1
      pairSpan$chromosome[i] <- 1e10
    }
  }
  
  
  pairSpan <- pairSpan %>% makeGRangesFromDataFrame()
  
  # count the number of annotated CUTs in each pair spanning region
  isect <- countOverlaps(pairSpan, cuts)
  
  # examine how many CUTs are found per pair spanning region 
  table(isect)
  
  # see if the presence of a cut changes pairing score
  wilcox.test(bplus[isect != 0],bplus[isect == 0]) %>% print()
  
  rm(temp,cuts,pairSpan,isect)
  
}


#####################################################
# EXCESS BY WHICH REAL DOUBLETS EXCEED PERMUTATIONS #
#####################################################

# histogram for observed doublets minus mean doublets 
  # the red bar indicates a difference of 5
  # this is not to be included as a figure but is included in the respnse to our reviewers 
if(make_plots){
  temp <- read.csv(file=paste0(wd,'/tables/','File_S2_01_doublets.csv'))
  
  excess <- temp$obsDoublets-temp$Median
  
  # num gene pairs affected by more hotspots than median of nominally significant excess
  sum(neighbor_grab(bplus,chrkey) >= 18)
  
  # how many exceed the 23 boundary set by the top 20
  sum(excess > 23)
  sum(excess > 23 & temp$pVal < 0.05)
  sum(excess <= 23)
  sum(excess <= 23 & temp$pVal < 0.05)
  
  # how many are in excess by greater than 5 or 20 (the reviewer bounds)
  sum(excess[temp$pVal < 0.05] > 5)/length((excess[temp$pVal < 0.05]))
  
  summary(excess[temp$pVal < 0.05])
  
  
  # 
  png(paste0(wd,'/plots/','revResponse_ExcessHist.png'),width=400,height=300)
  hist(excess[temp$pVal < 0.05],breaks=20,
       main='Histogram of excess doublets\n for nominally significant hotspots',
       xlab='Excess of doublets over permutation median',
       cex.lab=1.3,
       cex.axis=1.3)
  abline(v=5,lwd=3,col='red')
  abline(v=20,lwd=3,col='red')
  abline(v=median(excess[temp$pVal < 0.05]),lwd=3,col='blue')
  legend(40, 10, legend=c("5 & 20 bounds", "Median"),
         col=c("red", "blue"), lwd=3,cex=1.2,bty='n')
  dev.off()
  
  # plot the excess versus the -log10(p)
  # # above the red line are significant entries at nominal 0.05
  # # the black vertical line is the median excess of those significantly in excess
  png(paste0(wd,'/plots/','revResponse_ExcessVsP.png'),width=400,height=300)
  plot(excess,-log10(temp$pVal),
       main='Excess versus significance',
       xlab = 'Excess of doublets over permutation median',
       ylab= expression('-log'[10]*'(p-value)'),
       cex.lab=1.4,
       cex.axis=1.4)
  abline(h=-log10(0.05),lwd=3,col='red')
  abline(v=median(excess[temp$pVal < 0.05]),lwd=3,col='blue')
  legend(18, 1, legend=c("p = 0.05", "Median of nominally significant"),
         col=c("red", "blue"), lwd=3,cex=1.1,bty='n')
  dev.off()
  
  # the ggplot version
  templot <- cbind(excess,-log10(temp$pVal)) %>% as.data.frame() %>% 
    `names<-` (c('Excess','pVal')) 
  templot$pVal[templot$pVal == Inf] <- -log10(1/100000)
  templot <- templot %>% ggplot() +
    geom_vline(aes(xintercept=median(excess[temp$pVal < 0.05]), color = 'Median of nominally significant hotspots'), size=1) + 
    geom_hline(aes(yintercept=-log10(0.05), color = 'p < 0.05'), size=1) +
    geom_point(aes(x=Excess,y=pVal)) + theme_bw() + 
    ylab(expression('-log'[10]*'(p-value)')) + xlab('Excess of doublets over permutation median') + 
    scale_color_viridis(name="statistics",discrete=TRUE) +
    theme(legend.position="bottom",legend.title = element_blank(),
          axis.text.x = element_text(color='black'),
          axis.text.y = element_text(color='black'))

  ggsave(paste0(wd,'/plots/','Fig_S3_ExcessVsP.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 7.08,height = 7,units = c("in"),useDingbats=F)
  
  rm(templot)
  
  # histogram of how many hotspots affect gene pairs
  png(paste0(wd,'/plots/','revResponse_HotspotsOnPairs.png'),width=400,height=300)
  hist(neighbor_grab(bplus,chrkey),
       breaks=20,
       main = 'Histogram of doublet counts per gene pair',
       xlab='Number of hotspots affecting gene pair',
       cex.lab=1.4,
       cex.axis=1.4)
  dev.off()
  
  # remove unnecessary objects
  rm(temp,excess)  
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
  
  # get the permutations summaries for each hotspot 
  if(make_plots){
    apply(top20_permutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = apply(b_mod,2,function(x)sum(x!=0))) %>% add_column (obsDoublets = n20_nCounts) %>%
      add_column(pVal = n20_pVals) %>% write.csv(file=paste0(wd,'/tables/','File_S2_07_noTop20.csv'))
  }
  
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
    write.csv(n20Perm,file=paste0(wd,'/tables/','File_S2_08_noTop20Summary.csv'))
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
  Proximity <- log(1/neighbor_grab(intrachr_dist,chrkey))[z]
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
  factMat <- cbind(divergent, tandem, convergent, TFInventory,Proximity,BaselineChromatin,
                   DeltaChromatin,geneInteractionSimilarity,geneOntologySimilarity) 
  factCorr <- factMat %>% rcorr(type='spearman')
  
  # get the tables of these
  write.csv(as.data.frame(factCorr$r), file=paste0(wd,'/tables/','File_S3_factCorrR.csv'))
  write.csv(as.data.frame(factCorr$P), file=paste0(wd,'/tables/','File_S3_factCorrP.csv'))
  
  
  rm(factMat,factCorr,divergent,tandem,convergent)
  
  # provide exposure variables for both genes in each pair 
  exposure <- apply(b_mod,1,function(x)sum(x!=0))
  g1exposure <- c(NA,exposure[2:length(exposure)])[-chrkey][z]
  g2exposure <- c(NA,exposure[1:(length(exposure)-1)])[-chrkey][z]
  
  # build the model
  fxNorm_model <- glm.nb(doubletScore ~ 
                          PairOrientation + 
                           sdize(TFInventory) +
                           sdize(Proximity) +
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
  
  # remove unnecessary objects
  rm(z,Proximity,doubletScore,PairOrientation,TFInventory,BaselineChromatin,DeltaChromatin, 
     geneOntologySimilarity, geneInteractionSimilarity, 
     exposure, g1exposure, g2exposure,fxNorm_model,fxNorm_t3)

}

####################
# ANALYSIS OF CUTS #
####################

if(do_analyses){
  
  # remove zero and top 20
  d <- neighbor_grab(bplus,chrkey)
  z <- d != 0 & !(d %in% sort(d,decreasing=T)[1:20])
  
  # provide appropriate names to use in the model 
  Proximity <- log(1/neighbor_grab(intrachr_dist,chrkey))[z]
  doubletScore <- neighbor_grab(bplus,chrkey)[z]
  PairOrientation <- pOrient[-chrkey][z]
  TFInventory <- neighbor_grab(jacTF,chrkey)[z]
  BaselineChromatin <- baseWeiner[z]
  DeltaChromatin <- deltaWeiner[z]
  
  # functional similarity 
  geneInteractionSimilarity <- neighbor_grab(giSim,chrkey)[z]
  geneOntologySimilarity <- neighbor_grab(ssWang,chrkey)[z]
  
  # provide exposure variables for both genes in each pair 
  exposure <- apply(b_mod,1,function(x)sum(x!=0))
  g1exposure <- c(NA,exposure[2:length(exposure)])[-chrkey][z]
  g2exposure <- c(NA,exposure[1:(length(exposure)-1)])[-chrkey][z]
  
  # build the model
  noTop20Model <- glm.nb(doubletScore ~ 
                           PairOrientation +
                           TFInventory +
                           Proximity + 
                           BaselineChromatin +
                           DeltaChromatin + 
                           offset(log(g1exposure)) + 
                           offset(log(g2exposure)),
                         na.action = na.exclude,
                         maxit=1000)
  
  noTop20_t3 <- noTop20Model %>% 
    Anova(type='III') %>%
    `names<-`(c('LR_Chisq','Df','pVals'))
  
  
  # format it
  noTop20_t3 <- noTop20_t3 %>% add_column(Factor=rownames(noTop20_t3)) %>%
    mutate(Significance=cut(pVals, breaks=c(-Inf,0.001,0.05,Inf), 
                            labels=c("p < 0.001","p < 0.05","p > 0.05")))
  
  write.csv(noTop20_t3[-(which(names(noTop20_t3) == 'Significance'))] %>% `rownames<-`(c(NULL)),
            file=paste0(wd,'/tables/','revResponse_noTop20Model.csv'))
  
}

################################
# CLEAR ANY REMAINING OBJECTS  #
################################

rm(list=ls())
