########################################################################################
########################################################################################
## CODE FOR: "TRANS-ACTING GENETIC VARIATION AFFECTS THE EXPRESSION OF ADJACENT GENES ##
## CODE WRITTEN BY: BY KRISNA VAN DYKE                                                ##
########################################################################################
########################################################################################

# README
# do not forget to install all the libraries before trying to run this script
# all intermediary objects are saved in their own folder
# all analyses are enclosed in logical satements allow you to run one of three ways
# # generate intermediary objects
# # generate intermediary objects and do analyses
# # generate intermediary objects, do analyses, and generate plots 
# you cannot do analyses without intermediary objects or generate plots without analyses
# do the analyses you want, or turn them all on by setting anlyses to T
# cartoons for the paper were generated separately 
# some figures were superficially modified in post explicitly to
# # change certain font colors
# # change spacing
# # panel figures
# # append p-values generated in this script
# time estimates are provided for some slow parts of the script
# # these are generally listes as (~x min)
# # time-estimates are super crude and for a Microsoft Surface Laptop 2 
# not all analyses are properly fitted with a print statement
# # some need to be added if you want the script to spit them out properly
# you must select a working directory to give the script
# # under wd you must replace SELECT_DIRECTORY with the directory you want to house this script
# you will need the following subfolders in your selected working directory
# # raw
# # intermediary
# # plots
# # tables 
# note that raw should contain all the starting files provided
# # the other three subfolders should only have a placeholder document in them
# old chromatin analyses on Chereji et al. 2018, an Schep et al. 2015, are not included
# analyses on Bergenhold overlap is not included
# # data are held in the raw folder under the name inline-supplementary-material-4.xlsx

#############################
# DIRECTORIES AND LIBRARIES #
#############################

# SET UP WORKING DIRECTORY
wd <-'C:/Users/vandy147/Documents/albert_lab/p01_trans_spreading/code'
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
library(Biostrings) # for string manipulation
library(car) # for the Type III ANOVA
library(GenomicRanges) # used for efficient accessing of basepair data
library(GOSemSim) # used to get pairwise functional similarity
library(Hmisc) # various useful functions
library(magrittr) # for the ever useful pipe operator
library(MASS) # for negative binomial regression
library(org.Sc.sgd.db) # used for SacCer annotations
library(preprocessCore) # used for efficient quantile normalization
library(Rcpp) # allows C++ code to be run in R
library(reshape2) # for melting the correlation matrix
library(rtracklayer) # for the import function and getting things as granges
library(tidyverse) # for some data processing
library(viridis) # for manual coloring using viridis

##################
# RCPP FUNCTIONS #
##################

# NOTE ON N-LET COUNTING
# these were never made into a single clean function or any n-let size
# the idea of non-doublet counting being important came towards the end
# the redundant nature of the counting functions was purely to save time
# had strong reason to believe n-lets couldn't get that large in our case 

# COUNTS HOW MANY DOUBLETS APPEAR IN A VECTOR
# v holds effect scores for a given hotspot
# key holds all relevant values that do not cross a chromosomal boundary
# only considers a doublet if effect is in the same direction
# returns a doublet count
# the key is generally passed from R such as ((1:nrow(b_mod))[-chrkey])-1

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

# GENERATES THE FREQUENCY MATRIX
# if the signs agree you receive a positive value
# if the signs disagree you receive a negative value

cppFunction('NumericMatrix bFrequency(NumericMatrix b) {
  int m = b.nrow();
  int n = b.ncol();
  NumericMatrix s(m,m);
  float d = 2*n;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      for(int k = 0; k < m; k++){
        if(b(j,i)*b(k,i) == 0){
          s(j,k) += 0;
        } else {
          float q;
          if(b(j,i)*b(k,i) > 0){
            q = 1;
          } else {
            q = -1;
          }
          s(j,k) += q;
        }
      }
    }
  }
  
  return(s);
}')

# GENERATES THE BPLUS MATRIX
# counts number of times two genes are both affected by a hotspot
# only counts them if they go in the same direction

cppFunction('NumericMatrix bPlus(NumericMatrix b) {
  int m = b.nrow();
  int n = b.ncol();
  NumericMatrix s(m,m);
  float d = 2*n;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      for(int k = 0; k < m; k++){
        if(b(j,i)*b(k,i) > 0){
          s(j,k) += 1;
        }
      }
    }
  }
  
  return(s);
}')

# COMPLETES PERMUTATION OF HOTSPOTS AND COUNTS DOUBLETS
# shuffling and counting algorithm built in RCpp for speed
# returns a vector of permuted doublet counts
# only considers a doublet if effect is in the same direction
# key holds all relevant values that do not cross a chromosomal boundary
# the key is generally passed from R as # ((1:nrow(b_mod))[-chrkey])-1

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

# TRIPLET VERSIONS OF THE DOUBLET COUNTING FUNCTIONS

cppFunction('int triCount(NumericVector v, IntegerVector key) {
  int m = key.size();
  int q = 0;
    
  for (int j = 0; j < m; j++){
  
    if (((v[key[j]] > 0) & (v[key[j]-1] > 0) & (v[key[j]-2] > 0)) | 
    ((v[key[j]] < 0) & (v[key[j]-1] < 0) & (v[key[j]-2] < 0))){
      q++;
    }
  }
  
  return q;
}')

cppFunction('IntegerVector triPerm(int t, NumericVector hs, IntegerVector key) {
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
        if (((hs[key[k]] > 0) & (hs[key[k]-1] > 0) & (hs[key[k]-2] > 0)) | 
        ((hs[key[k]] < 0) & (hs[key[k]-1] < 0) & (hs[key[k]-2] < 0))){
            out[i]++;
        }
    }
    
  }

  return out;
}')

# QUADRUPLET VERSIONS OF THE DOUBLET COUNTING FUNCTIONS

cppFunction('int quadCount(NumericVector v, IntegerVector key) {
  int m = key.size();
  int q = 0;
    
  for (int j = 0; j < m; j++){
  
    if (((v[key[j]] > 0) & (v[key[j]-1] > 0) & (v[key[j]-2] > 0) & (v[key[j]-3] > 0)) |
    ((v[key[j]] < 0) & (v[key[j]-1] < 0) & (v[key[j]-2] < 0) & (v[key[j]-3] < 0))){
      q++;
    }
  }
  
  return q;
}')

cppFunction('IntegerVector quadPerm(int t, NumericVector hs, IntegerVector key) {
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
        if (((hs[key[k]] > 0) & (hs[key[k]-1] > 0) & (hs[key[k]-2] > 0) & (hs[key[k]-3] > 0)) |
        ((hs[key[k]] < 0) & (hs[key[k]-1] < 0) & (hs[key[k]-2] < 0) & (hs[key[k]-3] < 0))){
            out[i]++;
        }
    }
    
  }

  return out;
}')

# QUINTUPLET VERSIONS OF THE DOUBLET COUNTING FUNCTIONS

cppFunction('int quinCount(NumericVector v, IntegerVector key) {
  int m = key.size();
  int q = 0;
    
  for (int j = 0; j < m; j++){
  
    if (((v[key[j]] > 0) & (v[key[j]-1] > 0) & (v[key[j]-2] > 0) & (v[key[j]-3] > 0) & (v[key[j]-4] > 0)) |
    ((v[key[j]] < 0) & (v[key[j]-1] < 0) & (v[key[j]-2] < 0) & (v[key[j]-3] < 0) & (v[key[j]-4] < 0))){
      q++;
    }
  }
  
  return q;
}')

cppFunction('IntegerVector quinPerm(int t, NumericVector hs, IntegerVector key) {
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
        if (((hs[key[k]] > 0) & (hs[key[k]-1] > 0) & (hs[key[k]-2] > 0) & (hs[key[k]-3] > 0) & (hs[key[k]-4] > 0)) |
        ((hs[key[k]] < 0) & (hs[key[k]-1] < 0) & (hs[key[k]-2] < 0) & (hs[key[k]-3] < 0) & (hs[key[k]-4] < 0))){
            out[i]++;
        }
    }
    
  }

  return out;
}')

cppFunction('int sextCount(NumericVector v, IntegerVector key) {
  int m = key.size();
  int q = 0;
    
  for (int j = 0; j < m; j++){
  
    if (((v[key[j]] > 0) & (v[key[j]-1] > 0) & (v[key[j]-2] > 0) & 
    (v[key[j]-3] > 0) & (v[key[j]-4] > 0) & (v[key[j]-5] > 0)) |
    ((v[key[j]] < 0) & (v[key[j]-1] < 0) & (v[key[j]-2] < 0) & 
    (v[key[j]-3] < 0) & (v[key[j]-4] < 0) & (v[key[j]-5] < 0))){
      q++;
    }
  }
  
  return q;
}')

cppFunction('IntegerVector sextPerm(int t, NumericVector hs, IntegerVector key) {
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
        if (((hs[key[k]] > 0) & (hs[key[k]-1] > 0) & (hs[key[k]-2] > 0) & 
        (hs[key[k]-3] > 0) & (hs[key[k]-4] > 0) & (hs[key[k]-5] > 0)) |
        ((hs[key[k]] < 0) & (hs[key[k]-1] < 0) & (hs[key[k]-2] < 0) & 
        (hs[key[k]-3] < 0) & (hs[key[k]-4] < 0) & (hs[key[k]-5] < 0))){
            out[i]++;
        }
    }
    
  }

  return out;
}')

###############
# R FUNCTIONS #
###############

# assumes that x and y are the same length
jaccard <- function(x,y){
  z <- x+y
  return(sum(z == 2)/sum(z > 0))
}

# make sure each dataset has the same columns (NA where appropriate)
conform_universe <- function(df,key){
  df <- df[rownames(df) %in% key,]
  missing <- key[!(key %in% rownames(df))]# exclaimation debacle???
  nulls <- data.frame(matrix(NA,nrow=length(missing),ncol=ncol(df)))
  rownames(nulls) <- missing
  colnames(nulls) <- colnames(df)
  df <- rbind(df,nulls)
  df <- df[match(key, row.names(df)),]
  return(df)
}

# annotate this
norm_wrapper <- function(df){
  row_hold <- rownames(df)
  col_hold <- colnames(df)
  df <- normalize.quantiles(as.matrix(df))
  colnames(df) <- col_hold
  rownames(df) <- row_hold
  return(df)
}

# gets Spearman correlation matrices 
spear_cor <- function(dir,x,r=F){
  if(file.exists(paste0(dir,x,'_cor.RData'))){
    load(paste0(dir,x,'_cor.RData'))
  } else {
    load(paste0(dir,x,'.RData'))
    assign(paste0(x,'_cor'),rcorr(get(x), type='spearman'))
    save(list = paste0(x,'_cor'),file=paste0(dir,x,'_cor.RData'))
  }
  if(r){
    return(get(paste0(x,'_cor'))$r)
  }
  rm(list = c(paste0(x,'_cor'),paste0(x)))
}

# changes the format of strandedness
strand_formatter <- function(df){
  for(i in 1:nrow(df)){
    if(df$strand[i] > 0){
      df$strand[i] <- '+'
    } else {
      df$strand[i] <- '-'
    }
  }
  return(df)
}

# assures that negative strand genes are properly flipped
strand_flip <- function(df){
  for(i in 1:nrow(df)){
    if((df$strand[i] == -1 | df$strand[i] == '-') & (df$start[i] < df$end[i])){
      temp <- df$start[i] 
      df$start[i] <- df$end[i]
      df$end[i] <- temp
    }
  }
  return(df)
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

# create a key that says which genes are interchromosomal neighbors (not true neighbors)
chrkey_gen <- function(x){
  a <- table(x)
  b <- c(rep(0,length(a)))
  for(i in 1:length(a)){
    b[i] <- sum(a[1:i])
  }
  b <- b + 1
  b <- c(1,b[-length(b)])
  return(b)
}

# gets the adjacent-to-diagonal elements of a matrix
diagonal_grabber <- function(x){
  d <- c(rep(NA,ncol(x)))
  for (i in 2:length(d)){
    d[i] <- x[i,i-1]
  }
  return(d)
}

# gets elements of a matrix beyond the adjacent-to-diagonal
# these can also be thought of being two from the center diagonal (tc)
tc_grab <- function(x){
  y <- x[-(nrow(x)),-1][upper.tri(x[-(nrow(x)),-1])]
  return(y)
}

# retrieves neighbors and removes false neighbors
neighbor_grab <- function(x,key){
  diagonal_grabber(x)[-key]
}

# retrieves nonneighbors and adds false neighbors at the start
nonneighbor_grab <-function(x,key){
  c(diagonal_grabber(x)[key[-1]],tc_grab(x))
}

# makes elements of a desired sign into NA
matsign <- function(x,s){
  if(s > 0){
    x[x < 0] <- NA
  } else if (s == 0){
    x
  } else {
    x[x > 0] <- NA
  }
  return(x)
}

# R-BASED DOUBLET COUNTS (SLOW)
# can act as a sanity check for the Rcpp counting logic
# gathers targets uniquely contributing to at least 1 doublet
# returns a named vector with the doublet count and number unique contributors 
if(check_sanity){
  dUniq <- function(v,c){
    u <- c()
    q <- 0
    for(t in 1:length(c)){
      i <- c[t]
      if((v[i] > 0 & v[i-1] > 0) | (v[i] < 0 & v[i-1] < 0)){
        q <- q + 1
        u <- append(u,c(i,i-1))
      }
    }
    r <- c(q,length(unique(u)))
    names(r) <- c('Doublets','Unique')
    return(r)
  } 
  
}


#########################
# PART I. PREPROCESSING # 
#########################
########################################
# PARSE POSITIONS OF VARIANTS IN CROSS #
########################################

# Holds the raw counts and BY/RM allele information
load(paste0(wd,'/raw/Counts.RData'))

# reformats certain information about variant position 
load(paste0(wd,'/intermediary/relation_table.RData'))
if(!exists('relation_table')){
  relation_table <- as.data.frame(matrix(nrow=ncol(counts$gdata),ncol=4))
  colnames(relation_table) <- c('original','chromosome','position','variant')
  relation_table$original <- colnames(counts$gdata)
  
  for (i in 1:nrow(relation_table)){
    relation_table[i,2:4] <- strsplit(relation_table[i,1], "_|:")[[1]]
  }
  
  chromosome_dictionary <- data.frame(matrix(ncol=16,nrow=1))
  chromosome_dictionary[1,] <- 1:16
  colnames(chromosome_dictionary) <- c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII',
                                       'chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
  
  for (i in 1:nrow(relation_table)){
    relation_table[i,2] <- chromosome_dictionary[relation_table[i,2]]
  }
  
  relation_table[,2] <- as.integer(relation_table[,2])
  relation_table[,3] <- as.integer(relation_table[,3])
  
  
  rm(chromosome_dictionary)
  save(relation_table, file=paste0(wd,'/intermediary/relation_table.RData'))
}

######################################
# BUILD GENE MODELS METADATA OBJECTS #
######################################

# this metadata object holds the gene models and is frequently called
load(paste0(wd,'/intermediary/metaNT.RData'))
load(paste0(wd,'/intermediary/chrkey.RData'))
load(paste0(wd,'/intermediary/grBodies.RData'))
if(!(exists('metaNT') & exists('chrkey') & exists('grBodies'))){
  
  # load in the metadata for all SGD verified ORFs
  metaNT <- read.csv(paste0(wd,'/raw/vORFs.csv'), check.names = F)
  
  # order correctly 
  metaNT <- metaNT [order(metaNT$chromosome,metaNT$start),]
  
  # Change factor levels
  metaNT <- mutate(metaNT,chromosome = 
                     recode_factor(chromosome,
                                   chrI='1',chrII='2',chrIII='3',chrIV='4',
                                   chrV='5',chrVI='6',chrVII='7',chrVIII='8',
                                   chrIX = '9',chrX ='10',chrXI='11',chrXII='12',
                                   chrXIII='13',chrXIV='14',chrXV='15',chrXVI='16')) %>% 
    mutate(chromosome=as.numeric(as.character(chromosome)))
  
  # remove the tiled ASP3 genes and other genes found in the rRNA cluster
  metaNT <- metaNT[-((which(metaNT$systematic_name =='YLR154C')+1):(which(metaNT$systematic_name =='YLR163C')-1)),]
  
  # load the tpm data
  load(paste0(wd,'/raw/tpm.RData'))
  
  # restrict the list of genes to verified ORFs expressed in the Albert & Bloom, et al. 2018 dataset
  metaNT <- metaNT[metaNT$systematic_name %in% colnames(tpm),]
  
  # remove the tpm data... for now
  rm(tpm)
  
  # generate the marker cutoffs telling us what is not telomeric  
  marker_cutoffs <- sapply(1:max(relation_table$chromosome),function(x)c(
    min(relation_table[relation_table$chromosome == x,]$position),
    max(relation_table[relation_table$chromosome == x,]$position)))
  marker_cutoffs <- as.data.frame(t(rbind(1:max(relation_table$chromosome),marker_cutoffs)))
  colnames(marker_cutoffs) <- c('chromosome','left','right')
  
  # apply the conditions that the gene cannot overlap the telomere in any way
  metaNT <- do.call(rbind, apply(marker_cutoffs,1,function(x)metaNT[(
    metaNT$chromosome==x[1] 
    & metaNT$start > x[2]
    & metaNT$end > x[2]
    & metaNT$start < x[3]
    & metaNT$end < x[3]
  ),]))
  
  # add the skipped rRNA region to the chromosome breaks (should be 3267)
  chrkey <- c(chrkey_gen(metaNT$chromosome),which(metaNT$systematic_name =='YLR163C')) %>% sort()
  save(chrkey,file=paste0(wd,'/intermediary/chrkey.RData'))
  
  # prepare genomic ranges objects needed for nucleosome occupancy data
  grBodies <- metaNT[c('start','end','strand','chromosome')] %>% 
    strand_formatter() %>% 
    makeGRangesFromDataFrame()
  save(grBodies,file=paste0(wd,'/intermediary/grBodies.RData'))
  
  # flip stars of negative genes
  metaNT <- strand_flip(metaNT)
  
  # prevent any matching tasks from trying to match factors
  metaNT$systematic_name <- as.character(metaNT$systematic_name)
  metaNT$standard_name <- as.character(metaNT$standard_name)
  
  save(metaNT,file=paste0(wd,'/intermediary/metaNT.RData'))
  rm(marker_cutoffs)
  
}

rm(relation_table)

# FIGURE OUT PAIR ORIENTATIONS
# we already know these are sorted by position of start codon
# from this know pair orientation based on the "lead" and "lag" being compared 
# these are therefore, only defined by the genes EXPRESSED in the BY/RM cross
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

rm(lead,lag)

# GET THE BASEPAIR DISTANCE MATRIX

# what is on the same chromosome is T
temp <- distancer(metaNT$chromosome) == 0
# make sure everything not on the same chromosome becomes NA
temp[!temp] <- NA
# get the values
intrachr_dist <- distancer(metaNT$start)*temp 
rm(temp)

##############################
# PREPARE THE CHROMATIN DATA #
##############################

# Weiner et al., 2015 ChIP-seq

# Table S2. Nucleosome Atlas (nucleosome positions)
nucPos <- read.delim(paste0(wd,'/raw/','mmc3.csv'),sep=',')

# Table S3. Normalized Modification Levels (chromatin marks over conditions)
norMark <- read.delim(paste0(wd,'/raw/','Table_S3.csv'),sep=',')

# keep only the normalized level of marks
norMark <- norMark[-1,2:157]

# make sure all elements of the matrix are properly numeric
norMark <- apply(norMark,2,as.numeric)

# correlation of marks on +1 nucleosomes of neighbors at t0
if(!file.exists(paste0(wd,'/intermediary/baseWeiner.RData'))){
  firstNucs <- nucPos[nucPos$gene_pos == 1 & !is.na(nucPos$gene_pos),]
  rownames(firstNucs) <- firstNucs$acc
  firstNucs <- conform_universe(firstNucs,metaNT$systematic_name)
  baseWeiner <- rep(NA,nrow(firstNucs))
  for(i in 2:nrow(firstNucs)){
    if(!is.na(firstNucs$nuc_id[i-1]) & !is.na(firstNucs$nuc_id[i])){
      baseWeiner[i] <- cor.test(norMark[firstNucs$nuc_id[i-1],-seq(1,151,6)],
                                norMark[firstNucs$nuc_id[i],-seq(1,151,6)],
                                method='spearman')$estimate
    }
  }
  baseWeiner <- baseWeiner[-chrkey]
  save(baseWeiner,file=paste0(wd,'/intermediary/baseWeiner.RData'))
  rm(baseWeiner, firstNucs)
}


# correlation of marks on +1 nucleosomes of neighbors at t15
if(!file.exists(paste0(wd,'/intermediary/deltaWeiner.RData'))){
  firstNucs <- nucPos[nucPos$gene_pos == 1 & !is.na(nucPos$gene_pos),]
  rownames(firstNucs) <- firstNucs$acc
  firstNucs <- conform_universe(firstNucs,metaNT$systematic_name)
  deltaWeiner <- rep(NA,nrow(firstNucs))
  for(i in 2:nrow(firstNucs)){
    if(!is.na(firstNucs$nuc_id[i-1]) & !is.na(firstNucs$nuc_id[i])){
      deltaWeiner[i] <- cor.test(abs(norMark[firstNucs$nuc_id[i-1],seq(4,154,6)]-norMark[firstNucs$nuc_id[i-1],seq(1,151,6)]),
                                 abs(norMark[firstNucs$nuc_id[i],seq(4,154,6)]-norMark[firstNucs$nuc_id[i],seq(1,151,6)]),
                                 method='spearman')$estimate
    }
  }
  deltaWeiner <- deltaWeiner[-chrkey]
  save(deltaWeiner,file=paste0(wd,'/intermediary/deltaWeiner.RData'))
  rm(deltaWeiner, firstNucs)
}

rm(nucPos,norMark)

########################################################
# CREATE A B-MATRIX WITH GENES AFFECTED BY NO HOTSPOTS #
########################################################

load(paste0(wd,'/intermediary/b_mod.RData'))
if(!exists('b_mod')){
  b_matrix <- read.csv(paste0(wd,'/raw/b_matrix.csv'), header=T, row.names=1, check.names = F)
  b_mod <- conform_universe(b_matrix,metaNT$systematic_name)
  rm(b_matrix)
  b_mod[is.na(b_mod)] <- 0
  save(b_mod, file=paste0(wd,'/intermediary/b_mod.RData'))
}

############################################
# PROCESS 2018 BY/RM CROSS EXPRESSION DATA #
############################################

# create the batch and OD corrected data
load(paste0(wd,'/intermediary/mAlbert.RData'))
if(!exists('mAlbert')){
  
  # load the transcript per million counts you are going to start with
  load(paste0(wd,'/raw/tpm.RData'))
  
  # load in batch and OD information with batch as factor
  covariates <- read.csv(paste0(wd,'/raw/covariates.csv'), check.names = F)
  covariates$batch <- as.factor(covariates$batch)
  
  # correct for the effects of batch and optical density 
  mAlbert <- apply(tpm,2,function(x)lm(x ~ covariates$batch + covariates$OD)$residuals)
  
  # susbet tpm by what is in metaNT, because we know these are non-telomeric
  mAlbert <- mAlbert[,(colnames(mAlbert) %in% metaNT$systematic_name)]
  
  # then order genes by postition on chromosome
  mAlbert <- mAlbert[,match(colnames(mAlbert),metaNT$systematic_name)]
  
  # save your new object
  save(mAlbert,file=paste0(wd,'/intermediary/mAlbert.RData'))
  
  # remove extraneous objects
  rm(covariates)
}

# create the eQTL corrected expression data (~4 min)
load(paste0(wd,'/intermediary/Albert18.RData'))
if(!exists('Albert18')){
  
  # you will need the eQTL data in addition to the counts data to get this properly done
  eQTL_data <- read.csv(paste0(wd,'/raw/eQTL_data.csv'),check.names=F)
  
  Albert18 <- mAlbert
  
  for(i in 1:ncol(mAlbert)){
    markers <- which(colnames(counts$gdata) %in% as.character(eQTL_data[eQTL_data$gene == (colnames(Albert18)[i]),]$pmarker))
    if(length(markers) != 0){
      Albert18[,i] <- lm(Albert18[,i] ~ counts$gdata[,markers])$residuals
    }
  }
  
  # center each distribution on zero
  Albert18 <- apply(Albert18,1,function(x)x-mean(x,na.rm=T))
  
  # quantile normalize the data
  Albert18 <- Albert18 %>% norm_wrapper %>% t
  
  save(Albert18,file=paste0(wd,'/intermediary/Albert18.RData'))
  rm(eQTL_data,markers)
  
}

rm(counts)

#################################################
# GET QUALITY CONTROL ON THE CORRECTION PROCESS #
#################################################

if(make_plots){
  templot <- rbind(
    mAlbert[,which(metaNT$standard_name=='HO')] %>% as.data.frame() %>% 
      `names<-`(c('Residuals')) %>% add_column(Correction='Only Batch + OD'),
    Albert18[,which(metaNT$standard_name=='HO')]  %>% as.data.frame() %>% 
      `names<-`(c('Residuals')) %>% add_column(Correction='Batch + OD + Genetics')) %>% 
    ggplot(aes(x=Residuals,fill=Correction)) + 
    geom_density(size=0.5, alpha=0.4) + theme_bw() + scale_fill_viridis_d() + 
    theme(legend.position="bottom")  + xlab('Residual Gene Expression') + 
    ylab('Density') + ggtitle('HO')
  
  ggsave(paste0(wd,'/plots/','SupplFigure04_ResidualExpressionHO.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 3.5,height = 3,units = c("in"))
  rm(templot)
  
  templot <- rbind(
    mAlbert[,which(metaNT$standard_name=='STE2')] %>% as.data.frame() %>% 
      `names<-`(c('Residuals')) %>% add_column(Correction='Only Batch + OD'),
    Albert18[,which(metaNT$standard_name=='STE2')]  %>% as.data.frame() %>% 
      `names<-`(c('Residuals')) %>% add_column(Correction='Batch + OD + Genetics')) %>% 
    ggplot(aes(x=Residuals,fill=Correction)) + 
    geom_density(size=0.5, alpha=0.4) + theme_bw() + scale_fill_viridis_d() +
    theme(legend.position="bottom") + xlab('Residual Gene Expression') + 
    ylab('Density') + ggtitle('STE2')
  
  ggsave(paste0(wd,'/plots/','SupplFigure04_ResidualExpressionSTE2.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 3.5,height = 3,units = c("in"))
  rm(templot) 
}

rm(mAlbert,Albert18)

##########################################
# PREPARE EXTERNAL COEXPRESSION DATASETS #
##########################################

# CLEAN UP AND PREPARE ALL EXTERNAL DATASETS FOR USE
# SPECIFICALLY THINK ABOUT THE FLEMING 99 BEING MORE LIKE FLEMING 02

if(!file.exists(paste0(wd,'/intermediary/Brem05.RData'))){
  
  Brem05_swap <- read.table(paste0(wd,"/raw/Brem05_dyeSwap.flt.knn.avg.pcl"),sep ="\t", header=T,row.names=1,check.names=F)
  rownames(Brem05_swap) <- Brem05_swap[,1]
  Brem05_swap <- Brem05_swap[-1,-1:-2]
  
  Brem05_orig <- read.table(paste0(wd,"/raw/Brem05_orig.flt.knn.avg.pcl"),sep="\t",header=T,row.names=1,check.names=F)
  rownames(Brem05_orig) <- Brem05_orig[,1]
  Brem05_orig <- Brem05_orig[-1,-1:-2]
  
  
  # get the mean of the original and the swap to ensure no "dye-effect'
  Brem05 <- Brem05_orig
  for(i in 1:nrow(Brem05)){
    for(j in 1:ncol(Brem05)){
      Brem05[i,j] <- mean(c(Brem05_orig[i,j],Brem05_swap[i,j]))
    }
  }
  
  # remove the dye swaps after averaging
  rm(Brem05_swap, Brem05_orig)
  
  # conform universes and normalize quantiles without forgetting names
  Brem05 <- conform_universe(Brem05,metaNT$systematic_name)
  
  # save with genes as columns
  Brem05 <- t(Brem05)
  save(Brem05,file=paste0(wd,'/intermediary/Brem05.RData'))
  rm(Brem05)
  
}

if(!file.exists(paste0(wd,'/intermediary/Fleming02.RData'))){
  
  # read in information and drop unnecessary information
  Fleming02 <- read.table(paste0(wd,"/raw/Fleming99.log.flt.knn.avg.pcl"),sep="\t",header=T,row.names=1,check.names=F)
  rownames(Fleming02) <- Fleming02[,1]
  Fleming02 <- Fleming02[-1,-1:-2]
  
  # conform universes and normalize quantiles without forgetting names
  Fleming02 <- conform_universe(Fleming02,metaNT$systematic_name)
  
  # save with genes as columns
  Fleming02 <- t(Fleming02)
  save(Fleming02,file=paste0(wd,'/intermediary/Fleming02.RData'))
  rm(Fleming02)
  
}

if(!file.exists(paste0(wd,'/intermediary/Hughes00.RData'))){
  
  # read in data
  Hughes00 <- read.table(paste0(wd,"/raw/Hughes00.flt.knn.avg.pcl"),sep="\t",quote="\"",header=T,row.names=1,check.names=F)
  Hughes00 <- Hughes00[-1,-1:-2]
  
  # conform universes and normalize quantiles without forgetting names
  Hughes00 <- conform_universe(Hughes00,metaNT$systematic_name)
  
  # save with genes as columns
  Hughes00 <- t(Hughes00)
  save(Hughes00,file=paste0(wd,'/intermediary/Hughes00.RData'))
  rm(Hughes00)
  
}

if(!file.exists(paste0(wd,'/intermediary/Klevecz04.RData'))){
  
  # load raw data
  Klevecz04 <- read.table(paste0(wd,'/raw/Klevecz04_SPELL.pcl'),sep="\t",header=T,row.names=1,check.names=F)
  
  # remove unnecessary information
  Klevecz04 <- Klevecz04[-1,-1:-2]
  
  # mean normalize data
  Klevecz04 <- apply(Klevecz04,2,function(x)x-mean(x))
  
  # conform universes and normalize quantiles withotu forgetting names
  Klevecz04 <- conform_universe(Klevecz04,metaNT$systematic_name)
  Klevecz04 <- norm_wrapper(Klevecz04)
  
  # save with genes as columns
  Klevecz04 <- t(Klevecz04)
  save(Klevecz04,file=paste0(wd,'/intermediary/Klevecz04.RData'))  
  rm(Klevecz04)
  
}

if(!file.exists(paste0(wd,'/intermediary/Knijnenburg09.RData'))){
  
  Knijnenburg09 <- read.table(paste0(wd,"/raw/Knijnenburg09_SPELL.pcl"),sep="\t",header=T,row.names=1,check.names=F,comment.char = "")
  rownames(Knijnenburg09) <- Knijnenburg09[,1]
  Knijnenburg09 <- Knijnenburg09[-1,-1:-2]
  
  # conform universes and normalize quantiles without forgetting names
  Knijnenburg09 <- conform_universe(Knijnenburg09,metaNT$systematic_name)
  Knijnenburg09 <- norm_wrapper(Knijnenburg09)
  
  # save with genes as columns
  Knijnenburg09 <- t(Knijnenburg09)
  save(Knijnenburg09,file=paste0(wd,'/intermediary/Knijnenburg09.RData'))
  rm(Knijnenburg09)
  
}

if(!file.exists(paste0(wd,'/intermediary/Lenstra11.RData'))){
  
  # load the table and get your systematic as characters
  Lenstra11 <- read.table(paste0(wd,"/raw/Lenstra11_raw.txt"), header=T, quote="\"", sep='\t',check.names=F)
  Lenstra11$systematicName <- as.character(Lenstra11$systematicName)
  
  # there are 11 duplicated items to be removed (22 rows total)
  Lenstra11 <- Lenstra11[!(Lenstra11$systematicName %in% Lenstra11$systematicName[duplicated(Lenstra11$systematicName)]),]
  
  # remove P-values so we are only left with fold-change and then remove that identifier row
  Lenstra11 <- Lenstra11[,Lenstra11[1,] != 'p_value']
  Lenstra11 <- Lenstra11[-1,-2]
  
  # find the names of the duplicated items
  duplicated_names <- Lenstra11$systematicName[duplicated(Lenstra11$systematicName)]
  
  # remove duplicated items by name
  Lenstra11 <- Lenstra11[!(Lenstra11$systematicName %in% duplicated_names),]
  # clear your duplicated names, no longer necessary
  rm(duplicated_names)
  
  # get names attached
  rownames(Lenstra11) <- Lenstra11[,1]
  Lenstra11 <- Lenstra11[,-1]
  
  # conform universes
  Lenstra11 <- conform_universe(Lenstra11,metaNT$systematic_name)
  
  # turn to floats
  for(i in 1:ncol(Lenstra11)){
    Lenstra11[,i] <- as.numeric(levels(Lenstra11[,i])[Lenstra11[,i]])
  }
  
  # normalize quantiles without forgetting names
  Lenstra11 <- norm_wrapper(Lenstra11)
  
  # save with genes as columns
  Lenstra11 <- t(Lenstra11)
  save(Lenstra11,file=paste0(wd,'/intermediary/Lenstra11.RData'))  
  rm(Lenstra11)
  
}

if(!file.exists(paste0(wd,'/intermediary/Myers19.RData'))){
  
  # read in the data
  Myers19 <- read.table(paste0(wd,"/raw/Myers19.tsv"),sep="\t",header=T,row.names=1,check.names=F,quote="\"")
  
  # remove unnecessary information
  Myers19 <- Myers19[,-1]
  
  # conform the universe without forgetting the names
  # # note this dataset required mean normalization due to its high center
  Myers19 <- conform_universe(Myers19,metaNT$systematic_name)
  Myers19 <- apply(Myers19,2,function(x)x-mean(x,na.rm=T))
  Myers19 <- norm_wrapper(Myers19)
  
  # save with genes as columns
  Myers19 <- t(Myers19)
  save(Myers19,file=paste0(wd,'/intermediary/Myers19.RData'))
  rm(Myers19)
  
}

if(!file.exists(paste0(wd,'/intermediary/Sameith15.RData'))){
  
  # get names in correct places
  Sameith15 <- read.table(paste0(wd,"/raw/Sameith15_SPELL.pcl"),sep="\t",header=T,row.names=1,check.names=F)
  rownames(Sameith15) <- Sameith15[,1]
  Sameith15 <- Sameith15[-1,-1:-2]
  
  # conform universes and normalize quantiles without forgetting names
  Sameith15 <- conform_universe(Sameith15,metaNT$systematic_name)
  Sameith15 <- norm_wrapper(Sameith15)
  
  # save with genes as columns
  Sameith15 <- t(Sameith15)
  save(Sameith15,file=paste0(wd,'/intermediary/Sameith15.RData'))
  rm(Sameith15)
  
}

if(!file.exists(paste0(wd,'/intermediary/Schurch16.RData'))){
  
  # read in data
  Schurch16 <- read.table(paste0(wd,'/raw/Schurch16_raw.tsv'),header=T,check.names=F,row.name=1)
  rownames(Schurch16) <- toupper(rownames(Schurch16))
  
  # read in the lengths of genes
  gene_lengths <- read.csv(paste0(wd,'/raw/gene_lengths.csv'),header=T,check.names=F)
  
  # conform universes to 2018 segregant data
  Schurch16 <- conform_universe(Schurch16,metaNT$systematic_name)
  
  # remove it if it does not appear as at least 1 read in half or more 
  Schurch16[!apply(Schurch16,1,function(x)sum(x > 0) > length(x)/2),] <- NA
  
  # get the gene lengths for each isogenic entry
  Schurch16_lengths <- gene_lengths[gene_lengths$Gene_Name %in% rownames(Schurch16),]
  
  # match the order of the gene length directory to the order of the isogenic data
  Schurch16_lengths <- gene_lengths[order(match(Schurch16_lengths$Gene_Name,rownames(Schurch16))),]
  
  # normalize reads for gene size
  for(i in 1:nrow(Schurch16)){
    Schurch16[i,] <- Schurch16[i,]/Schurch16_lengths$Gene_Length[i]
  }
  
  # normalize so that 0s are not improperly converted into log-space
  Schurch16 <- Schurch16 + 0.5
  
  # enter log space for all reads
  Schurch16 <- log(Schurch16)
  
  # center all items based on mean
  Schurch16 <- apply(Schurch16,2,function(x)x-mean(x,na.rm=T))
  
  # normalize quantiles without forgeting names
  Schurch16 <- norm_wrapper(Schurch16)
  
  # genes as columns - NOTE: data Schurch has a long right tail
  Schurch16 <- t(Schurch16)
  save(Schurch16,file=paste0(wd,'/intermediary/Schurch16.RData'))  
  rm(Schurch16_lengths, Schurch16)
}

if(!file.exists(paste0(wd,'/intermediary/Simola10.RData'))){
  
  # read in data and remove the Paradoxus samples
  Simola10 <- read.table(paste0(wd,"/raw/Simola10_NATURAL_IMPUTE.tsv"),sep="\t",header=T,row.names=1,check.names=F) %>%
    dplyr::select(-contains("YPS3395"))
  
  # conform universes and normalizes quantiles without forgetting names
  Simola10 <- conform_universe(Simola10,metaNT$systematic_name)
  Simola10 <- norm_wrapper(Simola10)
  
  # save with genes as columns
  Simola10 <- t(Simola10)
  save(Simola10,file=paste0(wd,'/intermediary/Simola10.RData'))
  rm(Simola10)
  
}

###########################################
# TRANSCRIPTION FACTOR PROFILE SIMILARITY #
###########################################

# Jaccard index for transcription factor similarity (~5 min)
if(!file.exists(paste0(wd,'/intermediary/jacTF.RData'))){
  binaryTF <- read.delim(paste0(wd,'/raw/RegulationMatrix_Documented_2020123_1648_427964113.csv'),
                         sep=';',row.names=1,check.names=F) %>% 
    t() %>% as.data.frame() %>% (function(x){
      for(i in 1:nrow(x)){
        if(rownames(x)[i] %in% metaNT$standard_name){
          temp <- which(metaNT$standard_name == rownames(x)[i])
          rownames(x)[i] <- toString(metaNT$systematic_name[temp])
        }
      }
      return(x)}) %>% 
    (function(x){conform_universe(x,metaNT$systematic_name)})
  jacTF <- apply(binaryTF,1,function(x)apply(binaryTF,1,function(y)jaccard(x,y)))
  save(jacTF,file=paste0(wd,'/intermediary/jacTF.RData'))
  rm(jacTF,binaryTF)
}

##############################
# FUNCTIONAL SIMILARITY (GO) #
##############################

# Semantic similarity matrix (~36 min)
if(!file.exists(paste0(wd,'/intermediary/ssWang.RData'))){
  
  # get the SC ontology data
  scGO <- godata('org.Sc.sgd.db', ont="BP")
  
  # get the sc map object
  scMap <- AnnotationDbi::select(
    org.Sc.sgd.db,
    keys = as.character(metaNT$systematic_name),
    columns = c("ENTREZID", "GENENAME", "SGD")
  )
  
  # get the gene similarity - this is the slow step 
  ssWang <- mgeneSim(genes=scMap$ENTREZID, semData=scGO, measure="Wang",verbose=F)
  
  # Match the ENTREZ IDs back to their names
  colnames(ssWang) <- scMap$ORF[match(colnames(ssWang),scMap$ENTREZID)]
  rownames(ssWang) <- scMap$ORF[match(rownames(ssWang),scMap$ENTREZID)]
  
  # Figure out what the mgeneSim algorithm dropped because it automatically kills NA
  dropped <- scMap$ORF[!(scMap$ORF %in% colnames(ssWang))]
  
  # Give the names of the dropped to what you will cbind
  temp1 <- as.data.frame(matrix(NA,ncol=length(dropped),nrow=nrow(ssWang)))
  colnames(temp1) <- dropped
  
  # Give the names of the dropped to what you will rbind
  temp2 <- as.data.frame(matrix(NA,ncol=nrow(scMap),nrow=length(dropped)))
  rownames(temp2) <- dropped
  
  # Bind first set of NA columns
  ssWang <- cbind(ssWang,temp1)
  
  # Provide names to make next binding possible
  colnames(temp2) <- colnames(ssWang)
  
  # Bind second set of NA columns
  ssWang <- rbind(ssWang,temp2)
  
  # sort the matrix row and column-wise back into its proper order
  ssWang <- ssWang[match(metaNT$systematic_name,rownames(ssWang)),]
  ssWang <- ssWang[,match(metaNT$systematic_name,colnames(ssWang))]
  
  # Return us to a matrix format
  ssWang <- as.matrix(ssWang)
  
  # Get rid of unnecessary objects
  rm(temp1,temp2, scMap, dropped)
  
  # save data
  save(ssWang,file=paste0(wd,'/intermediary/ssWang.RData'))
  rm(ssWang)
}

#################
# GI SIMILARITY #
#################

# load and format genetic interaction similarity matrix
if(!file.exists(paste0(wd,'/intermediary/giSim.RData'))){
  
  # read in the genetic interaction similarity matrix as a matrix
  unzip(zipfile = paste0(wd,'/raw/','cc_ALL.zip'), exdir = paste0(wd,'/raw'))
  giSim <- read.delim(paste0(wd,'/raw/cc_ALL.txt'),sep='\t', header=F,check.names=F,colClasses='character')
  giSim <- giSim %>% as.matrix
  
  # read in the genetic interaction similarity matrix as a matrix
  giSim <- read.delim(paste0(wd,'/raw/cc_ALL.txt'),sep='\t', header=F,check.names=F,colClasses='character')
  giSim <- giSim %>% as.matrix
  
  # figure out what is unique and get names for unique entries
  uniq <- !(giSim[,2] %in% giSim[(duplicated(giSim[,2])),2])
  nUniq <- giSim[uniq,2]
  
  # remove what is not unique
  giSim <- giSim[uniq,uniq]
  
  # change values from text to numeric
  giSim <- apply(giSim,2,as.numeric)
  
  # reassign column and row names
  colnames(giSim) <- nUniq
  rownames(giSim) <- nUniq
  
  # get rid of any row or column not in the metadata systematic names
  giSim <- giSim[(rownames(giSim) %in% metaNT$systematic_name),
                 (colnames(giSim) %in% metaNT$systematic_name)]
  
  # figure out how many items are missing between the matrix and metadata
  diff <- nrow(metaNT) - nrow(giSim)
  
  # add NA columns for items missing in the matrix that are in the metadata
  dropped <- metaNT$systematic_name[!(metaNT$systematic_name %in% rownames(giSim))]
  giSim <- cbind(giSim, matrix(NA,ncol=length(dropped),nrow=nrow(giSim)) %>% `colnames<-`((dropped)))
  giSim <- rbind(giSim, matrix(NA,ncol=ncol(giSim),nrow=length(dropped)) %>% `rownames<-`((dropped)))
  
  # order the matrix based on the order in the metadata
  giSim <- giSim[match(metaNT$systematic_name,rownames(giSim)),]
  giSim <- giSim[,match(metaNT$systematic_name,colnames(giSim))]
  save(giSim,file=paste0(wd,'/intermediary/giSim.RData'))
  
  rm(giSim, nUniq, uniq, diff, dropped)
  
}

####################################
# PART II. PROCESSING AND ANALYSIS #
####################################

###########################
# B-MATRIX & PERMUTATIONS #
###########################

# hotspot summary statistics 
if(do_analyses){
  
  # genes affected by hotspots
  temp <- apply(b_mod,1,function(x)sum(x!=0))
  min(temp)
  max(temp)
  median(temp)
  mean(temp)
  
  # how many genes are affected by at least one hotspot
  sum(apply(b_mod,1,function(x)sum(x!=0)) != 0)/nrow(b_mod)
  
  
  # looking at hotspots over genes
  temp <- apply(b_mod,2,function(x)sum(x!=0))
  min(temp)
  max(temp)
  median(temp)
  mean(temp)
  
}

# Figure out if there are any overlapping genes
# # crude and not optimized for bigger searches
# # there aren't enough to worry
if(check_sanity){
  
  load(paste0(wd,'/intermediary/grBodies.RData'))
  
  # find which genes overlaps each other (and themselves trivially)
  olaps <- findOverlaps(grBodies, grBodies) %>% as.data.frame
  
  # Remove overlaps of genes with themselves (duh)
  olaps <- olaps[apply(olaps, 1, function(x)x[1] != x[2]),]
  
  # nothing overlaps with anything but its neigbor (sanity check), should be 1 if true
  apply(olaps,1,function(x)x[1]-x[2]) %>% (function(x)sum(abs(x) == 1)/length(x))  
  
  # every other is redundant, e.g. gene 5 overlaps with 6 and 6 with 5
  olaps <- olaps[seq(1,nrow(olaps),2),]
  
  # plug these into the metadata
  metaNT[sort(unlist(olaps)),1:7]
  
  rm(olaps,grBodies)
  
}

# the Rcpp-style indices for all genes to be examined
numGenes <- 1:(nrow(b_mod)-1)

# this index holds all the indices which will NOT cross chromosomes
nIndex <- numGenes[!(numGenes %in% chrkey[-1])]

# get the number of genes affected by each hotspot
nAffected <- apply(b_mod,2,function(x)sum(x!=0))

# generate B-matrix permutations without omissions (~18.5 min)
load(paste0(wd,'/intermediary/b_permutations.RData'))
if(!exists('b_permutations')){
  b_permutations <- apply(b_mod,2,function(x)permute(1e5,x,nIndex))
  save(b_permutations,file=paste0(wd,'/intermediary/b_permutations.RData'))
}

if(do_analyses){
  
  # median doublet counts of permutations for each hotspot without omissions
  permMedians <- apply(b_permutations,2,median)
  
  # doublet counts for each hotspot without omissions
  b_nCounts <- apply(b_mod,2,function(x)nCount(x,nIndex))
  
  # The following is a sanity check to:
  # double check that doublet counts are consistent in R and Rcpp
  # all 102 should be the same
  # this function also allows us to see how many items uniquely contribute to doublets
  if(check_sanity){
    
    rS <- apply(b_mod,2,function(x)dUniq(x,nIndex+1))
    sum(b_nCounts == rS[1,]) %>% print()
    
    # OAF 1 targets uniquely contributing to at least 1 doublet
    rS[2,1] %>% print()
    
    rm(rS)
  }
  
  # p-values comparing real data doublets to permutation doublets
  # # any p-values appearing to be zero are actually p < number_permutations
  b_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(b_permutations),b_nCounts)
  
  allPerms <- as.data.frame(matrix(NA,nrow=6,ncol=6)) %>%
    `colnames<-`(c('Bonf_Count','Bonf_pVal','Nom_Count','Nom_pVal','Med_Count','Med_pVal')) %>%
    `rownames<-`(c('Doublets','NoDivergent','Triplets','Quadruplets','Quintuplets','Sextuplets'))
  
  # binomial tests at different significance thresholds
  allPerms['Doublets','Bonf_Count'] <- sum(b_pVals < 0.05/length(b_pVals))
  allPerms['Doublets','Bonf_pVal'] <- binom.test(sum(b_pVals < 0.05/length(b_pVals)),
                                                 length(b_pVals),0.05/length(b_pVals),alternative='greater')$p.value
  allPerms['Doublets','Nom_Count'] <- sum(b_pVals < 0.05)
  allPerms['Doublets','Nom_pVal'] <-binom.test(sum(b_pVals < 0.05),length(b_pVals),0.05,alternative='greater')$p.value
  allPerms['Doublets','Med_Count'] <- sum(b_pVals < 0.5)
  allPerms['Doublets','Med_pVal'] <- binom.test(sum(b_pVals < 0.5),length(b_pVals),0.5,alternative='greater')$p.value
  
  # an exmaination of OAF1
  nAffected[which(names(nAffected) == 'chrI:48890_C/T')]
  b_nCounts[which(names(b_nCounts) == 'chrI:48890_C/T')]
  b_pVals[which(names(b_pVals) == 'chrI:48890_C/T')]
  permMedians[which(names(permMedians) == 'chrI:48890_C/T')]
  
  rm(permMedians)
  
  # Generate Figure01 no omission permutation with visual normalization 
  if(make_plots){
    
    # min-median normalize everything
    temp <- rbind(b_nCounts,b_permutations) %>% (function(x)apply(x,2,function(y){return((y-median(y))/max(y))}))
    
    # separate the real doublets fromt the permuted doublets
    rDoublets <- temp[1,]
    
    # tabulate results with normalized doublets to help match to figure
    tempTab <- apply(b_permutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = nAffected) %>% add_column (obsDoublets = b_nCounts) %>%
      add_column(pVal = b_pVals) %>% 
      add_column(normalizedDoublets = rDoublets) %>%
      write.csv(file=paste0(wd,'/tables/','SupplTable01_01_doublets.csv'))
    rm(tempTab)
    
    # get the order of the highest real doublets with respect to permutation
    fOrder <- names(rDoublets)[order(rDoublets)]
    
    # get the permuted doublets into a plot-ready form
    pDoublets <- temp[2:ncol(temp),] %>% as.data.frame %>% 
      gather(1:length(b_nCounts),key='Hotspot',value='Doublets') %>%
      mutate(Hotspot = ordered(Hotspot, levels = fOrder))
    
    # get the real doublets into a plot-ready form
    rDoublets <- temp[1,] %>% as.data.frame %>% add_column(Hotspot=colnames(temp)) %>%
      `names<-`(c('Doublets','Hotspot')) %>%
      add_column(pVals = b_pVals) %>% 
      mutate(category=cut(pVals, breaks=c(-Inf,0.05/102,0.05,Inf), 
                          labels=c("p < 0.05/102","0.05/102 \u2264 p < 0.05","p \u2265 0.05"))) %>%
      mutate(Hotspot = ordered(Hotspot, levels = fOrder)) %>% rename(Significance = category)
    
    templot <- ggplot() + geom_boxplot(data=pDoublets,aes(x=Hotspot,y=Doublets),width=0.5,outlier.size=0.4, outlier.shape=15) +
      geom_point(data=rDoublets,aes(x=Hotspot,y=Doublets,colour=Significance),size=0.8) + 
      theme_bw() + scale_color_manual(values=viridis(6)[c(1,3,5)]) + 
      theme(axis.text.x=element_blank(),legend.position = 'bottom') +
      ylab('Normalized doublet count') + xlab('Ordered hotspots')
    
    
    ggsave(paste0(wd,'/plots/','Figure01_realVsPerm_b.pdf'),plot = templot,device = cairo_pdf,
           path = NULL,scale = 1,width = 7.08,height = 3.5,units = c("in"))
    
    rm(temp,rDoublets,fOrder,pDoublets,templot)
    
  }
  
  # Generate Figure01 no omission permutation with visual normalization 
  if(make_plots) {
    
    # min-median normalize everything
    temp <- rbind(b_nCounts,b_permutations)
    
    # separate the real doublets fromt the permuted doublets
    rDoublets <- temp[1,]
    
    # get the order of the highest real doublets with respect to permutation
    fOrder <- names(rDoublets)[order(rDoublets)]
    
    pDoublets <- temp[2:ncol(temp),] %>% as.data.frame %>% 
      gather(1:length(b_nCounts),key='Hotspot',value='Doublets') %>%
      mutate(Hotspot = ordered(Hotspot, levels = fOrder))
    
    pDoublets$Break <- pDoublets$Hotspot %in% fOrder[93:102]
    
    rDoublets <- temp[1,] %>% as.data.frame %>% add_column(Hotspot=colnames(temp)) %>%
      `names<-`(c('Doublets','Hotspot')) %>%
      add_column(pVals = b_pVals) %>% 
      mutate(category=cut(pVals, breaks=c(-Inf,0.05/102,0.05,Inf), 
                          labels=c("p < 0.05/102","0.05/102 \u2264 p < 0.05","p \u2265 0.05"))) %>%
      mutate(Hotspot = ordered(Hotspot, levels = fOrder)) %>% rename(Significance = category)
    
    rDoublets$Break <- rDoublets$Hotspot %in% fOrder[93:102]
    
    mapply(function(a,b,c){
      templot <- ggplot() + geom_boxplot(data=pDoublets[pDoublets$Break==a,],aes(x=Hotspot,y=Doublets),width=0.5,outlier.size=0.4, outlier.shape=15) +
        geom_point(data=rDoublets[rDoublets$Break==a,],aes(x=Hotspot,y=Doublets,colour=Significance),size=0.8) + 
        theme_bw() + scale_color_manual(values=viridis(6)[c(1,3,5)]) + 
        theme(axis.text.x=element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.position = 'bottom') +
        ylab('Doublet count') + xlab('Ordered hotspots')
      
      ggsave(b,plot = templot,device = cairo_pdf,
             path = NULL,scale = 1,width = c,height = 7,units = c("in"))
      
    },c(F,T),
    c(paste0(wd,'/plots/','SupplFigure01_Left.pdf'),paste0(wd,'/plots/','SupplFigure01_Right.pdf')),
    c(5.5,1.5)
    )
    
    rm(temp,rDoublets,fOrder,pDoublets)
  }
  
  rm(b_nCounts,b_pVals)
  
}

rm(b_permutations)

# Rcpp-style indices for gene pairs which are not divergent and do not span a chromosome break
nDiverge <- numGenes[!(numGenes %in% sort(unique(c(chrkey,which(pOrient == 'd')))))]

# GENERATE PERMUTATIONS MINUS THE DIVERGENT GENE PAIRS (~17 min)
load(paste0(wd,'/intermediary/nDiverge_permutations.RData'))
if(!exists('nDiverge_permutations')){
  nDiverge_permutations <- apply(b_mod,2,function(x)permute(1e5,x,nDiverge))
  save(nDiverge_permutations,file=paste0(wd,'/intermediary/nDiverge_permutations.RData'))
}

# doublet counts with the omission of divergent genes
if(do_analyses){
  
  # no divergent doublet counts for each hotspot
  nDiverge_counts <- apply(b_mod,2,function(x)nCount(x,nDiverge))
  
  # no divergent doublet p-values
  nDiverge_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(nDiverge_permutations),nDiverge_counts)
  
  # binomial tests at different significance thresholds
  allPerms['NoDivergent','Bonf_Count'] <- sum(nDiverge_pVals < 0.05/length(nDiverge_pVals))
  allPerms['NoDivergent','Bonf_pVal'] <- binom.test(sum(nDiverge_pVals < 0.05/length(nDiverge_pVals)),
                                                    length(nDiverge_pVals),0.05/length(nDiverge_pVals),alternative='greater')$p.value
  allPerms['NoDivergent','Nom_Count'] <- sum(nDiverge_pVals < 0.05)
  allPerms['NoDivergent','Nom_pVal'] <- binom.test(sum(nDiverge_pVals < 0.05),length(nDiverge_pVals),0.05,alternative='greater')$p.value
  allPerms['NoDivergent','Med_Count'] <- sum(nDiverge_pVals < 0.5)
  allPerms['NoDivergent','Med_pVal'] <-binom.test(sum(nDiverge_pVals < 0.5),length(nDiverge_pVals),0.5,alternative='greater')$p.value
  
  # tabulate results
  if(make_plots){
    apply(nDiverge_permutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = nAffected) %>% add_column (obsDoublets = nDiverge_counts) %>%
      add_column(pVal = nDiverge_pVals) %>% write.csv(file=paste0(wd,'/tables/','SupplTable01_02_noDiverge.csv'))
  }
  
  rm(nDiverge_counts, nDiverge_pVals)
  
}

rm(nDiverge_permutations, nDiverge)

# this index holds all the triplet indices which will NOT cross chromosomes
triIndex <- numGenes[!(numGenes %in% sort(c(chrkey,chrkey[-1]+1)))]

# generate triplet B-matrix permutations without omissions (~21 min)
load(paste0(wd,'/intermediary/triPermutations.RData'))
if(!exists('triPermutations')){
  triPermutations <- apply(b_mod,2,function(x)triPerm(1e5,x,triIndex))
  save(triPermutations,file=paste0(wd,'/intermediary/triPermutations.RData'))
}

if(do_analyses){
  
  # triplet counts for each hotspot without omissions
  triCounts <- apply(b_mod,2,function(x)triCount(x,triIndex))
  
  # triplet p-values
  tri_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(triPermutations),triCounts)
  
  # binomial tests at different significance thresholds
  allPerms['Triplets','Bonf_Count'] <- sum(tri_pVals < 0.05/length(tri_pVals))
  allPerms['Triplets','Bonf_pVal'] <- binom.test(sum(tri_pVals < 0.05/length(tri_pVals)),
                                                 length(tri_pVals),0.05/length(tri_pVals),alternative='greater')$p.value
  allPerms['Triplets','Nom_Count'] <- sum(tri_pVals < 0.05)
  allPerms['Triplets','Nom_pVal'] <- binom.test(sum(tri_pVals < 0.05),length(tri_pVals),0.05,alternative='greater')$p.value
  allPerms['Triplets','Med_Count'] <- sum(tri_pVals < 0.5)
  allPerms['Triplets','Med_pVal'] <-binom.test(sum(tri_pVals < 0.5),length(tri_pVals),0.5,alternative='greater')$p.value
  
  # tabulate results
  if(make_plots){
    tempTab <- apply(triPermutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = nAffected) %>% add_column (obsTriplets = triCounts) %>%
      add_column(pVal = tri_pVals) %>% write.csv(file=paste0(wd,'/tables/','SupplTable01_03_triplets.csv'))
    rm(tempTab)
  }
  
  rm(triCounts, tri_pVals)
  
}

rm(triPermutations, triIndex)

# QUADRUPLET PERMUTATIONS

# this index holds all the indices which will NOT cross chromosomes
quadIndex <- numGenes[!(numGenes %in% sort(c(chrkey,chrkey+1,chrkey[-1]+2)))]

# generate quadruplet B-matrix permutations without omissions (~19 min)
load(paste0(wd,'/intermediary/quadPermutations.RData'))
if(!exists('quadPermutations')){
  quadPermutations <- apply(b_mod,2,function(x)quadPerm(1e5,x,quadIndex))
  save(quadPermutations,file=paste0(wd,'/intermediary/quadPermutations.RData'))
}

# analyze quadruplet permutations
if(do_analyses){
  
  # quadruplet counts for each hotspot without omissions
  quadCounts <- apply(b_mod,2,function(x)quadCount(x,quadIndex))
  
  # quadruplet p-values
  quad_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(quadPermutations),quadCounts)
  
  # binomial tests at different significance thresholds
  allPerms['Quadruplets','Bonf_Count'] <- sum(quad_pVals < 0.05/length(quad_pVals))
  allPerms['Quadruplets','Bonf_pVal'] <- binom.test(sum(quad_pVals < 0.05/length(quad_pVals)),
                                                    length(quad_pVals),0.05/length(quad_pVals),alternative='greater')$p.value
  allPerms['Quadruplets','Nom_Count'] <- sum(quad_pVals < 0.05)
  allPerms['Quadruplets','Nom_pVal'] <- binom.test(sum(quad_pVals < 0.05),length(quad_pVals),0.05,alternative='greater')$p.value
  allPerms['Quadruplets','Med_Count'] <- sum(quad_pVals < 0.5)
  allPerms['Quadruplets','Med_pVal'] <-binom.test(sum(quad_pVals < 0.5),length(quad_pVals),0.5,alternative='greater')$p.value
  
  # tabulate results
  if(make_plots){
    tempTab <- apply(quadPermutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = nAffected) %>% add_column (obsQuadruplets = quadCounts) %>%
      add_column(pVal = quad_pVals) %>% write.csv(file=paste0(wd,'/tables/','SupplTable01_04_quadruplets.csv'))
    rm(tempTab)
  }
  
  rm(quadCounts, quad_pVals)
  
}

rm(quadPermutations, quadIndex)


# this index holds all the indices which will NOT cross chromosomes
quinIndex <- numGenes[!(numGenes %in% sort(c(chrkey,chrkey+1,chrkey+2,chrkey[-1]+3)))]

# generate quintuplet B-matrix permutations without omissions (~19 min)
load(paste0(wd,'/intermediary/quinPermutations.RData'))
if(!exists('quinPermutations')){
  quinPermutations <- apply(b_mod,2,function(x)quinPerm(1e5,x,quinIndex))
  save(quinPermutations,file=paste0(wd,'/intermediary/quinPermutations.RData'))
}

# analyze quintuplet permutations
if(do_analyses){
  
  # doublet counts for each hotspot without omissions
  quinCounts <- apply(b_mod,2,function(x)quinCount(x,quinIndex))
  
  # quintuplet p-values
  quin_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(quinPermutations),quinCounts)
  
  # binomial tests at different significance thresholds
  allPerms['Quintuplets','Bonf_Count'] <- sum(quin_pVals < 0.05/length(quin_pVals))
  allPerms['Quintuplets','Bonf_pVal'] <- binom.test(sum(quin_pVals < 0.05/length(quin_pVals)),
                                                    length(quin_pVals),0.05/length(quin_pVals),alternative='greater')$p.value
  allPerms['Quintuplets','Nom_Count'] <- sum(quin_pVals < 0.05)
  allPerms['Quintuplets','Nom_pVal'] <- binom.test(sum(quin_pVals < 0.05),length(quin_pVals),0.05,alternative='greater')$p.value
  allPerms['Quintuplets','Med_Count'] <- sum(quin_pVals < 0.5)
  allPerms['Quintuplets','Med_pVal'] <-binom.test(sum(quin_pVals < 0.5),length(quin_pVals),0.5,alternative='greater')$p.value
  
  # tabulate results
  if(make_plots){
    tempTab <- apply(quinPermutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = nAffected) %>% add_column (obsQuintuplet = quinCounts) %>%
      add_column(pVal = quin_pVals) %>% write.csv(file=paste0(wd,'/tables/','SupplTable01_05_quintuplets.csv'))
    rm(tempTab)
  }
  
  rm(quinCounts, quin_pVals)
  
}

rm(quinPermutations, quinIndex)

# this index holds all the indices which will NOT cross chromosomes
sextIndex <- numGenes[!(numGenes %in% sort(c(chrkey,chrkey+1,chrkey+2,chrkey+3,chrkey[-1]+4)))]

# generate sexttuplet B-matrix permutations without omissions (~19 min)
load(paste0(wd,'/intermediary/sextPermutations.RData'))
if(!exists('sextPermutations')){
  sextPermutations <- apply(b_mod,2,function(x)sextPerm(1e5,x,sextIndex))
  save(sextPermutations,file=paste0(wd,'/intermediary/sextPermutations.RData'))
}

# analyze sexttuplet permutations
if(do_analyses){
  
  # doublet counts for each hotspot without omissions
  sextCounts <- apply(b_mod,2,function(x)sextCount(x,sextIndex))
  
  # sexttuplet p-values
  sext_pVals <- mapply(function(x,y)sum(x >= y)/length(x),as.data.frame(sextPermutations),sextCounts)
  
  # binomial tests at different significance thresholds
  allPerms['Sextuplets','Bonf_Count'] <- sum(sext_pVals < 0.05/length(sext_pVals))
  allPerms['Sextuplets','Bonf_pVal'] <- binom.test(sum(sext_pVals < 0.05/length(sext_pVals)),
                                                   length(sext_pVals),0.05/length(sext_pVals),alternative='greater')$p.value
  allPerms['Sextuplets','Nom_Count'] <- sum(sext_pVals < 0.05)
  allPerms['Sextuplets','Nom_pVal'] <- binom.test(sum(sext_pVals < 0.05),length(sext_pVals),0.05,alternative='greater')$p.value
  allPerms['Sextuplets','Med_Count'] <- sum(sext_pVals < 0.5)
  allPerms['Sextuplets','Med_pVal'] <-binom.test(sum(sext_pVals < 0.5),length(sext_pVals),0.5,alternative='greater')$p.value
  
  # tabulate results
  if(make_plots){
    tempTab <- apply(sextPermutations,2,summary) %>% t %>% as.data.frame() %>% 
      add_column(nAffected = nAffected) %>% add_column (obsSextuplet = sextCounts) %>%
      add_column(pVal = sext_pVals) %>% write.csv(file=paste0(wd,'/tables/','SupplTable01_06_sextuplets.csv'))
    rm(tempTab)
  }
  
  rm(sextCounts, sext_pVals)
  
}

rm(sextPermutations, sextIndex)



# finish the table with all the binomial tests
if(make_plots){
  allPerms %>% `colnames<-`(c('p < 0.05/102 Counts','p < 0.05/102 Binomial p-Vals',
                              'p < 0.05 Counts','p < 0.05 Binomial p-Vals',
                              'p < 0.5 Counts','p < 0.5 Binomial p-Vals')) %>%
    write.csv(file=paste0(wd,'/tables/','Table01_permutationSummaries.csv'))
}

rm(nAffected, nIndex, allPerms)

###################################
# FREQUENCY AND DOUBLETS MATRICES #
###################################

# GENERATE TERNARY PARTNER SCORE
load(paste0(wd,'/intermediary/bsmpl.RData'))
if(!exists('bsmpl')){
  bsmpl <- bFrequency(as.matrix(b_mod))
  rownames(bsmpl) <- rownames(b_mod)
  colnames(bsmpl) <- rownames(b_mod)
  save(bsmpl,file=paste0(wd,'/intermediary/bsmpl.RData'))
}

# GENERATE DOUBLET PARTNER SCORE
load(paste0(wd,"/intermediary/bplus.RData"))
if(!exists('bplus')){
  bplus <- bPlus (as.matrix(b_mod))
  rownames(bplus) <- rownames(b_mod)
  colnames(bplus) <- rownames(b_mod)
  save(bplus,file=paste0(wd,"/intermediary/bplus.RData"))
}


# GIVE US THE PROPORTIONS OF POSITIVE, NEGATIVE, AND ZERO-VALUE PAIRS IN FREQUENCY MATRIX
if(do_analyses){
  
  # test for different proportions of positive values in neighbors and non-neighbors
  bnPosProps <- sapply(list(bsmpl,bplus),function(x)prop.test(
    x = c(sum(neighbor_grab(x,chrkey)>0,na.rm=T),sum(nonneighbor_grab(x,chrkey)>0,na.rm=T)),
    n = c(sum(!is.na(neighbor_grab(x,chrkey))),sum(!is.na(nonneighbor_grab(x,chrkey)))))$p.value)
  names(bnPosProps) <- c('Frequency','Doublet')
  bnPosProps
  
  # test for different proportions of negative values in neighbors and non-neighbors (bplus has no neg)
  freqNegProps <-prop.test(
    x = c(sum(neighbor_grab(bsmpl,chrkey)<0,na.rm=T),sum(nonneighbor_grab(bsmpl,chrkey)<0,na.rm=T)),
    n = c(sum(!is.na(neighbor_grab(bsmpl,chrkey))),sum(!is.na(nonneighbor_grab(bsmpl,chrkey)))))$p.value
  freqNegProps
  
  # test for different proportions of zero values in neighbors and non-neighbors
  bnZeroProps <- sapply(list(bsmpl,bplus),function(x)prop.test(
    x = c(sum(neighbor_grab(x,chrkey)==0,na.rm=T),sum(nonneighbor_grab(x,chrkey)==0,na.rm=T)),
    n = c(sum(!is.na(neighbor_grab(x,chrkey))),sum(!is.na(nonneighbor_grab(x,chrkey)))))$p.value)
  names(bnPosProps) <- c('Frequency','Doublet')
  bnZeroProps
  
  # wilcoxons (two-tailed)
  bnPvals <- sapply(list(bsmpl,bplus),function(x)
    wilcox.test(neighbor_grab(x,chrkey),nonneighbor_grab(x,chrkey))$p.value)
  names(bnPvals) <- c('Frequency','Doublet')
  bnPvals
  
  # get mean values
  bnMeans <- sapply(list(bsmpl,bplus),function(x)
    c(mean(neighbor_grab(x,chrkey)),mean(nonneighbor_grab(x,chrkey))))
  colnames(bnMeans) <- c('Frequency','Doublet')
  rownames(bnMeans) <- c('Neighbor','Non-Neighbor')
  bnMeans
  
  rm(bnPosProps,freqNegProps,bnZeroProps,bnPvals,bnMeans)
  
  # plot proportions of directions by position from frequency matrix
  if(make_plots){
    
    templot <- lapply(c(neighbor_grab,nonneighbor_grab),function(x)x(bsmpl,chrkey)) %>% 
      sapply(function(x){
        c(sum(x>0),sum(x==0),sum(x<0))
      }) %>% `colnames<-`(c('Adjacent','Non-adjacent')) %>%
      apply(2,function(x)x/sum(x)) %>% as.data.frame() %>% 
      add_column(Direction = c('Same Direction','Zero','Different Direction')) %>%
      mutate(Direction = ordered(Direction, levels = c('Same Direction','Different Direction','Zero'))) %>%
      gather(key='Position',value='Proportion',1:2) %>% 
      mutate(Position = ordered(Position, levels = c('Non-adjacent','Adjacent'))) %>%
      ggplot(aes(x=Position,y=Proportion,fill=Position)) +
      geom_bar(stat='Identity',  colour="black") + theme_bw() + ylim(c(0,0.6)) +
      facet_wrap('Direction') + scale_fill_grey() +
      theme(legend.position = "none", axis.title.x=element_blank(), 
            axis.text.x.bottom = element_text(vjust = 0.9,angle = 45, hjust = 1))
    
    ggsave(paste0(wd,'/plots/','Figure03_bProps.pdf'),plot = templot,device = NULL,
           path = NULL,scale = 1,width = 7.08,height = 3,units = c("in"))
    
    rm(templot) 
  }
  
}

#######################
# WITHIN COEXPRESSION #
#######################

# Build correlation matrices if not already done (~50 min)
cor_names <- sort(c('Fleming02','Hughes00', 'Klevecz04','Brem05', 'Knijnenburg09', 'Simola10', 'Lenstra11',
                    'Sameith15','Schurch16','Myers19','Albert18'))
cor_mats <- rep(list(NA),length(cor_names))
names(cor_mats) <- cor_names
for(i in 1:length(cor_names)){
  print(cor_names[i])
  cor_mats[[i]] <- spear_cor(paste0(wd,'/intermediary/'),cor_names[i],r=T)
}
rm(cor_names)

# example correlation for figure 2
if(make_plots){
  
  # get a subset of the Albert 2018 coexpression matrix from  YDL082W  to YDL078C
  miniCor <- melt(cor_mats$Albert18[674:677,674:677],na.rm=T) %>% `names<-`(c('Var1','Var2','Rho'))
  
  diverge <- scale_fill_gradient2(
    low = "#440154FF",
    mid = "white",
    high = "#2A788EFF",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1, 1)
  )
  
  templot <- ggplot(data=miniCor, aes(Var1, Var2, fill= Rho)) + 
    geom_tile() +
    coord_equal() + 
    geom_text(aes(label = round(Rho,2)),colour = "black",size=3) + 
    diverge + 
    labs(x = "Gene 1", y = "Gene 2") + 
    theme_minimal() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x.bottom = element_text(vjust = 1.1,angle = 90, hjust = 1),
          legend.position='bottom',
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + scale_y_discrete(position = 'right')
  
  ggsave(paste0(wd,'/plots/','Figure02_miniCor.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 3.2,height = 3.2,units = c("in"))
  
  rm(templot, Albert18, temp, miniCor)
  
}

# view distributions (~36 min)
if(make_plots){
  templot <- sapply(cor_mats,function(x)tri(x)) %>% as.data.frame %>% 
    (function(x)gather(x,1:length(x),key='Source',value='Rho')) %>% 
    ggplot(aes(x=Source,y=Rho)) + geom_violin(aes(fill=Source)) +
    geom_boxplot(width=0.1, outlier.shape=NA) + scale_fill_viridis_d() + 
    theme_bw() + theme(axis.text.x.bottom = element_text(vjust = 0.9,angle = 45, hjust = 1))
  
  ggsave(paste0(wd,'/plots/','SupplFigure02_corDistribs.pdf'),device = NULL,
         path = NULL,scale = 1,width = 6.5,height = 4,units = c("in"))
  
  dev.off()
}

# look at proportion of positive correlations
load(paste0(wd,'/intermediary/','nbrProps.RData'))
if(!exists('nbrProps')){
  nbrProps <- sapply(cor_mats,function(x)prop.test(
    x = c(sum(neighbor_grab(x,chrkey)>0,na.rm=T),sum(nonneighbor_grab(x,chrkey)>0,na.rm=T)),
    n = c(sum(!is.na(neighbor_grab(x,chrkey))),sum(!is.na(nonneighbor_grab(x,chrkey))))
  )$p.value)
  save(nbrProps,file=paste0(wd,'/intermediary/','nbrProps.RData'))
}

# complete neighbor testing looking at the follwing (~31 min):
# if neighbors or more or less than non-neighborings
# for each dataset
# for negative, all, and positive correlations
load(paste0(wd,'/intermediary/','nbrWilcoxes.RData'))
if(!exists('nbrWilcoxes')){
  
  nbrWilcoxes <- lapply(c('greater','less'),function(x)
    sapply(c(-1:1),function(y)
      sapply(cor_mats,function(z){
        
        # so you have some semblance of when you will finish
        
        wilcox.test(
          (neighbor_grab(z,chrkey) %>% matsign(y)),
          (nonneighbor_grab(z,chrkey) %>% matsign(y)),
          alternative = x
        )$p.value %>% return()
        
      }
      )))
  
  # assign names to prevent confusion
  names(nbrWilcoxes) <- c('Greater','Less')
  colnames(nbrWilcoxes$Greater) <- c('Negative','All','Positive')
  colnames(nbrWilcoxes$Less) <- c('Negative','All','Positive')
  nbrWilcoxes <- lapply(nbrWilcoxes,as.data.frame)
  
  # save ye test
  save(nbrWilcoxes,file=paste0(wd,'/intermediary/','nbrWilcoxes.RData'))
  
}

if(do_analyses){
  
  # are neighbors always more positive than expected?
  apply(nbrWilcoxes$Greater,2,function(x)binom.test(sum(x < 0.05),length(x)))
  
  # weakest of comparisons by group
  max(nbrWilcoxes$Greater$All)
  max(nbrWilcoxes$Greater$Negative)
  max(nbrWilcoxes$Greater$Positive)
  max(nbrProps) 
  
}

rm(nbrProps,nbrWilcoxes)

# display elevated correlation of neighbors for figure 2B
if(make_plots){
  
  # split by position and remove matrices that are not to be used
  pSplit <- sapply(cor_mats,function(x){
    temp <- list(neighbor_grab(x,chrkey),nonneighbor_grab(x,chrkey))
    return(temp)
  }) 
  
  # first get the median for negative, positive, and unfiltered rhos by neighbor or non-neighbor
  rhoSplit <- lapply(c(-1:1),function(y)apply(pSplit,c(1,2),function(x)median(matsign(unlist(x),y),na.rm=T))) %>% 
    lapply(as.data.frame) %>% # get into dataframe format
    lapply(function(x)add_column(x,Position=c('Adjacent','Non-adjacent'))) %>%
    lapply(function(x)gather(x,key='Source',value='Value',1:(ncol(x)-1))) %>%
    (function(a)mapply(function(x,y)add_column(x,Metric=y),a,c('Negative','All','Positive'))) %>%
    apply(2,cbind) %>%
    lapply(function(x)do.call(cbind,x)) %>%
    (function(x)do.call(rbind,x)) %>%
    as.data.frame %>%
    `names<-`(c('Position','Source','Value','Metric')) # assure your names are correct
  
  # make sure some values are not factors
  rhoSplit$Value <- rhoSplit$Value %>% as.character() %>% as.numeric()
  
  # split your proportions
  propSplit <- apply(pSplit,c(1,2),function(x)sum(unlist(x) > 0,na.rm=T)/sum(!is.na(unlist(x)))) %>% as.data.frame %>%
    add_column(Position=c('Adjacent','Non-adjacent')) %>%
    (function(x)gather(x,key='Source',value='Value',1:(ncol(x)-1))) %>% 
    add_column(Metric='Proportion')
  
  # combine the proportions are rhos into yourfinal plot-ready item
  allSplit <- rbind(rhoSplit,propSplit) %>% 
    mutate(Position = ordered(Position, levels = c('Non-adjacent','Adjacent'))) %>%
    mutate(Metric = ordered(Metric, levels = c('All','Positive','Negative','Proportion')))
  
  # remove the processing items used to create the above
  rm(pSplit,rhoSplit,propSplit)
  
  # shift the plot limits to make room for p-values
  limTemp <- data.frame(Metric = c("All", "All", "Positive", "Positive", 
                                   "Negative", "Negative","Proportion","Proportion"), 
                        Position = 'Adjacent', Value = c(-0.02, 0.15, 0.09, 0.5, -0.4, 0, 0.48, 0.7)) %>%
    mutate(Metric = ordered(Metric, levels = c('All','Positive','Negative','Proportion')))
  
  # plot your data
  templot <- ggplot(allSplit) + 
    geom_boxplot(aes(y=Value,x=Position),fill='grey',outlier.shape=NA) + 
    geom_line(aes(y=Value,x=Position,group=Source, colour=Source)) +
    geom_point(aes(y=Value,x=Position,colour=factor(Source)),size=2) +
    geom_blank(data = limTemp,aes(y=Value,x=Position)) +
    theme_bw() + theme(axis.title.x=element_blank(),
                       legend.title = element_blank(),
                       axis.text.x = element_text(angle=45,hjust=1.1,vjust=1)) + 
    facet_wrap(~Metric, scales='free_y') + 
    scale_colour_viridis_d()
  
  # save your plot
  ggsave(paste0(wd,'/plots/','Figure02_nbrEffect.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 4,height = 5,units = c("in"), useDingbats=F)
  
  # remove extraneous objects
  rm(allSplit,templot,limTemp)
  
}

# number of significant correlations within matrices and across them
if(do_analyses){
  
  corSig <- as.data.frame(matrix(NA,ncol=3,nrow=length(cor_mats)))
  colnames(corSig) <- c('Nominal','Bonferroni','Total')
  rownames(corSig) <- names(cor_mats)
  
  for(i in names(cor_mats)){
    load(paste0(wd,'/intermediary/',i,'_cor.RData'))
    temp <- get(paste0(i,'_cor'))$P %>% tri()
    rm(list = paste0(i,'_cor'))
    corSig[i,'Nominal'] <- sum(temp < 0.05,na.rm=T)
    corSig[i,'Bonferroni'] <- sum(temp < 0.05/sum(!is.na(temp)),na.rm=T)
    corSig[i,'Total'] <- sum(!is.na(temp))
    rm(temp)
  }
  
  corSig$fracNom <- corSig$Nominal/corSig$Total
  corSig$fracBon <- corSig$Bonferroni/corSig$Total
  
  print(corSig)
  
  print(sum(corSig$fracNom > 0.5))
  print(nrow(corSig))
  print(sum(corSig$Nominal)/sum(corSig$Total))
  print(sum(corSig$Bonferroni)/sum(corSig$Total))
  
  rm(corSig)
  
}

#######################
# ACROSS COEXPRESSION #
#######################

# pairwise spearman correlate each correlation matrix (~8 min)
load(paste0(wd,'/intermediary/','metaCors.RData'))
if(!(exists('metaCors'))){
  metaCors <- rep(list(NA),2) %>% `names<-` (c('nbr','non'))
  metaCors$nbr <- sapply(cor_mats,function(x)neighbor_grab(x,chrkey)) %>% rcorr(type='spearman') 
  metaCors$non <- sapply(cor_mats,function(x)nonneighbor_grab(x,chrkey)) %>% rcorr(type='spearman')
  save(metaCors,file=paste0(wd,'/intermediary/','metaCors.RData'))
}

# rudimentary analysis of meta correlations
if(do_analyses){
  
  # compare the neighbor structures to the non-neighbor structures
  wilcox.test(metaCors$nbr$r,metaCors$non$r,paired=T)
  
  # weakest p-value
  max(metaCors$nbr$P,na.rm=T) %>% print()
  max(metaCors$non$P,na.rm=T) %>% print()
  
  # medians
  median(metaCors$nbr$r) %>% print()
  median(metaCors$non$r) %>% print()
}

# create labelled heatmap of the correlation between coexpression matrices 
if(make_plots){
  
  # turn one matrix into an upper and the other into a lower triangle dataframe
  
  dfCor1 <- melt(metaCors$non$r %>% 
                   (function(x){
                     x[lower.tri(x,diag=T)]<- NA
                     return(x)
                   }), na.rm=T)
  
  dfCor2 <- melt(metaCors$nbr$r %>% 
                   (function(x){
                     x[upper.tri(x,diag=T)]<- NA
                     return(x)
                   }), na.rm=T)
  
  # merge the upper and lower triangles into a single dataframe and name Rho correctly
  dfCor <- rbind(dfCor1, dfCor2) %>% rename(Rho=value)
  
  # remove the objects used to generate dfCor
  rm(dfCor1,dfCor2)
  
  
  # Correct your order
  dfCor$Var1 <- ordered(dfCor$Var1, levels = sort(as.character(unique(dfCor$Var1))))
  dfCor$Var2 <- ordered(dfCor$Var2, levels = sort(as.character(unique(dfCor$Var2))))
  
  diverge <- scale_fill_gradient2(
    low = "#440154FF",
    mid = "white",
    high = "#2A788EFF",
    midpoint = 0,
    na.value = "grey50",
    limits=c(-1, 1)
  )
  
  # MAIN PLOT
  
  templot <- ggplot(data=dfCor, aes(Var1, Var2, fill= Rho)) + 
    geom_tile() +
    coord_equal() + 
    geom_text(aes(label = round(Rho,2)),colour = "black",size=3) + 
    diverge + 
    labs(x = "Adjacent gene pairs", y = "Non-adjacent gene pairs") + 
    theme_minimal() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x.bottom = element_text(vjust = 1.1,angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 45),
          legend.position='left')
  
  ggsave(paste0(wd,'/plots/','Figure02_metaCor.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 5,height = 5,units = c("in"))
  
  rm(templot,dfCor)
  
}

#####################################
# RELATING HOTSPOTS TO COEXPRESSION #
#####################################

# compare coexpression matrices to magnitude and frequency matrices (~10 min)
load(paste0(wd,'/intermediary/','bTests.RData'))
if(!exists('bTests')){
  
  bTests <- rep(list(NA),2)
  names(bTests) <- c('fNon','fNbr')
  
  # frequency non-neighbor
  bTests$fNon <- lapply(cor_mats,function(x)
    cor.test(nonneighbor_grab(x,chrkey),nonneighbor_grab(bsmpl,chrkey),method='spearman'))
  
  # frequency neighbors
  bTests$fNbr <- lapply(cor_mats,function(x)
    cor.test(neighbor_grab(x,chrkey),neighbor_grab(bsmpl,chrkey),method='spearman'))
  
  save(bTests,file=paste0(wd,'/intermediary/','bTests.RData'))
  
}

# get the p-values for your correlations
if(do_analyses){
  
  # get your p- and rho-values from the spearman correlations of frequency and coexpression
  bTestP <- sapply(bTests,function(x)sapply(x,function(y)y$p.value)) %>% as.data.frame()
  bTestR <- sapply(bTests,function(x)sapply(x,function(y)y$estimate)) %>% as.data.frame() 
  
  # compare spearman correlations for neighbors and non-neighbors
  print(bTestR)
  print(bTestP)
  median(bTestR$fNon) %>% print()
  median(bTestR$fNbr) %>% print()
  wilcox.test(bTestR$fNbr,bTestR$fNon,paired=T) %>% print()
  
}

rm(bTestP, bTestR, bTests)

# understand if strong correlations drive the relationship (~21 min)
load(paste0(wd,'/intermediary/','bcor_quant.RData'))
if(!exists('bcor_quant')){
  
  # a quantile based function for analyzing relationship in slices
  gated_bcor <- function(a,b,i) {
    # get the triangles of a coexpression matrix
    x <- tri(a)
    # get the triangle of your frequency matrix
    y <- tri(b)
    # step up by quantile
    z <-  quantile(abs(x), seq(from=0,to=(1-i),by=i), na.rm=T)
    # gate what is compared by how far from 0 coexpression values are
    r <- lapply(z,function(t)cor.test(y[abs(x) > t & abs(x)< (t+i)],x[abs(x) > t & abs(x) < (t+i)],method='spearman'))
    return(r)
  }
  
  # set step interval
  interval = 0.05
  
  # get gated correlation
  bcor_quant <- sapply(cor_mats,function(x)gated_bcor(x,bsmpl,interval)) %>%
    `rownames<-`(seq(from=0,to=1-interval,by=interval))
  
  # save and remove unnecessary objects
  save(bcor_quant,file=paste0(wd,'/intermediary/','bcor_quant.RData'))
  rm(gated_bcor,interval)
}

# if you want to extract something
bCorR <- apply(bcor_quant,c(1,2),function(x)x[[1]]$estimate)
bCorP <- apply(bcor_quant,c(1,2),function(x)x[[1]]$p.value)


if(do_analyses){
  
  # understand how much the strength of the relationship correlates with quantile
  apply(bCorR,2,function(y)
    cor.test(y,as.numeric(rownames(bCorR)),type='spearman')
  )
  
  # weakest rho value for frequency
  min(bCorR)
  
  # weakest p-value for frequency
  max(bCorP)
  
}

# graph the relationship of agreement to coexpression strength
if(make_plots){
  
  # prepare data for plotting (before I knew how to use dplyr)
  temp <- bCorR %>% (function(x){
    temp <- as.data.frame(c(x))
    names(temp) <- "Rho"
    temp$Source <- c(sapply(colnames(x),function(y)rep(y,nrow(x))))
    temp$Cutoff <- c(rep(rownames(x),ncol(x)))
    return(temp)
  }) %>% as.data.frame()
  
  # plot data
  templot <- ggplot(temp,aes(y=Rho,x=Cutoff)) + geom_boxplot(fill='grey',outlier.shape=NA) +
    theme_bw() + geom_line(aes(group=factor(Source), colour=Source),show.legend=F) + 
    theme(axis.text.x = element_text(angle = 90),legend.title = element_blank()) +
    xlab('Quantile Strength Bin (n to n + 0.05)') + scale_colour_viridis_d() 
  
  ggsave(paste0(wd,'/plots/','Figure03_agreeStrength.pdf'),plot = templot,device = NULL,
         path = NULL,scale = 1,width = 3,height = 3,units = c("in"))
  
  rm(templot)
  
}

rm(cor_mats, bcor_quant, bCorR, bCorP, bsmpl)

#########################
# INDIVIDUAL MECHANISMS #
#########################
if(do_analyses){
  
  # analyze TF similarity 
  load(paste0(wd,'/intermediary/','jacTF.RData'))
  wilcox.test(neighbor_grab(jacTF,chrkey),nonneighbor_grab(jacTF,chrkey))
  median(neighbor_grab(jacTF,chrkey),na.rm=T)
  median(nonneighbor_grab(jacTF,chrkey),na.rm=T)
  cor.test(nonneighbor_grab(bplus,chrkey),nonneighbor_grab(jacTF,chrkey),method='spearman')
  cor.test(neighbor_grab(bplus,chrkey),neighbor_grab(jacTF,chrkey),method='spearman')
  
  
  # analyze chromatin baseline similarity 
  load(paste0(wd,'/intermediary/baseWeiner.RData'))
  cor.test(baseWeiner,neighbor_grab(bplus,chrkey),method='spearman')
  
  # analyze chromatin change similarity 
  load(paste0(wd,'/intermediary/deltaWeiner.RData'))
  cor.test(deltaWeiner,neighbor_grab(bplus,chrkey),method='spearman')
  
  # analyze pair orienation and doublets
  cor.test(as.numeric(pOrient[-chrkey] == 'd'),neighbor_grab(bplus,chrkey),method='spearman')
  cor.test(as.numeric(pOrient[-chrkey] == 't'),neighbor_grab(bplus,chrkey),method='spearman')
  cor.test(as.numeric(pOrient[-chrkey] == 'c'),neighbor_grab(bplus,chrkey),method='spearman')
  
  scoresD <- neighbor_grab(bplus,chrkey)[pOrient[-chrkey] == 'd']
  scoresT <- neighbor_grab(bplus,chrkey)[pOrient[-chrkey] == 't']
  scoresC <- neighbor_grab(bplus,chrkey)[pOrient[-chrkey] == 'c']
  
  wilcox.test(scoresD,scoresT)
  wilcox.test(scoresD,scoresC)
  wilcox.test(scoresT,scoresC)
  
  Anova(lm(neighbor_grab(bplus,chrkey) ~ pOrient[-chrkey]),type='III')
  
  rm(scoresD,scoresT,scoresC)
  
  # analyze intragene distance and doublets
  summary(neighbor_grab(intrachr_dist,chrkey))
  cor.test(neighbor_grab(intrachr_dist,chrkey),neighbor_grab(bplus,chrkey),method='spearman')
  
  # semantic similarity from Gene Ontology 
  load(paste0(wd,'/intermediary/','ssWang.RData'))
  median(neighbor_grab(ssWang,chrkey),na.rm=T)
  median(nonneighbor_grab(ssWang,chrkey),na.rm=T)
  wilcox.test(neighbor_grab(ssWang,chrkey),nonneighbor_grab(ssWang,chrkey))
  
  # genetic interaction similarity 
  load(paste0(wd,'/intermediary/','giSim.RData'))
  median(neighbor_grab(giSim,chrkey),na.rm=T)
  median(nonneighbor_grab(giSim,chrkey),na.rm=T)
  wilcox.test(neighbor_grab(giSim,chrkey),nonneighbor_grab(giSim,chrkey))
  
  # are predictors colinear with distance
  cor.test(as.numeric(pOrient[-chrkey] == 'd'),neighbor_grab(intrachr_dist,chrkey),method='spearman')
  cor.test(as.numeric(pOrient[-chrkey] == 't'),neighbor_grab(intrachr_dist,chrkey),method='spearman')
  cor.test(as.numeric(pOrient[-chrkey] == 'c'),neighbor_grab(intrachr_dist,chrkey),method='spearman')
  cor.test(neighbor_grab(jacTF,chrkey),neighbor_grab(intrachr_dist,chrkey),method='spearman')
  cor.test(baseWeiner,neighbor_grab(intrachr_dist,chrkey),method='spearman')
  cor.test(deltaWeiner,neighbor_grab(intrachr_dist,chrkey),method='spearman')
  
  # plot predictors and doublet scores by orientation as density plots/histograms  
  if(make_plots){
    
    # doublet counts by orientation (histogram)
    templot <- pOrient[-chrkey] %>% as.data.frame %>% `names<-`(c('Orientation')) %>%  
      mutate(PairScore = neighbor_grab(bplus,chrkey)) %>%   
      mutate(Orientation = recode_factor(Orientation, c = "Convergent", t = "Tandem", d = "Divergent")) %>%
      mutate(Orientation = ordered(Orientation, levels = c('Divergent','Tandem','Convergent'))) %>% 
      ggplot(aes(x=PairScore,fill=Orientation,colour=Orientation)) + 
      geom_histogram(size=0.5, binwidth=1, position = 'identity', bins=30,alpha=0.4) + theme_bw() + scale_colour_viridis_d() +
      ylab('Count') + theme(legend.position="bottom") + xlab('Doublet Count')
    
    ggsave(paste0(wd,'/plots/','SupplFigure03_DoubletCount.pdf'),plot = templot,device = NULL,
           path = NULL,scale = 1,width = 3.5,height = 3,units = c("in"))
    rm(templot) 
    
    # distance by orientation (density)
    templot <- pOrient[-chrkey] %>% as.data.frame %>% `names<-`(c('Orientation')) %>%  
      mutate(Distance = neighbor_grab(intrachr_dist,chrkey)) %>%   
      mutate(Orientation = recode_factor(Orientation, c = "Convergent", t = "Tandem", d = "Divergent")) %>%
      mutate(Orientation = ordered(Orientation, levels = c('Divergent','Tandem','Convergent'))) %>% 
      ggplot(aes(x=Distance,fill=Orientation)) + 
      geom_density(size=0.5, alpha=0.4) + theme_bw() + xlim(0,1e4) + scale_fill_viridis_d() +
      ylab('Density') + theme(legend.position="right")
    
    
    ggsave(paste0(wd,'/plots/','Figure04_Distance.pdf'),plot = templot,device = NULL,
           path = NULL,scale = 1,width = 7.08,height = 2,units = c("in"))
    rm(templot) 
    
    # TFBS similarity by orientation (density)
    templot <- pOrient[-chrkey] %>% as.data.frame %>% `names<-`(c('Orientation')) %>%  
      mutate(jacTF = neighbor_grab(jacTF,chrkey)) %>%   
      mutate(Orientation = recode_factor(Orientation, c = "Convergent", t = "Tandem", d = "Divergent")) %>%
      mutate(Orientation = ordered(Orientation, levels = c('Divergent','Tandem','Convergent'))) %>% 
      ggplot(aes(x=jacTF,fill=Orientation,color=Orientation)) + 
      geom_density(size=0.5, alpha=0.4) + theme_bw() + scale_fill_viridis_d() +
      ylab('Density') + theme(legend.position="bottom") + xlab('TFBS')
    
    ggsave(paste0(wd,'/plots/','SupplFigure03_TFBS.pdf'),plot = templot,device = NULL,
           path = NULL,scale = 1,width = 3.5,height = 3,units = c("in"))
    rm(templot) 
    
    
    # Chromatin baseline similarity by orientation (density)
    templot <- pOrient[-chrkey] %>% as.data.frame %>% `names<-`(c('Orientation')) %>%  
      mutate(Permissiveness = baseWeiner) %>%   
      mutate(Orientation = recode_factor(Orientation, c = "Convergent", t = "Tandem", d = "Divergent")) %>%
      mutate(Orientation = ordered(Orientation, levels = c('Divergent','Tandem','Convergent'))) %>% 
      ggplot(aes(x=Permissiveness,fill=Orientation,color=Orientation)) + 
      geom_density(size=0.5, alpha=0.4) + theme_bw() + scale_fill_viridis_d() +
      ylab('Density') + theme(legend.position="bottom") + xlab('Chromatin Baseline')
    
    ggsave(paste0(wd,'/plots/','SupplFigure03_ChromatinBaseline.pdf'),plot = templot,device = NULL,
           path = NULL,scale = 1,width = 3.5,height = 3,units = c("in"))
    rm(templot) 
    
    
    # Chromatin change similarity by orientation (density)
    templot <- pOrient[-chrkey] %>% as.data.frame %>% `names<-`(c('Orientation')) %>%  
      mutate(ChromatinSimilarity = deltaWeiner) %>%   
      mutate(Orientation = recode_factor(Orientation, c = "Convergent", t = "Tandem", d = "Divergent")) %>%
      mutate(Orientation = ordered(Orientation, levels = c('Divergent','Tandem','Convergent'))) %>% 
      ggplot(aes(x=ChromatinSimilarity,fill=Orientation,color=Orientation)) + 
      geom_density(size=0.5, alpha=0.4) + theme_bw() + scale_fill_viridis_d() +
      ylab('Density') + theme(legend.position="bottom") + xlab('Chromatin Change')
    
    ggsave(paste0(wd,'/plots/','supplFigure03_ChromatinChange.pdf'),plot = templot,device = NULL,
           path = NULL,scale = 1,width = 3.5,height = 3,units = c("in"))
    rm(templot)
    
  }
  
  rm(scoresD, scoresT, scoresC)
}

#########################################
# NEIGHBOR LINEAR MODEL TYPE III ANOVAS #
#########################################
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
  
  # provide exposure variables for both genes in each pair 
  exposure <- apply(b_mod,1,function(x)sum(x!=0))
  g1exposure <- c(NA,exposure[2:length(exposure)])[-chrkey][z]
  g2exposure <- c(NA,exposure[1:(length(exposure)-1)])[-chrkey][z]
  
  if(check_sanity){
    
    # this is an alternate way to get the exposure numbers for genes 1 and 2
    g1 <- neighbor_grab(sqgen(apply(b_mod,1,function(x)sum(x!=0))),chrkey)[z]
    g2 <- neighbor_grab(t(sqgen(apply(b_mod,1,function(x)sum(x!=0)))),chrkey)[z]
    
    # the lengths of g1,g2,g1exposure,and g2exposure should all be the same
    length(g1exposure) == length(g1)
    length(g2exposure) == length(g2)
    length(g1) == length(g2)
    
    # these ratios should be 1 if g1 exposure and g2 are the same length and each entry is equal
    sum(g1exposure == g1)/length(g1)
    sum(g2exposure == g2)/length(g2)
  }
  
  rm(exposure)
  
  # Negative binomial model with functional similarity metrics
  fxSim_model <- glm.nb(doubletScore ~ 
                          PairOrientation +
                          TFInventory +
                          Closeness + 
                          BaselineChromatin +
                          DeltaChromatin + 
                          geneInteractionSimilarity + 
                          geneOntologySimilarity + 
                          offset(log(g1exposure)) + 
                          offset(log(g2exposure)),
                        maxit=1000)
  
  fxSim_t3 <- fxSim_model %>% 
    Anova(type='III') %>%
    `names<-`(c('LR_Chisq','Df','pVals'))
  
  # format it
  fxSim_t3 <- fxSim_t3 %>% add_column(Factor=rownames(fxSim_t3)) %>%
    mutate(Significance=cut(pVals, breaks=c(-Inf,0.001,0.05,Inf), 
                            labels=c("p < 0.001","p < 0.05","p > 0.05")))
  
  write.csv(fxSim_t3[-(which(names(fxSim_t3) == 'Significance'))] %>% `rownames<-`(c(NULL)),
            file=paste0(wd,'/tables/','supplTable02_FxSimilarity.csv'))
  
  # Negative binomial model with functional similarity metrics and no TF
  fxNoTF_model <- glm.nb(doubletScore ~ 
                           PairOrientation +
                           Closeness + 
                           BaselineChromatin +
                           DeltaChromatin + 
                           geneInteractionSimilarity + 
                           geneOntologySimilarity + 
                           offset(log(g1exposure)) + 
                           offset(log(g2exposure)),
                         maxit=1000)
  
  fxNoTF_t3 <- fxNoTF_model %>% 
    Anova(type='III') %>%
    `names<-`(c('LR_Chisq','Df','pVals'))
  
  # format it
  fxNoTF_t3 <- fxNoTF_t3 %>% add_column(Factor=rownames(fxNoTF_t3)) %>%
    mutate(Significance=cut(pVals, breaks=c(-Inf,0.001,0.05,Inf), 
                            labels=c("p < 0.001","p < 0.05","p > 0.05")))
  
  write.csv(fxNoTF_t3[-(which(names(fxNoTF_t3) == 'Significance'))] %>% `rownames<-`(c(NULL)),
            file=paste0(wd,'/tables/','SupplTable02_FxSimilarityNoTF.csv'))
  
  # Negative binomial model without functional similarity metrics 
  bGlmNb_model <- glm.nb(doubletScore ~ 
                           PairOrientation +
                           TFInventory +
                           Closeness + 
                           BaselineChromatin +
                           DeltaChromatin + 
                           offset(log(g1exposure)) + 
                           offset(log(g2exposure)),
                         maxit=1000)
  
  # Type III ANOVA of the above model
  bGlmNb_t3 <- bGlmNb_model %>% 
    Anova(type='III') %>%
    `names<-`(c('LR_Chisq','Df','pVals'))
  
  # format it
  bGlmNb_t3 <- bGlmNb_t3 %>% add_column(Factor=rownames(bGlmNb_t3)) %>%
    mutate(Significance=cut(pVals, breaks=c(-Inf,0.001,0.05,Inf), 
                            labels=c("p < 0.001","p < 0.05","p > 0.05")))
  
  # FIGURE OUT HOW MUCH THIS CORRECTS FOR
  print((bGlmNb_model$null.deviance-bGlmNb_model$deviance)/bGlmNb_model$null.deviance)
  
  # HOW SIGNIFICANT IS THIS MODEL
  nullModel <- (glm.nb(doubletScore ~ 1))
  print(anova(bGlmNb_model,nullModel))
  
  
  # ANOVAS FOR EACH ORIENTATION 
  ornts <- names(table(PairOrientation))
  oModels <- rep(list(NA),length(ornts))
  names(oModels) <- ornts
  
  # Look at each orientation individually 
  
  for(i in 1:length(ornts)){
    
    pO <- PairOrientation == ornts[i]
    
    g1e <- g1exposure[pO]
    g2e <- g2exposure[pO]
    
    doubletScore <- neighbor_grab(bplus,chrkey)[z][pO]
    TFInventory <- neighbor_grab(jacTF,chrkey)[z][pO]
    Closeness <- log(1/neighbor_grab(intrachr_dist,chrkey))[z][pO]
    BaselineChromatin <- baseWeiner[z][pO]
    DeltaChromatin <- deltaWeiner[z][pO]
    
    oModels[[i]] <- glm.nb(doubletScore ~
                             TFInventory +
                             Closeness + 
                             BaselineChromatin +
                             DeltaChromatin + 
                             offset(log(g1e)) + 
                             offset(log(g2e)),
                           maxit=1000)
    
  }
  
  oT3s <- lapply(oModels,function(x) x %>%
                   Anova(type='III') %>% 
                   as.data.frame())
  
  oT3s <- lapply(oT3s,function(x)x %>% `names<-`(c('LR_Chisq','Df','pVals')) %>% 
                   add_column(Factor=rownames(oT3s[[i]])) %>%
                   mutate(Significance=cut(pVals, breaks=c(-Inf,0.001,0.05,Inf), 
                                           labels=c("p < 0.001","p < 0.05","p > 0.05"))) %>%
                   as.data.frame)
  
  # add a column explaining orientation
  oT3s$c <- oT3s$c %>% add_column(Orientation = 'Convergent') 
  oT3s$t <- oT3s$t %>% add_column(Orientation = 'Tandem') 
  oT3s$d <- oT3s$d %>% add_column(Orientation = 'Divergent')
  oT3s$a <- bGlmNb_t3 %>% add_column(Orientation = 'All')
  
  oT3s <- do.call(rbind,oT3s)
  
  if(make_plots){
    
    write.csv(oT3s[-(which(names(oT3s) == 'Significance'))] %>% `rownames<-`(c(NULL)),
              file=paste0(wd,'/tables/','SupplTable02_MainModels.csv'))
    
    #  EXTRACT THE SIGNS OF THE GLM.NB MODEL
    aCoefs <- bGlmNb_model$coefficients %>% as.data.frame() %>% `names<-`(c('Coefficient')) %>% 
      add_column(Factor=names(bGlmNb_model$coefficients)) %>% add_column(Orientation = 'All')
    
    aCoefs <- aCoefs[c(2,3,1)]
    aCoefs$Coefficient <- sapply(aCoefs$Coefficient,function(x)if(x > 0){return(1)}else {return(-1)})
    
    oCoefs <- sapply(oModels,function(x)x$coefficients) %>% as.data.frame() %>% (function(x)x[-1,])
    oCoefs <- oCoefs %>% add_column(Factor=rownames(oCoefs)) %>% 
      gather(1:(ncol(oCoefs)),key=Orientation,value=Coefficient) %>%
      mutate(Orientation = recode_factor(Orientation, c = "Convergent", t = "Tandem", d = "Divergent")) 
    oCoefs$Coefficient <- sapply(oCoefs$Coefficient,function(x)if(x > 0){return(1)}else {return(-1)})
    
    coSigns <- rbind(aCoefs,oCoefs)
    
    rm(aCoefs,oCoefs)
    
    # Apply slope signs to likelihood ratios from the type 3 ANOVAS
    for(i in 1:nrow(oT3s)){
      fo <- oT3s[i,which(names(oT3s) == 'Factor')]
      oo <- oT3s[i,which(names(oT3s) == 'Orientation')]
      
      for(j in 1:nrow(coSigns)){
        fs <- coSigns[j,which(names(coSigns) == 'Factor')]
        os <- coSigns[j,which(names(coSigns) == 'Orientation')]
        
        if((fs == fo) & (oo == os)){
          oT3s[i,which(names(oT3s) == 'LR_Chisq')] <- oT3s[i,which(names(oT3s) == 'LR_Chisq')] * coSigns[j,which(names(coSigns) == 'Coefficient')]
        }
      }
      
      
    }
    
    rm(fo,oo,fs,os)
    
    templot <- oT3s %>% mutate(Factor = recode_factor(Factor, 
                                                      PairOrientation = 'Pair Orientation',
                                                      BaselineChromatin = 'Chromatin Baseline', 
                                                      DeltaChromatin = 'Chromatin Change',
                                                      TFInventory = 'TFBS')) %>%
      mutate(Factor = ordered(Factor, levels = c(
        'Pair Orientation',
        'TFBS',
        'Chromatin Baseline',
        'Chromatin Change',
        'Closeness'))) %>%
      mutate(Orientation = ordered(Orientation, levels = c('All','Divergent','Tandem','Convergent'))) %>% 
      mutate(Significance = recode(Significance, 
                                   `p < 0.001` = "p < 0.001", 
                                   `p < 0.05` = "0.001 \u2264 p < 0.05", 
                                   `p > 0.05` = "p \u2265 0.05")) %>% 
      ggplot(aes(Factor,y=LR_Chisq,fill=Significance)) + geom_bar(stat="identity") + theme_bw() +
      facet_wrap('Orientation', scales = "free_x") + scale_fill_grey() + ylab(bquote('Signed Likehlihood Ratio'~(X^2))) +
      theme(axis.text.x.bottom = element_text(vjust = 0.9,angle = 45, hjust = 1),
            legend.position = 'right') + xlab(element_blank())
    
    ggsave(paste0(wd,'/plots/','Figure04_typeIII_results.pdf'),plot = templot,device = cairo_pdf,
           path = NULL,scale = 1,width = 7.08,height = 4,units = c("in")) 
    
  }
  
}

rm(list=ls())