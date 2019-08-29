# Adapted from Jochen's script
# Changes made by Jason include:
# 1. Changed how statistics are done, now uses resampling
# 2. Changed how data is handled at many stages
# 3. Now confidence intervals are used to specify Bottleneck/Song's filter cutoff
# 4. Now the output contains confidence intervals, rather than standard error
# 5. Added graphs to help determine the minimum nonselect counts cutoff introduced by Rachel, by region
# 6. Calculates correlation between substitutions producing same AA, by region
# Changes made by Rachel include:
# 1. bottleneck and song's filter now filter out variants that are 0 in EITHER replicate (not the mean of the replicates)
# 2. add a plot for choosing position cutoff for calculating the median stop log(phi) 
# 3. multiple plots for visualization of data
# 4. can now turn off bottleneck (Select) filter
# 5. option to change the number of SDs to use in bottleneck and Song's filter (more or less conservative)
# 6. Include option to use a stop cutoff for calculating stop median
# 7. Can customize how many variants must pass filters for calculating medians
# 8. Can customize the quality threshhold for median calculation (and an alternative one in the case not enough variants pass the first one)
# 9. Added a nonselect count filter, remove any variants that are not abundant enough in nonselect condition

#Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of tileseqMave.
#
# tileseqMave is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tileseqMave is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with tileseqMave.  If not, see <https://www.gnu.org/licenses/>.

#'
#' This analysis function performs the following steps for each mutagenesis region:
#' 1. Construction of HGVS variant descriptor strings.
#' 2. Collapsing equivalent codons into amino acic change counts.
#' 3. Error regularization at the level of pre- and post-selection counts.
#' 4. Quality-based filtering filtering based on "Song's rule". AND BOTTLENECK FILTER
#' 5. Fitness score calculation and error propagation.
#' 6. Secondary error regularization at the level of fitness scores.
#' 7. Determination of synonymous and nonsense medians and re-scaling of fitness scores.
#' 8. Flooring of negative scores and adjustment of associated error.
#' 9. Output in MaveDB format.
#' 
#' Plus several new changes:
#' 1. bottleneck and song's filter now filter out variants that are 0 in EITHER replicate (not the mean of the replicates)
#' 2. add a plot for choosing position cutoff for calculating the median stop log(phi) 
#' 3. multiple plots for visualization of data
#' 4. can now turn off bottleneck (Select) filter
#' 5. option to change the number of SDs to use in bottleneck and Song's filter (more or less conservative)
#' 6. Include option to use a stop cutoff for calculating stop median
#' 7. Can customize how many variants must pass filters for calculating medians
#' 8. Can customize the quality threshhold for median calculation (and an alternative one in the case not enough variants pass the first one)
#' 9. Added a nonselect count filter, remove any variants that are not abundant enough in nonselect condition

#' 
#' @param countfile the path to the "rawData.txt" file produced by the legacy pipeline.
#' @param regionfile the path to a csv file describing the mutagenesis regions. Must contain columns 
#'  'region', start', 'end', 'syn', 'stop', i.e. the region id, the start position, end position, and
#'  and optional synonymous and stopm mean overrides.
#' @param outdir path to desired output directory
#' @param logger a yogilogger object to be used for logging (or NULL for simple printing)
#' @param inverseAssay a boolean flag to indicate that the experiment was done with an inverse assay
#'       i.e. protein function leading to decreased fitness. Defaults to FALSE
#' @param min_nonselect_counts Vector. Nonselective counts filter cutoff (in reads/million).
#'       Filter all variants below this count out of the entire region. Default is c(0,0,0,...).
#' @param stop_cutoff Integer. The amino position above which to not include stop mutations in calculation
#'       of the stop median for scaling fitness scores. Default is NULL, in which case the length of the
#'       gene is chosen. This script outputs a plot in outdir to help choose this parameter in
#'       future runs.
#' @param sdCutoff Numeric. The standard deviation cutoff to use to choose high
#'       confidence variants to calculate stop and syn medians for fitness score scaling.
#'       Default is 0.3.
#' @param sdCutoffAlt Numeric. Should be higher than sdCutoff. The alternative standard deviation cutoff to use to choose high
#'       confidence variants to calculate stop and syn medians if not enough variants pass the
#'       initial cutoff sdCutoff.
#'       Default is 1.    
#' @param min_variants_to_choose_median Integer. The minimum number of variants that must pass
#'       pass the quality filters in order to calculate stop and syn medians. Note that
#'       if this amount is not met after the alternative filter, then no filter will be used.  
#' @param ciSong Vector of two percentages between 0 and 100.
#'       The confidence threshold to use for the nonselect filter (aka Song's filter).
#'       Variant's within this confidence level of the control mean 
#'       will be termed "likely 0" and filtered out. (only the upper bound is actually used, 
#'       the lower bound is included for clarity of function)
#' @param ciBottle Vector of two percentages between 0 and 100.
#'       The confidence threshold to use for the select filter (aka bottleneck filter).
#'       Variant's within this confidence level of the control mean 
#'       will be termed "likely 0" and filtered out. (only the upper bound is actually used, 
#'       the lower bound is included for clarity of function)
#' @param ci Vector of two percentages between 0 and 100.
#'       The confidence interval to be displayed with the output.
#' @return nothing. output is written to various files in the output directory
#' 
#' @import ggplot2
#' 
#' @export
#' 
my_analyzeLegacyTileseqCounts <- function(countfile,
                                          regionfile,
                                          outdir,
                                          logger=NULL,
                                          inverseAssay=FALSE,
                                          ns_filt_num_sd=3,
                                          select_filt_num_sd=3, 
                                          select_filt=T, 
                                          min_nonselect_counts=c(0,0,0,0,0,0,0,0,0,0), #Minimum nonselect count for each region. Length is 10 here, but should be set to a vector of length equal to number of regions
                                          stop_cutoff=NULL, 
                                          sdCutoff=3, 
                                          sdCutoffAlt=4, 
                                          min_variants_to_choose_median=10,
                                          nb=100, #Number of bootstrap samples to take. When testing, set to 10 to save time; otherwise set to 100 or more.
                                          ciSong = c(0,97.5), #Confidence interval to use for Song's filter
                                          ciBottle = c(0,97.5), #Confidence interval to use for Bottlenecking filter
                                          ci = c(2.5,97.5) #Confidence interval to display with final scores
) {
  options(warn=-1)
  #countfile <- "Pipeline/countfile.txt"
  #regionfile <- "Pipeline/regions.txt"
  #outdir <- "Test2/"
  # logger <- NULL
  # inverseAssay <- F
  # pseudoObservations <- 2
  # conservativeMode <- T
  # num_sd <- 3
  
  library(hgvsParseR)
  library(yogilog)
  library(yogitools)
  library(ggplot2)
  library(histogram)
  
  options(stringsAsFactors=FALSE)
  
  if (!is.null(logger)) {
    stopifnot(inherits(logger,"yogilogger"))
  }
  
  logInfo <- function(...) {
    if (!is.null(logger)) {
      logger$info(...)
    } else {
      do.call(cat,c(list(...),"\n"))
    }
  }
  logWarn <- function(...) {
    if (!is.null(logger)) {
      logger$warning(...)
    } else {
      do.call(cat,c("Warning:",list(...),"\n"))
    }
  }
  
  ##############
  # Read and validate input data
  ##############
  
  #countfile <- "/home/jweile/projects/ccbr2hgvs/HMGCR_S_resultfile/rawData.txt"
  #regionfile <- "/home/jweile/projects/ccbr2hgvs/HMGCR_S_resultfile/regions.txt"
  canRead <- function(filename) file.access(filename,mode=4) == 0
  stopifnot(
    canRead(countfile),
    canRead(regionfile)
  )
  
  rawCounts <- read.delim(countfile)
  regions <- read.delim(regionfile)
  
  stopifnot(
    c(
      "wt_codon","pos","mut_codon","wt_aa","mut_aa",
      "nonselect1","nonselect2","select1","select2",
      "controlNS1","controlNS2","controlS1","controlS2"
    ) %in% colnames(rawCounts),
    all(apply(regions[,1:3],2,class)=="integer"),
    all(apply(regions[,4:5],2,class) %in% c("integer","numeric","logical")),
    c("region","start","end","syn","stop") %in% colnames(regions)
  )
  #make sure outdir ends with a "/"
  if (!grepl("/$",outdir)) {
    outdir <- paste0(outdir,"/")
  }
  #and if it doesn't exist, create it
  if (!dir.exists(outdir)) {
    dir.create(outdir,recursive=TRUE)
  }
  
  
  ##################
  # Build HGVS variant descriptor strings
  ##################
  
  #for nucleotide level descriptors
  cbuilder <- new.hgvs.builder.c()
  #for protein-level descriptors
  pbuilder <- new.hgvs.builder.p(aacode=3)
  
  #for nucleotide level descriptors
  hgvsc <- apply(rawCounts,1,function(row) {
    pos <- as.numeric(row["pos"])
    #codon start indices
    cstart <- pos*3-2
    wt <- row["wt_codon"]
    mut <- row["mut_codon"]
    #calculate differences between codons
    diffs <- sapply(1:3,function(i)substr(wt,i,i)!=substr(mut,i,i))
    ndiff <- sum(diffs)
    if (ndiff == 1) { #one difference is a SNV
      offset <- which(diffs)
      wtbase <- substr(wt,offset,offset)
      mutbase <- substr(mut,offset,offset)
      snvpos <- cstart+offset-1
      return(cbuilder$substitution(snvpos,wtbase,mutbase))
    } else if (ndiff > 1) { #multiple differences is a delIns
      return(cbuilder$delins(cstart,cstart+2,mut))
    } else {
      stop("mutation must differ from wt!")
    }
  })
  
  hgvsp <- apply(rawCounts,1,function(row) {
    pos <- as.numeric(row["pos"])
    wt <- row["wt_aa"]
    mut <- row["mut_aa"]
    if (mut=="_") mut <- "*" #correct stop character
    if (wt == mut) {
      return(pbuilder$synonymous(pos,wt))
    } else {
      return(pbuilder$substitution(pos,wt,mut))
    }
  })
  
  logInfo(sprintf(
    "Parsed data for %d variants covering %d amino acid changes",
    length(hgvsc),length(unique(hgvsp))
  ))
  
  #####################
  # Pre-processing and formatting input
  ####################
  
  #extract and examine condition names and replicates
  conditions <- colnames(rawCounts)[-c(1:6)]
  condStruc <- extract.groups(conditions,"^(\\w+)(\\d+)$")
  condNames <- unique(condStruc[,1])
  repNames <- unique(condStruc[,2])
  condMatrix <- do.call(rbind,lapply(condNames,paste0,repNames))
  dimnames(condMatrix) <- list(condNames,repNames)
  repMatrix <- head(condMatrix,length(condNames)/2)
  #TODO: Throw error if conditions are missing
  
  logInfo(sprintf(
    "Detected %d replicates for %d conditions: %s",
    length(repNames), length(unique(condNames)), paste(condNames,collapse=", ")
  ))
  
  #Mark Regions
  Region <- matrix(0, nrow(rawCounts), 1)
  colnames(Region) <- "Region"
  rawCounts <- cbind(rawCounts,Region)
  for (i in 1:nrow(regions)) {
    for (j in 1:nrow(rawCounts))
      if (regions[i,2]<= rawCounts[j,"pos"] & rawCounts[j,"pos"]<= regions[i,3])
        rawCounts[j,"Region"] <- i
  }
  
  # =====================================================
  # apply pre-filter on variants
  # =====================================================

  #Initialize vectors containing the result of each bootstrap sample
  theta_Song <- replicate(nb,0)
  theta_Bottle <- replicate(nb,0)
  
  #Get mean and maximum likelihood Standard Deviation for bootstrapping
  mcv <- do.call(cbind,lapply(condNames,function(cond) {
    mcv <- t(apply(rawCounts[,condMatrix[cond,]],1,function(xs){
      m <- mean(xs)
      sdev <- sd(xs,na.rm=TRUE)
      sdevML <- sdev*sqrt((length(xs)-1)/length(xs))
      c(m,sdevML)
    }))
    colnames(mcv) <- paste0(cond,c(".mean",".sdevML"))
    mcv
  }))
  mcv <- cbind(mcv,matrix(0, dim(mcv)[1], 4))
  for (i in c(1:dim(mcv)[1])) {
    if (mcv[i,5]==0) {theta_Song <- replicate(nb,0)}
    else {
      for (j in 1:nb){
        BootSamples <- c(rnorm(length(repNames), mcv[i,5], mcv[i,6]))
        theta_Song[j] <- mean(BootSamples[1],BootSamples[2])
    }}
    if (mcv[i,7]==0) {theta_Bottle <- replicate(nb,0)}
    else {
      for (j in 1:nb){
        BootSamples <- c(rnorm(length(repNames), mcv[i,7], mcv[i,8]))
        theta_Bottle[j] <- mean(BootSamples[1],BootSamples[2])
      }}
    ConfSong <- quantile(theta_Song,c(ciSong[1]/100,ciSong[2]/100),na.rm=TRUE)
    ConfBottle <- quantile(theta_Bottle,c(ciBottle[1]/100,ciBottle[2]/100),na.rm=TRUE)
    mcv[i,9] <- ConfSong[1]
    mcv[i,10] <- ConfSong[2]
    mcv[i,11] <- ConfBottle[1]
    mcv[i,12] <- ConfBottle[2]
  }
  colnames(mcv)[9] <- "NonSelect.CILower"
  colnames(mcv)[10] <- "NonSelect.CIUpper"
  colnames(mcv)[11] <- "Select.CILower"
  colnames(mcv)[12] <- "Select.CIUpper"
  
  rawCounts <- cbind(rawCounts,mcv)
  
  # Sequencing error filter (song's filter) / Nonselect filter
  flagged1 <- with(as.data.frame(rawCounts), {
    (NonSelect.CIUpper >= nonselect1) | (NonSelect.CIUpper >= nonselect2)
  })
  # Bottleneck filter / Select filter
  flagged2 <- with(as.data.frame(rawCounts), {
    (Select.CIUpper >= select1) | (Select.CIUpper >= select2)
  })
  if (select_filt==F) { # turn off the bottleneck filter if select_filt is false
    flagged2 <- rep.int(FALSE, times = nrow(rawCounts))
  }
  #Minimum nonselect counts filter
  flagged3=replicate(nrow(rawCounts),0)
  counter <- 0
  for (j in 1:nrow(rawCounts)){
    reg <- as.numeric(rawCounts[j,"Region"])
    flagged3[j] <- rawCounts[j,"nonselect.mean"] < min_nonselect_counts[reg]}
  #Maximum nonselect counts filter (use sparingly)
  flagged4 <- with(as.data.frame(rawCounts), {
    nonselect.mean > quantile(rawCounts[,"nonselect.mean"],0.99)
  })
  logInfo(sprintf(
    "Filtering out %d variants (=%.02f%%):
%d (=%.02f%%) due to likely sequencing error.
%d (=%.02f%%) due to likely bottlenecking.
%d (=%.02f%%) due to insufficient representation in nonselect condition.
%d (=%.02f%%) due to over-representation in nonselect condition.",
    sum(flagged1|flagged2|flagged3|flagged4), 100*sum(flagged1|flagged2|flagged3|flagged4)/length(flagged1|flagged2|flagged3|flagged4),
    sum(flagged1), 100*sum(flagged1)/length(flagged1),
    sum(flagged2), 100*sum(flagged2)/length(flagged2),
    sum(flagged3), 100*sum(flagged3)/length(flagged3),
    sum(flagged4), 100*sum(flagged4)/length(flagged4)
  ))
  rawCounts <- rawCounts[!(flagged1|flagged2|flagged3|flagged4),]
  hgvsc <- hgvsc[!(flagged1|flagged2|flagged3|flagged4)]
  hgvsp <- hgvsp[!(flagged1|flagged2|flagged3|flagged4)]
  logInfo(sprintf(
    "Data remains for for %d variants covering %d amino acid changes",
    length(hgvsc),length(unique(hgvsp))
  ))
  
  
  #collapse codons into unique AA changes
  logInfo("Collapsing variants by outcome...")
  rawCounts <- as.data.frame(cbind(
    hgvsp,
    hgvsc,
    rawCounts
  ))
  
  
  #########
  # Compute bootstrap sample matrix
  #########
  
  #Reset mcv
  mcv <- do.call(cbind,lapply(condNames,function(cond) {
    mcv <- t(apply(rawCounts[,condMatrix[cond,]],1,function(xs){
      m <- mean(xs)
      sdev <- sd(xs,na.rm=TRUE)
      sdevML <- sdev*sqrt((length(xs)-1)/length(xs))
      c(m,sdevML)
    }))
    colnames(mcv) <- paste0(cond,c(".mean",".sdevML"))
    mcv
  }))

  #Initialize Bootstrap Matrix
  BootMatNuc <- matrix(0,nb*nrow(rawCounts),14)
  CrudeScores <- replicate(0,2)
  #Generate Bootstrap simulations
  colnames(BootMatNuc) <- c("nonselect1","nonselect2","select1","select2",
                         "controlNS1","controlNS2","controlS1","controlS2",
                         "CrudeScore1","CrudeScore2","MeanScore",
                         "hgvsp","hgvsc","Region")
  for (i in 1:nb) {
    for (j in ((i-1)*nrow(rawCounts)+1):((i-1)*nrow(rawCounts)+nrow(rawCounts))){
      k <- j-(i-1)*nrow(rawCounts)
      nonselect.b <- c(rnorm(2, mcv[k,1], mcv[k,2]))
      select.b <- rnorm(2, mcv[k,3], mcv[k,4])
      nonselectC.b <- rnorm(2, mcv[k,5], mcv[k,6])
      selectC.b <- rnorm(2, mcv[k,7], mcv[k,8])
      BootMatNuc[j,1] <- nonselect.b[1]
      BootMatNuc[j,2] <- nonselect.b[2]
      BootMatNuc[j,3] <- select.b[1]
      BootMatNuc[j,4] <- select.b[2]
      BootMatNuc[j,5] <- nonselectC.b[1]
      BootMatNuc[j,6] <- nonselectC.b[2]
      BootMatNuc[j,7] <- selectC.b[1]
      BootMatNuc[j,8] <- selectC.b[2]
      BootMatNuc[j,9] <- (select.b[1]-selectC.b[1])/(nonselect.b[1]-nonselectC.b[1])
      BootMatNuc[j,10] <- (select.b[2]-selectC.b[2])/(nonselect.b[2]-nonselectC.b[2])
      BootMatNuc[j,11] <- exp(mean(c(log(as.numeric(BootMatNuc[j,9])),log(as.numeric(BootMatNuc[j,10])))))
      BootMatNuc[j,12] <- rawCounts[k,"hgvsp"]
      BootMatNuc[j,13] <- rawCounts[k,"hgvsc"]
      BootMatNuc[j,14] <- rawCounts[k,"Region"]}}
  #Calculate Combined AA Scores
  BootMatAA <- BootMatNuc[,10:13]
  hgvsp <- unique(rawCounts[,"hgvsp"])
  counter <- 0
  for (amino in hgvsp){
    aminoRows <- which(rawCounts$hgvsp == amino)
    aminoRowsBoot <- numeric(nb*length(aminoRows))
    for (i in 1:length(aminoRows)){
      for (j in 1:nb){
        aminoRowsBoot[j+(i-1)*nb] <- aminoRows[i]+(j-1)*nrow(rawCounts)}}
    BootLocal <- BootMatNuc[aminoRowsBoot,]
    CrudeScores <- numeric(2*dim(BootLocal)[1]/nb)
    for (i in 1:nb){
      for (j in 1:(dim(BootLocal)[1]/nb)) {
        CrudeScores[j] <- as.numeric(BootLocal[i+(j-1)*nb,"CrudeScore1"])
        CrudeScores[j+dim(BootLocal)[1]/nb] <- as.numeric(BootLocal[i+(j-1)*nb,"CrudeScore2"])
      }
      BootMatAA[i+counter,1] <- amino
      BootMatAA[i+counter,2] <- exp(mean(log(CrudeScores)))
      BootMatAA[i+counter,3] <- sd(CrudeScores)
      BootMatAA[i+counter,4] <- BootLocal[1,"Region"]
    }
    counter <- counter+nb
  }
  BootMatAA <- BootMatAA[1:counter,]
  colnames(BootMatAA) <- c("hgvsp","MeanScore","SDev","Region")
  
  #Turn our bootstrap matrices into data frames
  BootMatAA <- data.frame(BootMatAA)
  BootMatNuc <- data.frame(BootMatNuc[,c("hgvsc","MeanScore","CrudeScore2","Region","hgvsp","nonselect1")]) 
  #^ CrudeScore2 is an arbitrary choice for placeholder data, crude scores are not actually used past this point
  
  #Change some data to numeric
  BootMatAA[,"MeanScore"] <- as.numeric(as.character(BootMatAA[,"MeanScore"]))
  BootMatAA[,"SDev"] <- as.numeric(as.character(BootMatAA[,"SDev"]))
  BootMatAA[,"Region"] <- as.numeric(as.character(BootMatAA[,"Region"]))
  BootMatNuc[,"MeanScore"] <- as.numeric(as.character(BootMatNuc[,"MeanScore"]))
  BootMatNuc[,"CrudeScore2"] <- as.numeric(as.character(BootMatNuc[,"CrudeScore2"]))
  BootMatNuc[,"Region"] <- as.numeric(as.character(BootMatNuc[,"Region"]))
  BootMatNuc[,"nonselect1"] <- as.numeric(as.character(BootMatNuc[,"nonselect1"]))

  OutBootMatAA <- BootMatAA[,c("hgvsp","MeanScore","Region")]
  OutBootMatNuc <- BootMatNuc[,c("hgvsc","MeanScore","Region")]
  #################
  # Estimate Modes of synonymous and stop
  #################
  
  #################
  #TODO: Try gaussian mixture models with two underlying distributions?
  #################
  
  # funciton to get the aa position from hgvsp
  get_pos_from_hgvsp <- function(variant) {
    variant <- unlist(strsplit(variant, split = ""))
    nums <- grep(pattern = "[0-9]", x = variant, value = TRUE)
    pos <- as.numeric(paste(nums, collapse = ""))
    return(pos)}
  for (j in 1:nrow(regions)) {
    SampleBootMatAA <- BootMatAA[BootMatAA[,"Region"]==j,]
    SampleBootMatNuc <- BootMatNuc[BootMatNuc[,"Region"]==j,]
    region.syn <- regions[j,"syn"]
    region.stop <- regions[j,"stop"]
    for (i in 1:nb){
      LocalRows1 <- seq(i,nrow(SampleBootMatAA),nb)
      LocalBootMatAA <- SampleBootMatAA[LocalRows1,]
      LocalRows2 <- (1+(i-1)*nrow(SampleBootMatNuc)/nb):(i*nrow(SampleBootMatNuc)/nb)
      LocalBootMatNuc <- SampleBootMatNuc[LocalRows2,]
      # GET SD CUTOFF FOR MODES ESTIMATE
      #if we can't find enough syn/stop below the cutoff, increase the cutoff to 1
      if (with(LocalBootMatAA,sum(grepl("Ter$",hgvsp)&SDev<sdCutoff)<min_variants_to_choose_median ||
               sum(grepl("=$",hgvsp)&SDev<sdCutoff)<min_variants_to_choose_median )) {
        sdCutoff <- sdCutoffAlt}
      #if we still can't find enough below the new cutoff, get rid of it altogether
      if (with(LocalBootMatAA,sum(grepl("Ter$",hgvsp)&SDev<sdCutoff)<min_variants_to_choose_median ||
               sum(grepl("=$",hgvsp)&SDev<sdCutoff)<min_variants_to_choose_median )) {
        sdCutoff <- Inf}
      # GET STOP CUTOFF FOR MODES ESTIMATE
      # get the positions 
      LocalBootMatAA$positions <- unlist(lapply(LocalBootMatAA$hgvsp, get_pos_from_hgvsp))
      # length of the gene
      length <- max(LocalBootMatAA$positions)
  
      # set the stop cutoff to the length of the gene if not provided
      if (is.null(stop_cutoff)) {
        stop_cutoff <- length}
      
      
      #Plot for manually choosing stop cutoff (only does this for one iteration to save time)
      if(1==i&j==nrow(regions)){
        tiffFile <- paste0(outdir,"region",j,"_choose_stop_cutoff.tiff")
        tiff(tiffFile, units="in", width=4, height=3, res=300)
        stop_cut_plot_2 <- 
          ggplot(data = LocalBootMatAA[grepl("Ter$",LocalBootMatAA$hgvsp),],mapping = aes(x = positions, y = MeanScore)) +
          geom_point() +
          theme_linedraw() +
          geom_smooth() +
          ylab(expression(phi)) +
          xlab('AA position')
          print(stop_cut_plot_2)
          dev.off()}
      # now we have chosen the cutoff, calculate the medians to use to scale the fitness scores
      if (j==nrow(regions)){stops <- LocalBootMatAA$MeanScore[which(grepl("Ter$",LocalBootMatAA$hgvsp) 
                                         & LocalBootMatAA$SDev < sdCutoff 
                                         & LocalBootMatAA$positions <= stop_cutoff )]}
      else {stops <- LocalBootMatAA$MeanScore[which(grepl("Ter$",LocalBootMatAA$hgvsp) 
                                                    & LocalBootMatAA$SDev < sdCutoff)]}
      syns <- LocalBootMatAA$MeanScore[which(grepl("=$",LocalBootMatAA$hgvsp) & LocalBootMatAA$SDev < sdCutoff)]
      modes <- c(stop=median(stops,na.rm=TRUE),syn=median(syns,na.rm=TRUE)) # keep this for jochen's plots but i will use a different variable name for mine
      median_syn <- round(median(syns,na.rm=TRUE), digits = 2)
      median_stop <- round(median(stops,na.rm=TRUE), digits = 2)
  
      #################
      # Use syn/stop medians to scale scores 
      #################
      
      #Only log the first iteration to save time
      if (i==1) {logInfo(sprintf(
        "Scaling to synonymous (log(phi)=%.02f) and nonsense (log(phi)=%.02f) medians.",log10(modes[["syn"]]),log10(modes[["stop"]])))
        #if manual overrides for the synonymous and stop modes were provided, use them
        if (!is.na(region.syn)) {
          logInfo(sprintf(
            "Using manual override (=%.02f) for synonmous mode instead of automatically determined value (=%.02f).",region.syn,log10(modes[["syn"]])))
            modes[["syn"]] <- 10^(region.syn)}
        if (!is.na(region.stop)) {
          logInfo(sprintf(
            "Using manual override (=%.02f) for stop mode instead of automatically determined value (=%.02f).",region.stop,log10(modes[["stop"]])))
            modes[["stop"]] <- 10^(region.stop)}}
      #apply the scaling
      denom <- log10(modes[["syn"]])-log10(modes[["stop"]])
      for (k in 1:nrow(LocalBootMatAA)){
        LocalBootMatAA[k,"MeanScore"] <- (log10(max(0.0001,LocalBootMatAA[k,"MeanScore"]))-log10(modes[["stop"]]))/denom}
      for (k in 1:nrow(LocalBootMatNuc)){
        LocalBootMatNuc[k,"MeanScore"] <- (log10(max(0.0001,LocalBootMatNuc[k,"MeanScore"]))-log10(modes[["stop"]]))/denom}
      SampleBootMatAA[LocalRows1,] <- LocalBootMatAA
      SampleBootMatNuc[LocalRows2,] <- LocalBootMatNuc}
    OutBootMatAA[OutBootMatAA[,"Region"]==j,1:2] <- SampleBootMatAA[,1:2]
    OutBootMatNuc[OutBootMatNuc[,"Region"]==j,1:2] <- SampleBootMatNuc[,1:2]
    #Caclulate bootstrap averages
    for (i in 1:(nrow(SampleBootMatAA)/nb)){
      combined <- SampleBootMatAA[(nb*(i-1)+1):(nb*i),2]
      confint <- quantile(combined,ci/100,na.rm=TRUE)
      SampleBootMatAA[i,1] <- SampleBootMatAA[(nb*(i-1)+1),1]
      SampleBootMatAA[i,2] <- mean(combined,na.rm=TRUE)
      SampleBootMatAA[i,3] <- confint[1]
      SampleBootMatAA[i,4] <- confint[2]}
    SampleBootMatAA <- SampleBootMatAA[1:(nrow(SampleBootMatAA)/nb),]
    colnames(SampleBootMatAA) <- c("hgvsp","score","CILower","CIUpper")
    for (i in 1:(nrow(SampleBootMatNuc)/nb)){
      combined <- SampleBootMatNuc[seq(i,nrow(SampleBootMatNuc),(nrow(SampleBootMatNuc)/nb)),2]
      confint <- quantile(combined,ci/100,na.rm=TRUE)
      SampleBootMatNuc[i,2] <- mean(combined,na.rm=TRUE)
      SampleBootMatNuc[i,3] <- confint[1]
      SampleBootMatNuc[i,4] <- confint[2]}
    SampleBootMatNuc <- SampleBootMatNuc[1:(nrow(SampleBootMatNuc)/nb),]
    colnames(SampleBootMatNuc) <- c("hgvsc","score","CILower","CIUpper","hgvsp","nonselect1")
    TestMat <- SampleBootMatNuc[,c("score","CILower","hgvsp","nonselect1")]
    PreTest <- TestMat
    localhgvsp <- unique(PreTest[,"hgvsp"])
    counter <- 1
    for (amino in localhgvsp){
      localbooter <- PreTest[PreTest[,"hgvsp"]==amino,]
      if (nrow(localbooter)>=2){
        TestMat[counter,1] <- (localbooter[1,"score"]-localbooter[2,"score"])^2
        TestMat[counter,2] <- min(c(localbooter[1,"nonselect1"],localbooter[2,"nonselect1"]))
        TestMat[counter,3] <- amino
        TestMat[counter,4] <- localbooter[1,"score"]
        TestMat[counter,5] <- localbooter[2,"score"]
        counter <- counter + 1}}
    TestMat <- TestMat[1:(counter-1),]
    logInfo(sprintf(
      "Correlation between substitutions producing same AA is: %f",
      cor(TestMat[,4:5])[1,2]
    ))
    tiffFile <- paste0(outdir,"region",j,"_choose_nonselect_cutoff.tiff")
    tiff(tiffFile, units="in", width=4, height=3, res=300)
    stop_cut_plot_2 <- 
      ggplot(data = TestMat,mapping = aes(x = CILower, y = score)) +
      geom_point() +
      theme_linedraw() +
      geom_smooth() +
      ylab('error') +
      xlab('nonselect count')
    print(stop_cut_plot_2)
    dev.off()
    
    #################
    # Write output to file
    #################
    
    #protein-level MaveDB output
    outfile <- paste0(outdir,"mavedb_scores_perAA_region",j,".csv")
    write.csv(SampleBootMatAA,outfile,row.names=FALSE)
    
    #nucleotide-level MaveDB output
    outfile <- paste0(outdir,"mavedb_scores_perNt_region",j,".csv")
    write.csv(SampleBootMatNuc,outfile,row.names=FALSE)
    
    # get the positions 
    SampleBootMatAA$positions <- unlist(lapply(SampleBootMatAA$hgvsp, get_pos_from_hgvsp))
    length <- max(SampleBootMatAA$positions)
    
    #Find Stops and Syns
    stops <- SampleBootMatAA$score[which(grepl("Ter$",SampleBootMatAA$hgvsp))]
    syns <- SampleBootMatAA$score[which(grepl("=$",SampleBootMatAA$hgvsp))]
    modes <- c(stop=median(stops,na.rm=TRUE),syn=median(syns,na.rm=TRUE)) # keep this for jochen's plots but i will use a different variable name for mine
    
    median_syn <- round(median(syns,na.rm=TRUE), digits = 2)
    median_stop <- round(median(stops,na.rm=TRUE), digits = 2)
    
    # plot only the syns and stops used in the median calculation
    tiffFile <- paste0(outdir,"region",j,"_syn_stop_for_scaling.tiff")
    tiff(tiffFile, units="in", width=4, height=2.5, res=300)
    hist1 <- ggplot() + 
      geom_histogram(data = data.frame(score=stops), 
                     mapping = aes(x = score, y = ..density..),  
                     fill = 'darkred', alpha = 0.5, col = 'black') + 
      geom_histogram(data = data.frame(score=syns), 
                     mapping = aes(x = score, y = ..density..),  
                     fill = 'darkgreen', alpha = 0.5, col = 'black') + 
      theme_linedraw() + 
      geom_vline(xintercept = median_stop, size=1, linetype=2, color='darkred') + 
      geom_vline(xintercept = median_syn, size=1, linetype=2, color='darkgreen') + 
      geom_text(aes(x=median_syn + 0.15, label=median_syn, y=6)) + 
      geom_text(aes(x=median_stop + 0.15, label=median_stop, y=6)) + 
      xlab(expression(log(phi)))
    print(hist1)
    dev.off()
    
    # plot ALL of the syns and stops and their medians
    stop.is <- which(grepl("Ter$",SampleBootMatAA$hgvsp))
    syn.is <- which(grepl("=$",SampleBootMatAA$hgvsp))
    miss.is <- setdiff(1:nrow(SampleBootMatAA),c(stop.is,syn.is))
    all_stop_median <- round(median(SampleBootMatAA$score[stop.is]), digits = 2)
    all_syn_median <- round(median(SampleBootMatAA$score[syn.is]), digits = 2)
    tiffFile <- paste0(outdir,"region",j,"_all_syn_stop.tiff")
    tiff(tiffFile, units="in", width=4, height=2.5, res=300)
    hist2 <- ggplot() + 
      geom_histogram(data = SampleBootMatAA[stop.is,], mapping = aes(x = score, y = ..density..),  
                     fill = 'darkred', alpha = 0.5, col = 'black') + 
      geom_histogram(data = SampleBootMatAA[syn.is,], mapping = aes(x = score, y = ..density..),  fill = 'darkgreen', alpha = 0.5, col = 'black') + 
      theme_linedraw() + 
      geom_vline(xintercept = median(SampleBootMatAA$score[stop.is]), size=1, linetype=2, color='darkred') +
      geom_vline(xintercept = median(SampleBootMatAA$score[syn.is]), size=1, linetype=2, color='darkgreen') +
      geom_text(aes(x=all_syn_median + 0.15, label=all_syn_median, y=6)) +
      geom_text(aes(x=all_stop_median + 0.15, label=all_stop_median, y=6)) +
      xlab(expression(log(phi)))
    print(hist2)
    dev.off()
    # plot all syns and stops but with the medians used for scaling
    tiffFile <- paste0(outdir,"region",j,"_all_syn_stop_with_filtered_medians.tiff")
    tiff(tiffFile, units="in", width=4, height=2.5, res=300)
    hist4 <- ggplot() +
      geom_histogram(data = SampleBootMatAA[stop.is,], mapping = aes(x = score, y = ..density..),  fill = 'darkred', alpha = 0.5, col = 'black') +
      geom_histogram(data = SampleBootMatAA[syn.is,], mapping = aes(x = score, y = ..density..),  fill = 'darkgreen', alpha = 0.5, col = 'black') +
      theme_linedraw() +
      geom_vline(xintercept = median_stop, size=1, linetype=2, color='darkred') +
      geom_vline(xintercept = median_syn, size=1, linetype=2, color='darkgreen') +
      geom_text(aes(x=median_syn + 0.15, label=median_syn, y=6)) +
      geom_text(aes(x=median_stop + 0.15, label=median_stop, y=6)) +
      xlab(expression(log(phi)))
    print(hist4)
    dev.off()
    
    # plot the missense and medians used for scaling
    tiffFile <- paste0(outdir,"region",j,"_all_missense.tiff")
    tiff(tiffFile, units="in", width=4, height=2.5, res=300)
    hist3 <- ggplot() +
      geom_histogram(data = SampleBootMatAA[miss.is,], mapping = aes(x = score),  fill = 'grey', alpha = 0.5, col = 'black') +
      theme_linedraw() +
      geom_vline(xintercept = median_stop, size=1, linetype=2, color='darkred') +
      geom_vline(xintercept = median_syn, size=1, linetype=2, color='darkgreen') +
      geom_text(aes(x=median_syn + 0.2, label=median_syn, y=250)) +
      geom_text(aes(x=median_stop + 0.2, label=median_stop, y=250)) +
      xlab(expression(log(phi)))
    print(hist3)
    dev.off()
    
    #################
    #Plot Syn vs stop
    #################
    
    plotStopSyn <- function(data,title,modes) {
      with(data,{stop.is <- which(grepl("Ter$",hgvsp))
      syn.is <- which(grepl("=$",hgvsp))
      miss.is <- setdiff(1:nrow(data),c(stop.is,syn.is))
      br <- seq(floor(min(score)),ceiling(max(score)),.1)
      hist(score[miss.is],col="gray",breaks=br,xlab=expression(log(phi)),border=NA,main=title)
      hist(score[stop.is],add=TRUE,col=colAlpha("firebrick3",.5),breaks=br)
      hist(score[syn.is],add=TRUE,col=colAlpha("darkolivegreen3",.5),breaks=br)
      abline(v=modes,col=c("firebrick3","darkolivegreen3"),lty="dashed")})}
    pdfFile <- paste0(outdir,"scalingQC_region",j,".pdf")
    pdf(pdfFile)
    op <- par(mfrow=c(2,1))
    plotStopSyn(SampleBootMatAA,"Unfiltered",modes)
    legend("topright",c("missense","synonymous","stop"),fill=c("gray","darkolivegreen3","firebrick3"))
    par(op)
    invisible(dev.off())
    if(j==1){
      OutputAA <- SampleBootMatAA[,1:4]
      OutputNuc <- SampleBootMatNuc[,1:6]}
    else {
      OutputAA <- rbind(OutputAA,SampleBootMatAA[,1:4])
      OutputNuc <- rbind(OutputNuc,SampleBootMatNuc[,1:6])}}
  colnames(OutputAA) <- c("hgvsp","score","CILower","CIUpper")
  colnames(OutputNuc) <- c("hgvsc","score","CILower","CIUpper","hgvsp","nonselect1")

  
  #################
  # Write output to file
  #################
  
  
  #protein-level MaveDB output
  outfile <- paste0(outdir,"mavedb_scores_perAA",".csv")
  write.csv(OutputAA[,1:4],outfile,row.names=FALSE)
  
  #protein-level MaveDB Bootstrap matrix
  outfile <- paste0(outdir,"mavedb_scores_perAA_boot",".csv")
  write.csv(OutBootMatAA,outfile,row.names=FALSE)
  
  #nucleotide-level MaveDB output
  outfile <- paste0(outdir,"mavedb_scores_perNt",".csv")
  write.csv(OutputNuc[,1:4],outfile,row.names=FALSE)}

#Uncomment below for testing

my_analyzeLegacyTileseqCounts(countfile,
                              regionfile,
                              outdir,
                              logger=NULL,
                              inverseAssay=FALSE,
                              min_nonselect_counts=c(600,700,700,400,500,400),
                              stop_cutoff=100000, 
                              sdCutoff=3, 
                              sdCutoffAlt=4,
                              min_variants_to_choose_median=10,
                              nb=30)