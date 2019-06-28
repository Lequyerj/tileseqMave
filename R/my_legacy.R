
# Adapted from Jochen's script
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
#' @param pseudoObservations The number of pseudoObservations to use for the Baldi&Long regularization.
#'       Defaults to 2.
#' @param conservativeMode Boolean flag. When turned on, pseudoObservations are not counted towards 
#'       standard error and the first round of regularization uses pessimistic error estimates.
#' @param ns_filt_num_sd Integer (between -3 and 3 suggested, other values unlikely to be meaningful).
#'       The number of standard deviations to be used in the nonselect filter (aka Song's filter, aka sequencing error filter).
#'       Variant's counts less than this number of standard deviations from the control mean 
#'       will be termed "likely 0" and filtered out. 3 filters out the most data (anything that could be 0),
#'       and -3 filters out the least data (anything that is definitely 0). Default is 3.
#' @param select_filt_num_sd Integer (between -3 and 3 suggested, other values unlikely to be meaningful).
#'       The number of standard deviations to be used in the select filter (aka bottleneck filter, aka backwards ORF filter).
#'       Variant's counts less than this number of standard deviations from the control mean 
#'       will be termed "likely 0" and filtered out. 3 filters out the most data (anything that could be 0),
#'       and -3 filters out the least data (anything that is definitely 0). Default is 3. This parameter
#'       is ignored if select_filt is FALSE.
#' @param select_filt Boolean flag. Indicates whether to use the select filter (aka bottleneck filter).
#'       Default is TRUE.
#' @param min_nonselect_counts Numeric. Nonselective counts filter cutoff (in reads/million).
#'       Filter all variants below this count out of the entire data set. Default is -Inf.
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
#'
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
                                          pseudoObservations=2,
                                          conservativeMode=TRUE, 
                                          ns_filt_num_sd=3,
                                          select_filt_num_sd=3, 
                                          select_filt=T, 
                                          min_nonselect_counts=-Inf,
                                          stop_cutoff=NULL, 
                                          sdCutoff=0.3, 
                                          sdCutoffAlt=1, 
                                          min_variants_to_choose_median=10
                                          ) {
  
  options(warn=-1)
  # countfile <- "../R_DMS/data/rawData_GDI1_2016Q20_MAVEtileseq_input.txt"
  # regionfile <- "../R_DMS/data/GDI1_regionfile.txt"
  # outdir <- "../R_DMS/data/"
  # logger <- NULL
  # inverseAssay <- F
  # pseudoObservations <- 2
  # conservativeMode <- T
  # num_sd <- 3

	#library(hgvsParseR)
	# library(yogilog)
	#library(yogitools)

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

	#TODO: Throw error if conditions are missing

	logInfo(sprintf(
		"Detected %d replicates for %d conditions: %s",
		length(repNames), length(unique(condNames)), paste(condNames,collapse=", ")
	))
	
	
	# =====================================================
	# apply pre-filter on variants
	# =====================================================
	rawmsd <- do.call(cbind,lapply(condNames,function(cond) {
		msd <- t(apply(rawCounts[,condMatrix[cond,]],1,function(xs){
			c(mean(xs,na.rm=TRUE), sd(xs,na.rm=TRUE))
		}))
		colnames(msd) <- paste0(cond,c(".mean",".sd"))
		msd
	}))
	# Sequencing error filter (song's filter) / Nonselect filter
	flagged1 <- with(as.data.frame(rawmsd), {
	  (controlNS.mean + ns_filt_num_sd*controlNS.sd >= rawCounts$nonselect1) | (controlNS.mean + ns_filt_num_sd*controlNS.sd >= rawCounts$nonselect2)
	})
	# Bottleneck filter / Select filter
	flagged2 <- with(as.data.frame(rawmsd), {
		(controlS.mean + select_filt_num_sd*controlS.sd >= rawCounts$select1) | (controlS.mean + select_filt_num_sd*controlS.sd >= rawCounts$select2)
	})
	if (select_filt==F) { # turn off the bottleneck filter if select_filt is false
	  flagged2 <- rep.int(FALSE, times = nrow(rawmsd))
	}
	# Nonselect counts filter
	flagged3 <- with(as.data.frame(rawmsd), {
	  nonselect.mean < min_nonselect_counts
	})
	logInfo(sprintf(
"Filtering out %d variants (=%.02f%%):
%d (=%.02f%%) due to likely sequencing error.
%d (=%.02f%%) due to likely bottlenecking.
%d (=%.02f%%) due to insufficient representation in nonselect condition.",
		sum(flagged1|flagged2|flagged3), 100*sum(flagged1|flagged2|flagged3)/length(flagged1|flagged2|flagged3),
		sum(flagged1), 100*sum(flagged1)/length(flagged1),
		sum(flagged2), 100*sum(flagged2)/length(flagged2),
    sum(flagged3), 100*sum(flagged3)/length(flagged3)
	))
	rawCountsFiltered <- rawCounts[!(flagged1|flagged2|flagged3),]
	hgvsc <- hgvsc[!(flagged1|flagged2|flagged3)]
	hgvsp <- hgvsp[!(flagged1|flagged2|flagged3)]

	logInfo(sprintf(
		"Data remains for for %d variants covering %d amino acid changes",
		length(hgvsc),length(unique(hgvsp))
	))


	#collapse codons into unique AA changes
	logInfo("Collapsing variants by outcome...")
	combiCounts <- as.df(tapply(1:nrow(rawCountsFiltered),hgvsp,function(is) {
		# mut <- unique(hgvsp[is])
		hgvsp <- hgvsp[is[[1]]]
		hgvsc <- paste(unique(hgvsc[is]),collapse=" ")
		c(list(hgvsp=hgvsp,hgvsc=hgvsc),colSums(rawCountsFiltered[is,conditions]))
	},simplify=FALSE))

	logInfo("Parsing variant strings...")
	combiCountMuts <- parseHGVS(combiCounts$hgvsp)
	combiCountMuts$type[which(combiCountMuts$variant=="Ter")] <- "nonsense"
	logInfo("Parsing complete.")

	##############
	# Iterate over regions
	##############
	for (region.i in 1:nrow(regions)) {

		regionStart <- regions[region.i,"start"]
		regionEnd <- regions[region.i,"end"]
		region.rows <- with(combiCountMuts,which(start >= regionStart & start <= regionEnd))
		region.syn <- regions[region.i,"syn"]
		region.stop <- regions[region.i,"stop"]

		if (length(region.rows) < 1) {
			logWarn(sprintf("\n\nNo variants found for region %d! Skipping...",region.i))
			next
		} else {
			logInfo(sprintf("\n\nProcessing region %d. (pos. %d-%d)",region.i,regionStart,regionEnd))
		}

		localCombiCounts <- combiCounts[region.rows,]
		localCombiCountMuts <- combiCountMuts[region.rows,]

		#########
		# Regularize SD of individual counts
		#########

		#Baldi & Long's formula
		bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
			sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
		}

		#a helper function to construct an appropriately sized useful pseudocount 
		# pseudocount <- function(xs) if (all(xs==0)) 1e-4 else min(xs[xs!=0],na.rm=TRUE)/10
		ps <- 1e-4

		#calculate replicate means and coefficient of variation (CV) for each condition
		# (we use CV here instead of stdev, because that's what anti-correlates with read depth
		#  so, this will be the subject of our regression below)
		mcv <- do.call(cbind,lapply(condNames,function(cond) {
			mcv <- t(apply(localCombiCounts[,condMatrix[cond,]],1,function(xs){
				m <- mean(xs,na.rm=TRUE)
				# ps <- pseudocount(m)
				m <- m+ps
				c(m, (ps+sd(xs,na.rm=TRUE))/m)
			}))
			colnames(mcv) <- paste0(cond,c(".mean",".cv"))
			mcv
		}))

		#draw scatterplots for means vs CV
		pdfFile <- paste0(outdir,"region",region.i,"_regularizationInput.pdf")
		pdf(pdfFile,6,9)
		op <- par(mfrow=c(3,2))
		with(as.data.frame(mcv),{
			hist(log10(select.mean),breaks=100,col="gray",border=NA,main="")
			hist(log10(nonselect.mean),breaks=100,col="gray",border=NA,main="")
			plot(nonselect.mean,nonselect.cv,log="x",pch=".")
			plot(nonselect.mean,select.cv,log="x",pch=".")
			plot(nonselect.mean,controlNS.cv,log="x",pch=".")
			plot(controlNS.mean,controlNS.cv,log="x",pch=".")
		})
		par(op)
		invisible(dev.off())

		exportCoefficients <- function(z) {
			coef <- coefficients(z)
			coef <- coef[order(abs(coef),decreasing=TRUE)]
			paste(mapply(function(x,y)sprintf("%s = %.02f",x,y),x=names(coef),y=coef),collapse="; ")
		}

		#build regression input matrix by transforming to logspace and adding pseudocounts
		model.in <- log10(mcv)
		model.in <- as.data.frame(cbind(
			model.in,
			log.ratio=model.in[,"select.mean"]-model.in[,"nonselect.mean"]
		))

		#for each condition, perform the regularization
		regul <- do.call(cbind,lapply(condNames,function(cond) {
			#construct the regression formula
			input.column <- paste0(cond,".cv")
			mean.column <- paste0(cond,".mean")
			#extract empirical cv
			empiric.cv <- mcv[,input.column]
			input.column.id <- which(colnames(model.in)==input.column)
			regr.formula <- as.formula(paste0(input.column,"~."))
			#identify duplicated input columns, so they can be excluded
			dupli <- setdiff(which(apply(model.in,2,function(x) {
				all(x == model.in[,input.column])
			})),input.column.id)
			.model.in <- if (length(dupli) > 0) model.in[,-dupli] else model.in
			#perform regression
			z <- lm(regr.formula,data=.model.in)
			logInfo(paste("Regression coefficents for",cond,":",exportCoefficients(z)))
			#calculate prior (log(cv)) from regressor
			prior.lcv <- predict(z)
			logInfo(sprintf("Prior PCC=%.02f",cor(10^prior.lcv,empiric.cv)))
			#calculate bayesian regularization
			observations <- length(repNames)
			bayes.cv <- bnl(pseudoObservations,observations,10^prior.lcv,empiric.cv)
			#calculate pessimistic regularization
			pessim.cv <- mapply(max,10^prior.lcv,empiric.cv)
			#calculate the 1% quantile stdev (excluding the pseudocounts) as a minimum floor
			minCut <- quantile(empiric.cv[empiric.cv > 1e-5],0.01)
			#apply that minimum flooring to the pessimisting value
			pessim.cv <- sapply(pessim.cv,function(x) if (x < minCut) minCut else x)
			#return the result
			means <- mcv[,mean.column]
			out <- cbind(
				10^prior.lcv, (10^prior.lcv)*means,
				bayes.cv, bayes.cv*means,
				pessim.cv, pessim.cv*means
			)
			colnames(out) <- paste0(cond,c(
				".prior.cv",".prior.sd",
				".bayes.cv",".bayes.sd",
				".pessim.cv",".pessim.sd"
			))
			return(out)
		}))

		combiCountsReg <- cbind(localCombiCounts,mcv,regul)


		###################
		#Plot regularization results
		###################
		pdfFile <- paste0(outdir,"region",region.i,"_regularizationRound1.pdf")
		pdf(pdfFile,6,12)
		op <- par(mfrow=c(4,2))
		for (cond in condNames) {
			topoScatter(
				log10(mcv[,paste0(cond,".cv")]),
				log10(regul[,paste0(cond,".prior.cv")]),
				resolution=80,
				xlab=expression("Empiric"~log[10](sigma/mu)),
				ylab=expression("Prior"~log[10](sigma/mu)),
				main=cond
			)
			abline(0,1,col="gray",lty="dashed")
			topoScatter(
				log10(mcv[,paste0(cond,".cv")]),
				log10(regul[,paste0(cond,".bayes.cv")]),
				resolution=80,
				xlab=expression("Empiric"~log[10](sigma/mu)),
				ylab=expression("Regularized"~log[10](sigma/mu))
			)
			abline(0,1,col="gray",lty="dashed")
		}
		par(op)
		invisible(dev.off())


		#####################
		# Quality checks
		####################

		#TODO: plot spatial distribution and histogram of sequencing error
		#


		#############
		#Filter based on high sequencing error
		# In this version, we just try to catch any straggles that weren't identified before regularization
		############
		#Song's rule: Filter out anything where the nonselect count is smaller than the WT control plus three SDs.
		#(That is, where the nonselect count could be explained by sequencing error)
		if (conservativeMode) {
			flagged <- with(combiCountsReg, controlNS.mean + 3*controlNS.pessim.sd >= nonselect.mean)
		} else {
			flagged <- with(combiCountsReg, controlNS.mean + 3*controlNS.bayes.sd >= nonselect.mean)
		}
		logInfo(sprintf(
			"Filtering out %d variants (=%.02f%%) due to likely sequencing error.",
			sum(flagged), 100*sum(flagged)/nrow(combiCountsReg)
		))
		combiCountsFiltered <- combiCountsReg[!flagged,]


		############
		#Calculate raw scores
		############

		#Helper function: Taylor approximation for variance of ratio
		approx.ratio.var <- function(mr,ms,vr,vs,covrs) {
			(mr^2)/(ms^2) * (vr/(mr^2) - 2*(covrs/(mr*ms)) + (vs/ms^2))
		}

		#pseudocounts
		c.ps <- 1e-4
		phivar.ps <- 1e-4

		rawScores <- as.df(lapply(1:nrow(combiCountsFiltered),function(i) {
			#mean numerator
			mnum <- combiCountsFiltered[i,"select.mean"]-combiCountsFiltered[i,"controlS.mean"]
			if (mnum < c.ps) mnum <- c.ps #apply flooring, so the wt control can't result in negatives counts
			#mean denominator
			mden <- combiCountsFiltered[i,"nonselect.mean"]-combiCountsFiltered[i,"controlNS.mean"]
			#variances
			if (conservativeMode) {
				#variance numerator
				vnum <- combiCountsFiltered[i,"select.pessim.sd"]^2+combiCountsFiltered[i,"controlS.pessim.sd"]^2
				#variance denominator
				vden <- combiCountsFiltered[i,"nonselect.pessim.sd"]^2+combiCountsFiltered[i,"controlNS.pessim.sd"]^2
			} else {
				vnum <- combiCountsFiltered[i,"select.bayes.sd"]^2+combiCountsFiltered[i,"controlS.bayes.sd"]^2
				vden <- combiCountsFiltered[i,"nonselect.bayes.sd"]^2+combiCountsFiltered[i,"controlNS.bayes.sd"]^2
			}
			#covariance of numerator and denominator
			covnumden <- cov(
				unlist(combiCountsFiltered[i,c("select1","select2")])
				-unlist(combiCountsFiltered[i,c("controlNS1","controlNS2")]),
				unlist(combiCountsFiltered[i,c("nonselect1","nonselect2")])
				-unlist(combiCountsFiltered[i,c("controlS1","controlS2")])
			)
			#Use helper function to estimate variance for phi
			phivar <- approx.ratio.var(mnum,mden,vnum,vden,covnumden)
			#As this is based on a Taylor approximation, we can sometimes get freaky negative values
			if (phivar < phivar.ps) phivar <- phivar.ps

			phi <- mnum/mden
			phisd <- sqrt(phivar)
			list(
				hgvsp=combiCountsFiltered[i,"hgvsp"],
				hgvsc=combiCountsFiltered[i,"hgvsc"],
				phiPrime=min(combiCountsFiltered[i,c("nonselect1","nonselect2")]),
				mean.c=mnum,
				mean.phi=phi,
				sd.phi=phisd,
				mean.lphi=log10(phi),
				sd.lphi=abs(phisd/(log(10)*phi))
			)
		}))

		##############
		# 2nd round of regularization: This time for SD of scores
		##############

		# #regression input matrix
		# splinemat <- with(rawScores,data.frame(
		# 	logsd=log10(sd.lphi),
		# 	logminbc=log10(phiPrime+1),
		# 	lphi=mean.lphi
		# ))
		# #run regression
		# z <- lm(logsd ~.,data=splinemat)
		#calculate prior
		# priorSD <- 10^predict(z)

		#regression input matrix
		splinemat <- with(rawScores,data.frame(
			logcv=log10(sd.phi/mean.phi),
			logminbc=log10(phiPrime+.1)#,
			# lphi=mean.lphi
		))
		#run regression
		z <- lm(logcv ~.,data=splinemat)
		#calculate prior
		priorSD <- 10^predict(z)*rawScores$mean.phi

		#apply regularization
		observations <- length(repNames)
		bayesSD <- bnl(pseudoObservations,observations,priorSD,rawScores$sd.phi)
		rawScores$bsd.phi=bayesSD
		rawScores$bsd.lphi=with(rawScores,abs(bsd.phi/(log(10)*mean.phi)))
		rawScores$df=if(conservativeMode) observations else observations+pseudoObservations

		#####################
		# Filter out broken values if any
		#####################
		broken <- which(is.na(rawScores$mean.lphi) | is.infinite(rawScores$mean.lphi))
		if (length(broken) > 0) {
			logWarn(sprintf("Values for %d variants could not be calculated!",
				length(broken)
			))
			rawScores <- rawScores[-broken,]
		}

		################
		# Plot results of regularization
		################
		pdfFile <- paste0(outdir,"region",region.i,"_regularizationRound2.pdf")
		pdf(pdfFile,7,7)
		op <- par(mfrow=c(2,2))
		#Plot phiPrime vs SD
		with(rawScores[rawScores$sd.lphi < 1,],topoScatter(phiPrime+1,sd.lphi,log="x",maxFreq=35,thresh=3,
			resolution=40, xlab="Non-select count (per M.)", ylab=expression(sigma)
		))
		#Plot score vs SD
		with(rawScores,topoScatter(mean.lphi,sd.lphi,log="y",pch=20,resolution=40, 
			xlab="Fitness score",ylab=expression(sigma),maxFreq=35,thresh=3
		))
		#Plot phiPrime vs regSD
		if (sum(rawScores$bsd.lphi < 1) > 1) {
			with(rawScores[rawScores$bsd.lphi < 1,],topoScatter(phiPrime+1,bsd.lphi,log="x",maxFreq=35,thresh=3,
				resolution=40, xlab="Non-select count (per M.)", 
				ylab=expression("Bayesian Regularized"~sigma)
			))
		}
		# abline(0,1,col="gray",lty="dashed")
		#Plot Empiric vs regularized
		with(rawScores,topoScatter(sd.lphi,bsd.lphi,resolution=60,maxFreq=30,log="xy",
			xlab=expression("Empiric"~sigma),ylab=expression("Bayesian Regularized"~sigma)
		))
		abline(0,1,col="gray",lty="dashed")
		par(op)
		invisible(dev.off())


		#################
		# If the functional assay is based on inverse fitness, invert the scores
		#################
		if (inverseAssay) {
			rawScores$mean.lphi <- -rawScores$mean.lphi
			#stdev remains unchanged as sqrt(((-1)^2)) = 1
		}

		#################
		# Estimate Modes of synonymous and stop
		#################

		####################3
		#TODO: Require minimum amount of filter passes rather than just 1
		#TODO: Try gaussian mixture models with two underlying distributions?
		#########################
		
		# GET SD CUTOFF FOR MODES ESTIMATE
		#if we can't find enough syn/stop below the cutoff, increase the cutoff to 1
		if (with(rawScores,sum(grepl("Ter$",hgvsp)&bsd.lphi<sdCutoff)<min_variants_to_choose_median ||
			sum(grepl("=$",hgvsp)&bsd.lphi<sdCutoff)<min_variants_to_choose_median )) {
			sdCutoff <- sdCutoffAlt
			#if we still can't find enough below the new cutoff, get rid of it altogether
			if (with(rawScores,sum(grepl("Ter$",hgvsp)&bsd.lphi<sdCutoff)<min_variants_to_choose_median ||
				sum(grepl("=$",hgvsp)&bsd.lphi<sdCutoff)<min_variants_to_choose_median )) {
				sdCutoff <- Inf
			}
		}

		# GET STOP CUTOFF FOR MODES ESTIMATE
		# funciton to get the aa position from hgvsp
		get_pos_from_hgvsp <- function(variant) {
		  variant <- unlist(strsplit(variant, split = ""))
		  nums <- grep(pattern = "[0-9]", x = variant, value = TRUE)
		  pos <- as.numeric(paste(nums, collapse = ""))
		  return(pos)
		}
		
		
		# get the positions 
		rawScores$positions <- unlist(lapply(rawScores$hgvsp, get_pos_from_hgvsp))
		# length of the gene
		length <- max(rawScores$positions)
		
		# make plot to choose stop cutoff manually
		tiffFile <- paste0(outdir,"region",region.i,"_choose_stop_cutoff_2.tiff")
		tiff(tiffFile, units="in", width=4, height=3, res=300)
		stop_cut_plot_2 <- ggplot(data = rawScores[grepl("Ter$",rawScores$hgvsp),],
		                        mapping = aes(x = positions, y =mean.lphi)) +
		  geom_point() +
		  theme_linedraw() +
		  geom_smooth() +
		  ylab(expression(log(phi))) +
		  xlab('AA position')
		print(stop_cut_plot_2)
		dev.off()
		
		
		# set the stop cutoff to the length of the gene if not provided
		if (is.null(stop_cutoff)) {
		  stop_cutoff <- length
		}
		# now we have chosen the cutoff, calculate the medians to use to scale the fitness scores
		stops <- rawScores$mean.lphi[which(grepl("Ter$",rawScores$hgvsp) 
		                                   & rawScores$bsd.lphi < sdCutoff 
		                                   & rawScores$positions <= stop_cutoff 
		                                     )]
		syns <- rawScores$mean.lphi[which(grepl("=$",rawScores$hgvsp) & rawScores$bsd.lphi < sdCutoff)]
		modes <- c(stop=median(stops,na.rm=TRUE),syn=median(syns,na.rm=TRUE)) # keep this for jochen's plots but i will use a different variable name for mine
		median_syn <- round(median(syns,na.rm=TRUE), digits = 2)
		median_stop <- round(median(stops,na.rm=TRUE), digits = 2)

		# plot only the syns and stops used in the median calculation
		tiffFile <- paste0(outdir,"region",region.i,"_syn_stop_for_scaling.tiff")
		tiff(tiffFile, units="in", width=4, height=2.5, res=300)
		hist1 <- ggplot() +
		  geom_histogram(data = data.frame(mean.lphi=stops), mapping = aes(x = mean.lphi, y = ..density..),  fill = 'darkred', alpha = 0.5, col = 'black') +
		  geom_histogram(data = data.frame(mean.lphi=syns), mapping = aes(x = mean.lphi, y = ..density..),  fill = 'darkgreen', alpha = 0.5, col = 'black') +
		  theme_linedraw() +
		  geom_vline(xintercept = median_stop, size=1, linetype=2, color='darkred') +
		  geom_vline(xintercept = median_syn, size=1, linetype=2, color='darkgreen') +
		  geom_text(aes(x=median_syn + 0.15, label=median_syn, y=6)) +
		  geom_text(aes(x=median_stop + 0.15, label=median_stop, y=6)) +
		  xlab(expression(log(phi)))
		print(hist1)
		dev.off()
		
		# plot ALL of the syns and stops and their medians
		stop.is <- which(grepl("Ter$",rawScores$hgvsp))
		syn.is <- which(grepl("=$",rawScores$hgvsp))
		miss.is <- setdiff(1:nrow(rawScores),c(stop.is,syn.is))
		all_stop_median <- round(median(rawScores$mean.lphi[stop.is]), digits = 2)
		all_syn_median <- round(median(rawScores$mean.lphi[syn.is]), digits = 2)
		
		tiffFile <- paste0(outdir,"region",region.i,"_all_syn_stop.tiff")
		tiff(tiffFile, units="in", width=4, height=2.5, res=300)
		hist2 <- ggplot() +
		  geom_histogram(data = rawScores[stop.is,], mapping = aes(x = mean.lphi, y = ..density..),  fill = 'darkred', alpha = 0.5, col = 'black') +
		  geom_histogram(data = rawScores[syn.is,], mapping = aes(x = mean.lphi, y = ..density..),  fill = 'darkgreen', alpha = 0.5, col = 'black') +
		  theme_linedraw() +
		  geom_vline(xintercept = median(rawScores$mean.lphi[stop.is]), size=1, linetype=2, color='darkred') +
		  geom_vline(xintercept = median(rawScores$mean.lphi[syn.is]), size=1, linetype=2, color='darkgreen') +
		  geom_text(aes(x=all_syn_median + 0.15, label=all_syn_median, y=6)) +
		  geom_text(aes(x=all_stop_median + 0.15, label=all_stop_median, y=6)) +
		  xlab(expression(log(phi)))
		print(hist2)
		dev.off()
		
		# plot all syns and stops but with the medians used for scaling
		tiffFile <- paste0(outdir,"region",region.i,"_all_syn_stop_with_filtered_medians.tiff")
		tiff(tiffFile, units="in", width=4, height=2.5, res=300)
		hist4 <- ggplot() +
		  geom_histogram(data = rawScores[stop.is,], mapping = aes(x = mean.lphi, y = ..density..),  fill = 'darkred', alpha = 0.5, col = 'black') +
		  geom_histogram(data = rawScores[syn.is,], mapping = aes(x = mean.lphi, y = ..density..),  fill = 'darkgreen', alpha = 0.5, col = 'black') +
		  theme_linedraw() +
		  geom_vline(xintercept = median_stop, size=1, linetype=2, color='darkred') +
		  geom_vline(xintercept = median_syn, size=1, linetype=2, color='darkgreen') +
		  geom_text(aes(x=median_syn + 0.15, label=median_syn, y=6)) +
		  geom_text(aes(x=median_stop + 0.15, label=median_stop, y=6)) +
		  xlab(expression(log(phi)))
		print(hist4)
		dev.off()
		
		# plot the missense and medians used for scaling
		tiffFile <- paste0(outdir,"region",region.i,"_all_missense.tiff")
		tiff(tiffFile, units="in", width=4, height=2.5, res=300)
		hist3 <- ggplot() +
		  geom_histogram(data = rawScores[miss.is,], mapping = aes(x = mean.lphi),  fill = 'grey', alpha = 0.5, col = 'black') +
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
			with(data,{
				stop.is <- which(grepl("Ter$",hgvsp))
				syn.is <- which(grepl("=$",hgvsp))
				miss.is <- setdiff(1:nrow(data),c(stop.is,syn.is))
				br <- seq(floor(min(mean.lphi)),ceiling(max(mean.lphi)),.1)

				hist(mean.lphi[miss.is],col="gray",breaks=br,xlab=expression(log(phi)),border=NA,main=title)
				hist(mean.lphi[stop.is],add=TRUE,col=colAlpha("firebrick3",.5),breaks=br)
				hist(mean.lphi[syn.is],add=TRUE,col=colAlpha("darkolivegreen3",.5),breaks=br)
				abline(v=modes,col=c("firebrick3","darkolivegreen3"),lty="dashed")
			})
		}

		pdfFile <- paste0(outdir,"scalingQC_region",region.i,".pdf")
		pdf(pdfFile)
		op <- par(mfrow=c(2,1))
		plotStopSyn(rawScores,"Unfiltered",modes)
		legend("topright",c("missense","synonymous","stop"),fill=c("gray","darkolivegreen3","firebrick3"))
		if (any(rawScores$bsd.lphi < sdCutoff)) {
			plotStopSyn(rawScores[rawScores$bsd.lphi < sdCutoff,],
				bquote(sigma["regularized"] < .(sdCutoff)),
				modes
			)
		}
		par(op)
		invisible(dev.off())


		#################
		# Use syn/stop medians to scale scores
		#################

		logInfo(sprintf(
			"Scaling to synonymous (log(phi)=%.02f) and nonsense (log(phi)=%.02f) medians.",modes[["syn"]],modes[["stop"]]
		))

		#if manual overrides for the synonymous and stop modes were provided, use them
		if (!is.na(region.syn)) {
			logInfo(sprintf(
				"Using manual override (=%.02f) for synonmous mode instead of automatically determined value (=%.02f).",region.syn,modes[["syn"]]
			))
			modes[["syn"]] <- region.syn
		}
		if (!is.na(region.stop)) {
			logInfo(sprintf(
				"Using manual override (=%.02f) for stop mode instead of automatically determined value (=%.02f).",region.stop,modes[["stop"]]
			))
			modes[["stop"]] <- region.stop
		}

		#apply the scaling
		denom <- modes[["syn"]]-modes[["stop"]]
		scoreMat <- with(rawScores,{
			sd <- bsd.lphi/denom
			score <- (mean.lphi - modes[["stop"]])/denom
			cbind(
				score=score,sd=sd,se=sd/sqrt(df)
			)
		})
		scores <- cbind(rawScores,scoreMat)


		##################
		# Floor negatives and fix their excessive variances
		##################

		if (!inverseAssay) {
			#the target null-like score towards which we will shift these values
			targetScore <- 0
			#the quantile for which we want to keep the p-value fixed
			quantile <- 1
			#the row numbers containing the cases to be fixed
			toFix <- which(scores$score < targetScore)
		} else {
			#if we're dealing with an inverse assay, we have to apply a ceiling instead of flooring
			#the target functional (but dead) score towards which we will shift these values
			targetScore <- 1
			#the quantile for which we want to keep the p-value fixed
			quantile <- 0
			#the row numbers containing the cases to be fixed
			toFix <- which(scores$score > targetScore)
		}
		#the equivalent sds of a normal distribution with the target mean based on the above area
		equivalent.sds <- with(scores[toFix,], sd*(quantile-targetScore)/(quantile-score))
		#apply the fixed values to the table
		scores$score.unfloored <- scores$score
		scores$sd.unfloored <- scores$sd
		scores$se.unfloored <- scores$se
		scores$score[toFix] <- targetScore
		scores$sd[toFix] <- equivalent.sds
		scores$se[toFix] <- equivalent.sds/sqrt(scores$df[toFix])
		
		logInfo(sprintf(
			"Flooring adjusted the values of %d variant scores",length(toFix)
		))


		##############################
		# Draw a histogram of standard errors in the dataset
		#############################
		pdfFile <- paste0(outdir,"errorProfile_region",region.i,".pdf")
		pdf(pdfFile,5,5)
		hist(
			log10(scores$se),col="gray",border=NA,breaks=50,
			main="Standard Error distribution",
			xlab=expression(log[10](sigma[bar(x)]))
		)
		dev.off()

		
		#################
		# Write output to file
		#################

		#detailed output of all intermediate results
		outfile <- paste0(outdir,"detailed_scores_region",region.i,".csv")
		write.csv(scores,outfile,row.names=FALSE)

		#protein-level MaveDB output
		mavedbProtScores <- scores[,c("hgvsp","score","sd","se")]
		colnames(mavedbProtScores) <- c("hgvs_pro","score","sd","se")
		outfile <- paste0(outdir,"mavedb_scores_perAA_region",region.i,".csv")
		write.csv(mavedbProtScores,outfile,row.names=FALSE)

		#nucleotide-level MaveDB output
		mavedbNuclScores <- do.call(rbind,lapply(1:nrow(scores),function(i) {
			hgvsc <- strsplit(scores[i,"hgvsc"]," ")[[1]]
			data.frame(
				hgvs_nt=hgvsc,
				hgvs_pro = scores[i,"hgvsp"],
				score = scores[i,"score"],
				sd = scores[i,"sd"],
				se = scores[i,"se"]
			)
		}))
		outfile <- paste0(outdir,"mavedb_scores_perNt_region",region.i,".csv")
		write.csv(mavedbNuclScores,outfile,row.names=FALSE)

	}

	###################
	# Join the results into a single files
	###################

	joinFiles <- function(filenameBase) {
		joined <- do.call(rbind,lapply(
			paste0(outdir,filenameBase,"_region",regions$region,".csv"),
			function(rfile) {
				if (file.exists(rfile)) {
					return(read.csv(rfile))
				} else {
					return(NULL)
				}
			}
		))
		outfile <- paste0(outdir,filenameBase,".csv")
		write.csv(joined,outfile,row.names=FALSE)
	}

	joinFiles("detailed_scores")
	joinFiles("mavedb_scores_perNt")
	joinFiles("mavedb_scores_perAA")

}

