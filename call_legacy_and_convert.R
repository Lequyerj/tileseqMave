# devtools::install_github("jweile/hgvsParseR")
# devtools::install_github("jweile/yogilog")
# devtools::install_github("jweile/yogitools")
# devtools::install_github("RachelSilverstein/BioTools")

source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/R/my_legacy.R")
source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/R/my_legacy_no_bottleneck.R")
source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/R/legacy.R")
source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/my_convertForJoe.R")

# parameters: edit before use
#####################
countFile <- "/Users/Rachel/Desktop/Roth Lab/R_DMS/data/rawData_GDI1_2016Q20_MAVEtileseq_input.txt"
regionFile <- "/Users/Rachel/Desktop/Roth Lab/R_DMS/data/GDI1_regionfile.txt"
# a separate outdir will be made using this base for each of te 3 scripts
# so the outputs are not confused
outDirBase <- "/Users/Rachel/Desktop/Roth Lab/R_DMS/data/GDI1"
#######################

my_analyzeLegacyTileseqCounts(countFile,
                              regionFile,
                              outdir=paste0(outDirBase,"_NS200/"),
                              logger=NULL,
                              inverseAssay=FALSE,
                              pseudoObservations=2,
                              conservativeMode=TRUE, 
                              ns_filt_num_sd=3,
                              select_filt_num_sd=3, 
                              select_filt=T, 
                              min_nonselect_counts=200, # not default
                              stop_cutoff=350, # not default
                              sdCutoff=0.3, 
                              sdCutoffAlt=1, 
                              min_variants_to_choose_median=10)
my_convert_for_joe(infile = paste0(outDirBase,"_NS200/","detailed_scores.csv"),
                   outfile = paste0(outDirBase,"_NS200/","detailed_scores_joe_format.csv"),
                   df = 2)
