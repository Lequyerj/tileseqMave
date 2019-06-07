devtools::install_github("jweile/hgvsParseR")
devtools::install_github("jweile/yogilog")
devtools::install_github("jweile/yogitools")
devtools::install_github("RachelSilverstein/BioTools")

source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/R/my_legacy.R")
source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/R/my_legacy_no_bottleneck.R")
source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/R/legacy.R")
source("/Users/Rachel/Desktop/Roth Lab/tileseqMave/my_convertForJoe.R")

# parameters: edit before use
#####################
countFile <- "/Users/Rachel/Desktop/Roth Lab/R_DMS/data/rawData_LDLRAP1_O.txt"
regionFile <- "/Users/Rachel/Desktop/Roth Lab/R_DMS/data/LDLRAP1_override_regionfile.txt"
# a separate outdir will be made using this base for each of te 3 scripts
# so the outputs are not confused
outDirBase <- "/Users/Rachel/Desktop/Roth Lab/R_DMS/data/LDLRAP1"
num_sd <- 3
#######################


cat("\n PARSING COUNTS USING ORIGINAL SCRIPT \n")

analyzeLegacyTileseqCounts(countFile,
                           regionFile,
                           paste0(outDirBase,"_original/"),
                           conservativeMode = T,
                           logger = NULL
                           )
my_convert_for_joe(infile = paste0(outDirBase,"_original/","detailed_scores.csv"),
                   outfile = paste0(outDirBase,"_original/","detailed_scores_joe_format.csv"),
                   df = 2)

cat("\n PARSING COUNTS USING MY NEW SCRIPT \n")

my_analyzeLegacyTileseqCounts(countFile,
                              regionFile,
                              paste0(outDirBase,"_my_script/"),
                              conservativeMode = T,
                              logger = NULL,
                              num_sd = num_sd)
my_convert_for_joe(infile = paste0(outDirBase,"_my_script/","detailed_scores.csv"),
                   outfile = paste0(outDirBase,"_my_script/","detailed_scores_joe_format.csv"),
                   df = 2)


cat("\n PARSING COUNTS USING NO BOTTLENECK SCRIPT \n")

my_analyzeLegacyTileseqCounts_no_bottleneck(countFile,
                                            regionFile,
                                            paste0(outDirBase,"_no_bottleneck/"),
                                            conservativeMode = T,
                                            logger = NULL,
                                            num_sd = num_sd)
my_convert_for_joe(infile = paste0(outDirBase,"_no_bottleneck/","detailed_scores.csv"),
                   outfile = paste0(outDirBase,"_no_bottleneck/","detailed_scores_joe_format.csv"),
                   df = 2)

