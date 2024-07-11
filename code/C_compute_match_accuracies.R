rm(list=ls())
start_time <- Sys.time()

#-----------------#
# LOAD LIBRARIES  #
#-----------------#
#devtools::install("/home/users/tami/codis_snps/software/RecordMatching")
library(RecordMatching)
library(stringr)


#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the fraction of individuals that should be in the reference list", call.=FALSE)

} else if (length(args) != 0) {

  markers = as.character(args[1])
  bgl.jar = as.character(args[2])
  vcf.exe = as.character(args[3])

  dir_data_n_run = as.character(args[4])
  base.dir = as.character(args[5])

  # markers='CSF1PO,D10S1248,D12S391,D13S317,D18S51,D19S433,D1S1656,D22S1045,D2S1338'
  # base.dir='/home/users/tami/codis_snps/output/run_0.075/'
  # bgl.jar='/home/users/tami/codis_snps/software/beagle.22Jul22.46e.jar'
  # vcf.exe='/home/users/tami/codis_snps/software/vcftools_0.1.13/bin/vcftools'
  # dir_data_n_run='/home/users/tami/codis_snps/data/run_0.075/'

}

#--------#
# SETUP  #
#--------#
setup(base.dir=base.dir, bgl.jar=bgl.jar, vcf.exe=vcf.exe)


#-------------------------------#
# CALCULATE ALLELE FREQUENCIES  #
#-------------------------------#
n_markers = str_count(markers,",")+1
markers = unlist(strsplit(markers, ",", n_markers))

for (m in markers){

  ref.f = str_interp("${dir_data_n_run}/reference/${m}.recode.vcf")

  ## Calculate allele frequencies
  ref.al.freq(ref.f, m)

}

#---------------------------------#
# COMPUTE THE MATCH SCORE MATRIX  #
#---------------------------------#
cott.vec <- c(0,0,1)
llr.single.irt <- list()

for (i in 1:n_markers) {

  m = markers[i]
  print(m)
  str.f = str_interp("${dir_data_n_run}/STR/${m}_STR.GT.FORMAT")

  llr.single.irt[[i]] <- comp.match.mat(cott.vec=cott.vec, marker=m, str.f=str.f)
  names(llr.single.irt)[i] <- m
  print("Done!")

}

print("I've computed the MSM for all variants.")

mat <- Reduce('+', llr.single.irt)

print("I've created the matrix.")

#-----------------------#
# MAKE A FIGURE OF MSM  #
#-----------------------#
set.seed(312)
brks.crse <- c(-1000,-40,-30,-20,-10,0,10,20,1000)
rand.order <- sample.int(dim(mat)[1])
par(mar=c(1,1,0,0))

png(str_interp("${base.dir}/match_score_matrix.png"), width=500, height=500)
image(mat[rand.order, rand.order], axes=FALSE, asp=1, frame.plot=FALSE,
      col=rev(heat.colors(length(brks.crse) - 1)), breaks=brks.crse, xlab='STR profile', ylab='SNP Profile')
dev.off()

print("Done with computing the match score matrix!")


#---------------------------#
# COMPUTE MATCH ACCURACIES  #
#---------------------------#
print("I'm now going to calcualte match accuracies...")
match.acc <- comp.match.acc(mat)
print(match.acc)


#---------------#
# SAVE OUTPUTS  #
#---------------#
write.csv(mat, str_interp("${base.dir}/match_score_matrix.csv"))
write.csv(match.acc, str_interp("${base.dir}/match_accuracies.csv"), row.names=FALSE)

print(str_interp("Saved the match score matrix at ${base.dir}/match_score_matrix.csv."))
print(str_interp("Saved the match accuracies at ${base.dir}/match_accuracies.csv"))

# Print running time
end_time <- Sys.time()
time = end_time - start_time

print(str_interp("It took ${time} sec."))
