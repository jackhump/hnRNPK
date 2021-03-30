# Jack Humphrey
# 2019-2021
# Conservation scoring
# requires bigwigSummary from UCSC
library(data.table)
library(ggplot2)
library(dplyr)
library(optparse)
library(stringr)

# check bigWigSummary is in the $PATH
if( Sys.which("bigWigSummary") == "" ){
  stop("bigWigSummary is not on your PATH")
}

bed.files.list <- c("Buratti/all_cassettes.bed")
bed.names <- "all exons"


species <- "human"
code <- "hnRNPK"
outFolder <- "Buratti/conservation/"
#graph_title <- "F210I and M323K extreme splicing \nexon conservation (PhyloP 60 Way)"
graph_title <- "Buratti K"

conservation.outFolder <- outFolder


# create outFolder

if (! file.exists(conservation.outFolder)) dir.create(conservation.outFolder)



# make sure phyloP bigwig exists - if not then download it!

phyloP.bw <- "Buratti/conservation/hg38.phyloP30way.bw"

if( ! file.exists(phyloP.bw) ){
  message("downloading phyloP bigwig not found! download it to root folder of repository")
  message( "rsync http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP30way/hg38.phyloP30way.bw conservation/" )
}


message("PhyloP conservation!")
# get mean PhyloP for each exon

exon.conservation <- list()
for(i in 1:length(bed.names)){
  bed.out <- bed.files.list[i]
  conservation.tmp.out <- paste0(conservation.outFolder,"/temp.txt")
  cmd <- paste0(
    "cat ",
    bed.out,
    " | while read chr start end name exonID strand; do ",
    "~/Software/bigWigSummary ",
    phyloP.bw, 
    " $chr $start $end 1; done  2>&1 | awk \' $0 ~ \"data\"{print \"NA\"}$0 !~ \"data\" {print $0}\' | sed \'/^$/d\'",
    " > ", 
    conservation.tmp.out
  )
  system(cmd)
  conservation <- data.table::fread( conservation.tmp.out )
  names(conservation)[1] <- "phyloP.score"
  
  conservation$exon.type <- bed.names[i]
  
  conservation.out <- paste0(bed.names[i],".conservation")
  
  # read in bed file to annotate the conservation scores
  d <- as.data.frame(fread(bed.files.list[i]))
  
  conservation <- cbind( conservation,d)
  
  
  assign(conservation.out,conservation)
  exon.conservation[[i]] <- conservation
}


# prepare tables for graphing
exon.conservation.merge <- as.data.frame(do.call(args = c(exon.conservation, fill = TRUE) , what = rbind, ))

all_exons <- read_tsv("Buratti/all_cassette_metadata.tsv")

all_exons$mean_phyloP_score <- exon.conservation.merge$phyloP.score[match(all_exons$clusterID, exon.conservation.merge$V4)]

write_tsv(all_exons, path ="Buratti/all_cassette_metadata_phyloP.tsv" )
