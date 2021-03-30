# Jack Humphrey
# 2019-2021
# Prepare results 

library(tidyverse)
library(here)
#library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(biomaRt)
rerun <- FALSE

createGRange <- function(coords, strands, id){
  coord_split <- as.data.frame(stringr::str_split_fixed(coords,pattern = ":|-", n = 3), stringsAsFactors = FALSE  )
  names(coord_split) <- c("chr", "start", "end")
  coord_split$start <- as.numeric(coord_split$start)
  coord_split$end <- as.numeric(coord_split$end)
  grange <- GenomicRanges::GRanges(
    seqnames = coord_split$chr, 
    ranges = IRanges::IRanges(start = coord_split$start, end = coord_split$end),
    strand = strands,
    id = id
  )
  return(grange)
}


# get ensembl id to gene name map
meta <- 
  read_tsv("~/GENCODE/gencode.v30.tx2gene.tsv.gz") %>% janitor::clean_names() %>%
  dplyr::select(ensembl_id = geneid, gene = genename) %>%
  distinct()  %>%
  mutate(ensembl_id = gsub("\\.[0-9]+", "", ensembl_id) )


## RUN DIFFERENTIAL EXPRESSION WITH DESEQ2
 
# read in gene count and convert to CPM
genes <- read_tsv("Buratti/differential-pipeline/Buratti_gene_matrix.tsv") %>%
  janitor::clean_names() %>%
  mutate(ensembl_id = gsub("\\.[0-9]+", "", ensembl_id)) %>%
  column_to_rownames(var = "ensembl_id") 

genes <- floor(genes)

library(DESeq2)

if( rerun){
  
  col_data <- tibble(sample = colnames(genes), condition = c( rep("control", 3), rep("hnRNPK", 3))) %>%
    mutate(condition = factor(condition, levels = c("control", "hnRNPK")))
  
  dds <- DESeqDataSetFromMatrix(countData = genes, colData = col_data, design = ~condition )
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
  
  resLFC <- lfcShrink(dds, coef="condition_hnRNPK_vs_control", type="apeglm")
  
  deseq_res <- resLFC %>%
    as.data.frame() %>%
    rownames_to_column(var = "ensembl_id") %>%
    left_join(meta, by = "ensembl_id") %>%
    arrange(padj)
  
  
  save(dds,res,resLFC, deseq_res, file = "Buratti/DESeq2_results.RData")

}else{
  load("Buratti/DESeq2_results.RData")
}


genes <- edgeR::cpm(genes)


### Differential splicing results

# cassette exons
exons <- read_tsv("Buratti/leafcutter-pipeline/Buratti_hnRNPK_0.001//Buratti_hnRNPK_0.001_cassette_inclusion.tsv") 

#exons <- left_join(exons, SH_SY5Y_control_expression, by = "gene")

exons <- left_join(exons, deseq_res, by = c("gene") )

exons$strand <- stringr::str_split_fixed(exons$clusterID,pattern = "_", n =  3)[,3]


# create table
# master table for all exons 

df <- dplyr::select(exons, clusterID, strand,
                gene, ensembl_id, 
                #SH_SY5Y_control_expression, Frontal_Cortex_control_expression,  
                verdict, exon_coords, intron_coords, control_inclusion, case_inclusion, delta_inclusion, ds_fdr = FDR,
                deg_logfc = log2FoldChange, deg_p = pvalue, deg_padj = padj)



df_grange <- createGRange(coords = df$exon_coords, strands = df$strand, id = df$clusterID)

df$exon_width <- width(df_grange)

# finally filter by exon width < 250bp
all_exons <- filter(df, verdict != "complex", exon_width < 250, !is.na(ensembl_id))

all_exons_grange <- df_grange[ df_grange$id %in% all_exons$clusterID]

# add additional info
all_exons <- 
  all_exons %>% 
  mutate(exon_set = case_when(
    control_inclusion < 0.1 & delta_inclusion > 0.1 ~ "Cryptic",
    control_inclusion > 0.9 & delta_inclusion < -0.1 ~ "Skiptic",
    delta_inclusion > 0.1 ~ "Included",
    delta_inclusion < -0.1 ~ "Skipped",
    TRUE ~ "null set"
  )) %>%
  mutate(deg_direction = ifelse(deg_logfc > 0, "UP", "DOWN")) %>%
  mutate(direction = ifelse(delta_inclusion > 0, "Included", "Skipped")) %>%
  mutate(verdict = gsub("cassette", "annotated", verdict)) %>%
  mutate(annotation = str_split_fixed(verdict, ":", 2)[,1] ) %>%
  mutate(annotation = case_when(
    annotation == "annotated" ~ "Annotated",
    annotation == "novelexon" ~ "Novel\nexon",
    annotation == "novelskip" ~ "Novel\nskipping"
  ))

dim(all_exons)


## export for PhyloP conservation
write_tsv(all_exons, path = "Buratti/all_cassette_metadata_final.tsv")

# export BED file
all_exons_grange$name <- all_exons_grange$id
rtracklayer::export.bed(all_exons_grange, con = "Buratti/all_cassettes.bed")

strong_exons <- filter(all_exons, exon_set %in% c("cryptic", "skiptic"))
  
intron_grange <- createGRange(coords = strong_exons$intron_coords, strands = strong_exons$strand,id = strong_exons$clusterID)
rtracklayer::export.bed(intron_grange, con = "Buratti/strong_intron_coords.bed")

# flank by 10
start(intron_grange) <- start(intron_grange) - 10
end(intron_grange) <- end(intron_grange) + 10
rtracklayer::export.bed(intron_grange, con = "Buratti/strong_intron_coords_flank10.bed")

# # positive control - Klim TDP exons
# 
# klim <- read_tsv("Klim/Klim_covariate_cassette_inclusion.tsv")
# 
# klim_intron_grange <- createGRange(coords = klim$intron_coords, strands = NA,id = klim$clusterID)
# start(klim_intron_grange) <- start(klim_intron_grange) - 10
# end(klim_intron_grange) <- end(klim_intron_grange) + 10
# rtracklayer::export.bed(klim_intron_grange, con = "Buratti/Klim_intron_coords_flank10.bed")
# 
# # better - Ward TDP knockdown in neurons
# Ward <- read_tsv("Ward/ward_tdp_control_tdp_cassette_inclusion.tsv")
# Ward$strand <- str_split_fixed(Ward$clusterID, "_", 3)[,3]
# Ward_intron_grange <- createGRange(coords = Ward$intron_coords, strands = Ward$strand ,id = Ward$clusterID)
# 
# start(Ward_intron_grange) <- start(Ward_intron_grange) - 10
# end(Ward_intron_grange) <- end(Ward_intron_grange) + 10
# rtracklayer::export.bed(Ward_intron_grange, con = "Ward/Ward_intron_coords_flank10.bed")

