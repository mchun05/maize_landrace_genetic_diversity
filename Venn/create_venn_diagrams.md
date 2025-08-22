## Creating shared variant counts file
The code for the following chunk was taken from https://aabiddanda.github.io/geovar/notebooks/getting-started.html and edited to match the files used for analysis
```{python}
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pkg_resources
from geovar import *

# Filepath to the VCF File
vcf_file = "/group/jrigrp11/mchun/SeeDs/MAF_0.10/sorted.zea.maf_0.10_LD_pruned_in.vcf"

# Reading the population dataframe
pop_df = pd.read_csv("/group/jrigrp11/mchun/SeeDs/MAF_0.10/maf_0.10_filtered_pop_panel.txt", sep="\t")

# Writing out VCF to a Frequency Table
af_df = vcf_to_freq_table(vcf_file, pop_df=pop_df, outfile="/group/jrigrp11/mchun/SeeDs/MAF_0.10/maf.0.10.freq.csv".format(data_path), minor_allele=True)

# Print the beginning of the allele frequency table
af_df.head()

# Creating the GeoVar Object
geovar = GeoVar()

# Adding in the frequency file (all of it)
geovar.add_freq_mat(freq_mat_file="/group/jrigrp11/mchun/SeeDs/MAF_0.10/maf.0.10.freq.csv")

# Generate a geovar binning with the binning we used in our paper
geovar.bins = [(0, 0), (0, 0.10),(0.10, 1.0)]
geovar.geovar_binning()

# Output counts to file
counts_output = geovar.count_geovar_codes()
geovar_code_counts = pd.DataFrame({"codes":counts_output[0],"counts":counts_output[1]})
geovar_code_counts.to_csv("/group/jrigrp11/mchun/SeeDs/MAF_0.10/maf_0.10_geovar_code_counts.txt", sep=' ', header=False, index=False)
```
## Creating diagrams with the shared variant counts
Functions used to create the diagrams are a subset from James Kitchens' https://github.com/kitchensjn/visualizing-human-genetic-diversity/blob/main/create_euler_diagrams.R
```{r}
library(eulerr)
library(jsonlite)

sum_snps <- function(pop_combo, pops, counts, cutoff) {
  pop_tf <- rep(FALSE, length(pops))
  pop_tf[pop_combo] <- TRUE
  populations <- pops[pop_combo]
  if (length(pop_combo) > 1) {
    shared_common_snps <- which(rowSums(counts[,populations]>=cutoff)==length(populations))
  } else {
    shared_common_snps <- which(counts[,populations]>=cutoff)
  }
  snp_counts <- data.frame(t(c(pop_tf, sum(counts[shared_common_snps, ncol(counts)]))))
  colnames(snp_counts) <- c(pops, "shared_common_snps")
  return(snp_counts)
}

generate_euler_plot <- function(pairwise_file, poplist, selected_pops=c(), cutoff=3, common_pop="") {
  data <- read.table(
    file=pairwise_file,
    col.names=c("geovar_code", "counts"),
    colClasses=c("character", "numeric")
  )
  pops <- read.table(file=poplist)$V1
  if (length(selected_pops) < 1) {
    selected_pops <- pops
  }
  
  split.pops <- sapply(data[,1], function(a){strsplit(as.character(a),"")[[1]]})
  split.pops <- apply(split.pops,2,as.numeric)
  split.pops <- as.data.frame(t(split.pops))
  row.names(split.pops) <- 1:nrow(split.pops)
  colnames(split.pops) <- c(as.character(pops))
  data <- cbind(split.pops, data)
  data <- data[,c(selected_pops,"geovar_code","counts")]
  
  if (common_pop!=""){
    data <- data[which(data[common_pop]>=cutoff),]
  }
  
  combo_func <- Map(combn, list(1:length(selected_pops)), seq_along(1:length(selected_pops)), simplify=FALSE)
  combinations <- unlist(combo_func, recursive=FALSE)
  
  shared_common_snps <- do.call(rbind,lapply(combinations, FUN=sum_snps, pops=selected_pops, counts=data, cutoff=cutoff))
  shared_common_snps[nrow(shared_common_snps),"unique_snps"] <- shared_common_snps[nrow(shared_common_snps),"shared_common_snps"]
  
  for (i in (nrow(shared_common_snps)-1):1) {
    shared_common_snps[i,"unique_snps"] <- shared_common_snps[i,"shared_common_snps"] - sum(
      shared_common_snps[
        apply(
          shared_common_snps[,1:(ncol(shared_common_snps)-2)] - shared_common_snps[rep(i,nrow(shared_common_snps)),1:(ncol(shared_common_snps)-2)], 
          MARGIN=1, 
          FUN=function(x){!any(x<0)}
        ),"unique_snps"
      ], na.rm=TRUE
    )
  }
  
  euler_data <- shared_common_snps$unique_snps
  names(euler_data) <- apply(shared_common_snps[selected_pops], 1, function(x){paste(names(x)[which(x==1)],collapse="&")})
  return(list("sets"=euler_data,"euler"=eulerr::euler(euler_data, shape = "ellipse"),"size"=sum(shared_common_snps$unique_snps)))
}

geovar_counts = read_file("/Users/xian/JRILab/maf_0.05_fifteen_samples_geovar_code_counts.txt")
df_geovar_counts = data.frame((unlist(strsplit(geovar_counts, "\n"))))
print(df_geovar_counts)
pop_file = read_file("/Users/xian/JRILab/fifteen_plus_pop_panel.txt")



pops = c("TUXPEN", "OLOTIL", "TEPECI") #chosen because Jeff believes that tuxpeno should have lots of variation because it's so widespread and Olotillo and tepecintle is hypothesized to be the parents of the hybrid tuxpeno

the_plot <- generate_euler_plot("/Users/xian/JRILab/fifteen_geovar_code_counts.txt", "/Users/xian/JRILab/fifteen_plus_pop_panel.txt", selected_pops=pops, cutoff=1)
plot(your_plot$euler, quantities = TRUE, fills = c("#F0E442", "#0072B2", "#CC79A7"), edges = list(
  col = c("black", "black", "black"),        # consistent for clarity
  lwd = 5,
  lty = c("solid", "dashed", "dotted")       # distinguishable patterns
), main = "Hypothesized Parents (OLOTIL & TEPECI) and Hybrid (TUXPEN); GLOBAL MAF: 0.05")



v2.pops = c("TUXPEN", "COMITE", "REVENT") # chosen as a comparison to pops, pops should be unrelated if graphs looks the same it provides evidence to disprove the hybrid/parent theory

v2_the_plot <- generate_euler_plot("/Users/xian/JRILab/fifteen_geovar_code_counts.txt", "/Users/xian/JRILab/fifteen_plus_pop_panel.txt", selected_pops=v2.pops, cutoff=1)
plot(your_plot$euler, quantities = TRUE, fills = c("#E69F00", "#56B4E9", "#1EB489"), edges = list(
  col = c("black", "black", "black"),        # consistent for clarity
  lwd = 5,
  lty = c("solid", "dashed", "dotted")       # distinguishable patterns
),  main = "TUXPEN, COMITE, and REVENT (no hypothesized connection); GLOBAL MAF: 0.05")
```
