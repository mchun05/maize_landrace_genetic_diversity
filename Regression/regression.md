## Filtering out samples from singleton landraces
```{bash}
cd /group/jrigrp11/mchun/SeeDs/Linear_Regression
module load bcftools
module load R

bcftools view -S two_plus_sample_ids.txt no_maf_LD_pruned_for_lm.vcf > two_plus_no_maf_LD_pruned_for_lm.vcf
```

## Fitting a model
```{r}
options(java.parameters = c("-Xmx128g", "-Xms8g"))
library(rTASSEL)
genoPathVCF <- "/path/to/no_maf_LD_pruned_for_lm.vcf"
tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
    path = genoPathVCF
)
phenoDF <- read.table("dummy_phenotype.tsv", header = TRUE)
tasGenoPhenoDF <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPathVCF,
    phenoPathDFOrObj = phenoDF,
    taxaID = "Taxon",
    attributeTypes = NULL
)

normal_tasKin <- rTASSEL::kinshipMatrix(tasObj = tasGenoPhenoDF, method = "Normalized_IBS")
r_normal_kinship_matrix <- as.matrix(normal_tasKin)
write.csv(r_normal_kinship_matrix, file = "no_maf_normal_kinship_matrix.csv", row.names = FALSE) #kept for own records/future reference

library(tidyverse)
one_hot <- read_tsv("two_plus_0_1_landraces.tsv")
one_hot_extra_column_removed <- one_hot %>% select(-Landrace)
one_hot_matrix <- as.matrix(one_hot_extra_column_removed[, -1]) 
rownames(one_hot_matrix) <- one_hot_extra_column_removed[[1]]
shared_group_counts <- one_hot_matrix %*% t(one_hot_matrix)
similarity_matrix <- (shared_group_counts > 0) * 1
diag(similarity_matrix) <- 1
similarity_matrix = as.matrix(similarity_matrix)
similarity_matrix[lower.tri(similarity_matrix)] <- NA
stacked_similarity_matrix = stack(as.data.frame(similarity_matrix))


normal_no_maf_kinship = r_normal_kinship_matrix
normal_no_maf_matrix = as.matrix(normal_no_maf_kinship)
normal_no_maf_matrix[lower.tri(normal_no_maf_matrix)] <- NA
stacked_normal_no_maf_matrix = stack(as.data.frame(normal_no_maf_matrix))
normal_no_maf_data_for_lm = cbind(stacked_normal_no_maf_matrix$values, stacked_similarity_matrix$values)

col.names = c("kinship", "landrace similarity")
colnames(normal_no_maf_data_for_lm) = col.names
normal_no_maf_data_for_lm = as_data_frame(normal_no_maf_data_for_lm) 
normal_no_maf_kinship = normal_no_maf_data_for_lm$kinship
normal_no_maf_similarity = normal_no_maf_data_for_lm$'landrace similarity'

normal_no_maf_relationship = lm(normal_no_maf_kinship ~ normal_no_maf_similarity, data = normal_no_maf_data_for_lm)
summary(normal_no_maf_relationship)
```
