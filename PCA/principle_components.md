# Formatting the data before inputting to R
```{bash}
grep -v "^#" sorted.zea.maf_0.05_LD_pruned_in.vcf | wc -l #67767 variants before processing, 3248-9=3239 samples
bcftools query -l sorted.zea.maf_0.05_LD_pruned_in.vcf > zea.sample.ids.txt #pulls sampleIDs in a subset file
bcftools view -S zea.sample.ids.txt -v snps sorted.zea.maf_0.05_LD_pruned_in.vcf > zea.filtered.vcf
#uncompressed sorted > converted.chr9.vcf
grep "^#CHROM" zea.filtered.vcf > snps.for.pca.input.tsv
grep -v "^#CHROM" zea.filtered.vcf >> snps.for.pca.input.tsv
sed -i 's/0|0/0/g' snps.for.pca.input.tsv
sed -i 's/0|1/1/g' snps.for.pca.input.tsv
sed -i 's/1|0/1/g' snps.for.pca.input.tsv
sed -i 's/1|1/2/g' snps.for.pca.input.tsv
sed -i 's/\.|\./NA/g' snps.for.pca.input.tsv
```

# Running principle components analysis
```{r}
library(tidyverse)
library(vroom)
library(janitor)
snps <- read_tsv("snps.for.pca.input.tsv", col_names = FALSE)
header <- snps[1, ] #extracts proper header names
colnames(snps) <- header #assigns header as the column names
snps_cleaned <- snps[-c(1:28), ] #manually gets rid of metadata rows
transposed_snps<- select(snps_cleaned, ID, starts_with("SEEDGWAS")) %>% t() %>% row_to_names(row_number = 1)
row_names <- rownames(transposed_snps) #pulls row names from the original dataset
col_names = colnames(transposed_snps) #pulls column names from the original dataset
as.data.frame(transposed_snps[1:5, 1:5]) %>% mutate_if(is.character, as.numeric) #checking if mutate_if is doing its job (turning the "0" strings into 0 actual numbers)
for.pca <- as.data.frame(transposed_snps) %>% mutate_if(is.character, as.numeric)
impute_mean = function(X){
	X_mean = colMeans(X, na.rm = T)
	X[is.na(X)] = X_mean[rep(1:ncol(X), colSums(is.na(X)))]
	X
}
for.pca<-impute_mean(as.matrix(for.pca))

pca.snps<-prcomp(for.pca, center = T, scale = F, rank. = 10)
pca.snps$x %>% as_tibble() %>%
  ggplot(aes(x = PC1, y = PC2)) + geom_point()

# calculate pct variance explained per PC
pcvar <- pca.snps$sdev^2 # variance is standard deviation squared
pcvar.pct <- tibble(pctvar=pcvar/sum(pcvar) * 100,
                    PC=1:length(pcvar))
pcvar.pct # PC1 = 2.23%, PC2 = 1.03%
write_tsv(pcvar.pct, file = "pve_for_pca.tsv")


# export PCs
pca.snps$x %>% as.data.frame() %>% add_column(Unique.ID = rownames(for.pca)) %>% 
  write_tsv(file = "zea_v2_PC_outfile.tsv")

# load metadata and wrangle
accessions.and.landrace = read_tsv(file = "/Users/xian/JRILab/names_filtered.tsv")
names_w_colon = read_lines(file= "/Users/xian/JRILab/names_w_colon.txt")
print(names_w_colon)
names_w_colon = data.frame((unlist(strsplit(names_w_colon, "\n"))))
print(names_w_colon)
accessions.and.landrace.w.colon = cbind(accessions.and.landrace, names_w_colon)
print(accessions.and.landrace.w.colon)
extra.info <- accessions.and.landrace.w.colon %>% 
  select(`X.unlist.strsplit.names_w_colon....n....`, countries_country_name, PrimaryRace)
colnames(extra.info) <- c("Target.ID", "countries_country_name", "PrimaryRace")

# load PCs and wrangle
pcs <- read_tsv("zea_v2_PC_outfile.tsv") 
print(pcs)

# merge PCs and metadata
merged <- merge(pcs, extra.info, by=c("Target.ID"))
write_tsv(merged, "/Users/xian/Desktop/JRILab_Desktop/PCA/pcs_plus_id_country_and_landrace.tsv")

#plot colored by country
merged %>% ggplot(aes(x = PC1, y = PC2, color = countries_country_name)) + geom_point()

#plot colored by landrace
pca_plot = merged %>% ggplot(aes(x = PC1, y = PC2, color = PrimaryRace)) + geom_point()
print(pca_plot)
```
