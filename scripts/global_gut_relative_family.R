library(openxlsx)
library(tidyverse)
library(readr)
library(microbiome)

# Metadata
meta_data = read_tsv("data/global_gut_metadata.txt")
meta_data = meta_data %>% 
  transmute(sampleid = `#SampleID`, age = as.numeric(AGE), 
            age_group = ifelse(age <= 2, "Age \u2264 2 years old", "Age > 2 years old"),
            age_quant = cut(age, breaks = quantile(age, na.rm = T)),
            sex = SEX, country = COUNTRY)%>%
  arrange(sampleid)
meta_data = meta_data[-nrow(meta_data), ] # The last row is non-informative
meta_data$age_quant = recode(meta_data$age_quant, 
                             `(0.03,2]` = "Age \u2264 2", 
                             `(2,15]` = "2 < Age \u2264 15", 
                             `(15,33]` = "15 < Age \u2264 33",
                             `(33,83.2]` = "33 < Age \u2264 84")
meta_data$country = recode(meta_data$country, `GAZ:Malawi` = "MA", 
                           `GAZ:United States of America` = "US", `GAZ:Venezuela` = "VEN")
meta_data$country = factor(meta_data$country, levels = c("MA", "US", "VEN"))
meta_data$sex = factor(meta_data$sex, levels = c("female", "male"))
meta_data$age_group = factor(meta_data$age_group, 
                             levels = c("Age \u2264 2 years old", "Age > 2 years old"))
rownames(meta_data) = meta_data$sampleid # MANCOM-BC requires an id column

# Taxonomy
tax = read_tsv("data/global_gut_taxonomy.txt") %>% arrange(OTU_ID)
otu_id = tax$OTU_ID
tax = data.frame(tax[, -1], row.names = otu_id)
tax = apply(tax, 2, function(x) sapply(x, function(y) strsplit(y, "__")[[1]][2]))
tax = as.matrix(tax)

# OTU table
otu_table = read_tsv("data/global_gut_otu.txt")
otu_table = otu_table[, -ncol(otu_table)] # The last column is taxonomy
otu_id = otu_table$OTU_ID
otu_table = data.frame(otu_table[, -1], check.names = FALSE, row.names = otu_id)
otu_table = as.matrix(otu_table)

# OTU data
otu_data = phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), 
                    sample_data(meta_data), tax_table(tax))

# Family level
family_data = aggregate_taxa(otu_data, "Family")
family_ma = subset_samples(family_data, country == "MA")
family_us = subset_samples(family_data, country == "US")
family_ven = subset_samples(family_data, country == "VEN")

wb = createWorkbook()
addWorksheet(wb, sheetName = "MA")
writeData(wb, sheet = "MA", 
          x = t(abundances(family_ma, transform = "compositional")), 
          rowNames = TRUE)
addWorksheet(wb, sheetName = "US")
writeData(wb, sheet = "US", 
          x = t(abundances(family_us, transform = "compositional")), 
          rowNames = TRUE)
addWorksheet(wb, sheetName = "VEN")
writeData(wb, sheet = "VEN", 
          x = t(abundances(family_ven, transform = "compositional")), 
          rowNames = TRUE)

saveWorkbook(wb, file = "data/global_gut_relative_family.xlsx", overwrite = TRUE)
