conda activate qiime2

# Step 1: Import data
qiime tools import \
  --input-path ../intermediates/genus_table1.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ../intermediates/genus_table1.qza

qiime tools import \
 --type 'FeatureData[Taxonomy]' \
 --input-path ../intermediates/tax1.txt \
 --input-format TSVTaxonomyFormat \
 --output-path ../intermediates/tax1.qza

qiime tools import \
  --input-path ../intermediates/genus_table2.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ../intermediates/genus_table2.qza

qiime tools import \
 --type 'FeatureData[Taxonomy]' \
 --input-path ../intermediates/tax2.txt \
 --input-format TSVTaxonomyFormat \
 --output-path ../intermediates/tax2.qza

# Step 2: Creating balances
qiime gneiss correlation-clustering \
  --i-table ../intermediates/genus_table1.qza \
  --o-clustering ../intermediates/hierarchy1.qza

qiime gneiss correlation-clustering \
  --i-table ../intermediates/genus_table2.qza \
  --o-clustering ../intermediates/hierarchy2.qza

# Step 3: Building linear models using balances
qiime gneiss ilr-hierarchical \
  --i-table ../intermediates/genus_table1.qza \
  --i-tree ../intermediates/hierarchy1.qza \
  --o-balances ../intermediates/balances1.qza

qiime gneiss ols-regression \
  --p-formula "country" \
  --i-table ../intermediates/balances1.qza \
  --i-tree ../intermediates/hierarchy1.qza \
  --m-metadata-file ../intermediates/meta1.txt \
  --o-visualization ../figures/gneiss_summary1.qzv

qiime gneiss balance-taxonomy \
  --i-table ../intermediates/genus_table1.qza \
  --i-tree ../intermediates/hierarchy1.qza \
  --i-taxonomy ../intermediates/tax1.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y0' \
  --m-metadata-file ../intermediates/meta1.txt \
  --m-metadata-column country \
  --o-visualization ../figures/y0_taxa_summary1.qzv

qiime gneiss ilr-hierarchical \
  --i-table ../intermediates/genus_table2.qza \
  --i-tree ../intermediates/hierarchy2.qza \
  --o-balances ../intermediates/balances2.qza

qiime gneiss ols-regression \
  --p-formula "country" \
  --i-table ../intermediates/balances2.qza \
  --i-tree ../intermediates/hierarchy2.qza \
  --m-metadata-file ../intermediates/meta2.txt \
  --o-visualization ../figures/gneiss_summary2.qzv

qiime gneiss balance-taxonomy \
  --i-table ../intermediates/genus_table2.qza \
  --i-tree ../intermediates/hierarchy2.qza \
  --i-taxonomy ../intermediates/tax2.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y0' \
  --m-metadata-file ../intermediates/meta2.txt \
  --m-metadata-column country \
  --o-visualization ../figures/y0_taxa_summary2.qzv







