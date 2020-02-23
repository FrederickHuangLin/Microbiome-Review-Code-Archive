conda create -n songbird_env numpy=1.15.4 scikit-bio=0.5.5 seaborn pandas=0.23.4 -c conda-forge
conda activate songbird_env
conda install tensorflow=1.10 tqdm nomkl
conda install biom-format h5py -c conda-forge
conda install jupyter notebook
conda install songbird -c conda-forge

biom convert \
-i ../intermediates/genus_table1.txt \
-o ../intermediates/genus_table1.biom \
--table-type="OTU table" --to-hdf5

biom convert \
-i ../intermediates/genus_table2.txt \
-o ../intermediates/genus_table2.biom \
--table-type="OTU table" --to-hdf5

songbird multinomial \
	--formula "country" \
	--input-biom ../intermediates/genus_table1.biom \
	--metadata-file ../intermediates/meta1.txt \
	--summary-dir ../data/dr1

songbird multinomial \
	--formula "country" \
	--input-biom ../intermediates/genus_table2.biom \
	--metadata-file ../intermediates/meta2.txt \
	--summary-dir ../data/dr2

