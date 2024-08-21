# Fetch zenodo repository Bakir et al
wget https://zenodo.org/records/7649257/files/metsFigures.zip?download=1 -P data/TRACERx421_supplement_Bakir/
# Unzip archive
unzip 'data/TRACERx421_supplement_Bakir/metsFigures.zip?download=1'
# move file
mv data/TRACERx421_supplement_Bakir/metsFigures/data/mutTableAll.cloneInfo.20220726.txt data/TRACERx421_supplement_Bakir/

# Fetch zenodo repository Frankell et al
wget https://zenodo.org/records/7822002/files/figurecode.zip?download=1 -P data/TRACERx421_supplement_Frankell/
# Unzip archive
unzip 'data/TRACERx421_supplement_Frankell/figurecode.zip?download=1'
# move files
mv data/TRACERx421_supplement_Frankell/figurecode/data/20221109_TRACERx421_mutation_table.fst data/TRACERx421_supplement_Frankell/
mv data/TRACERx421_supplement_Frankell/figurecode/data/20221109_TRACERx421_all_patient_df.rds data/TRACERx421_supplement_Frankell/
mv data/TRACERx421_supplement_Frankell/figurecode/data/20221109_TRACERx421_all_tumour_df.rds data/TRACERx421_supplement_Frankell/

