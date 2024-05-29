# howto run it
### i.e. how i did it
----

0. prerequisites
```bash
cd /proj/snic2022-22-894/nobackup/43_amniotes/
module load bioinfo-tools conda snakemake 
conda create -n model python=3.9

conda activate model
```

1. fold data
```bash
mv chr1_simulated.columns.csv chr1_simulated.npz.columns.csv
mv chr1_derived.columns.csv chr1_derived.npz.columns.csv
mv chr1_derived.meta.csv.gz chr1_derived.npz.meta.csv.gz
mv chr1_simulated.meta.csv.gz chr1_simulated.npz.meta.csv.gz

python3 ../fold_data.py \
	-i chr1_derived.npz chr1_simulated.npz \
	-o fold_1.npz fold_2.npz fold_3.npz fold_4.npz fold_5.npz \
	-n 4 \
	-m ../data_helper.py 
```

2. train data
```bash
# five fold cross validation, by changing what is train and what is test below, 
# so run 5 times in this case once per fold, and exclude it in train

python3 ../train_model.py \
	-m ../data_helper.py \
	--train fold_1.npz fold_2.npz fold_3.npz fold_4.npz \
	--test fold_5.npz \
	--columns All \
	-c 0.1 1.0 10.0 \
	-i 100 \
	--file-pattern fold_5_[C]C_[ITER]iter.mod \
	-n 6 \
	--save-weights \
	--save-scaler fold_5.scaler.pickle
```