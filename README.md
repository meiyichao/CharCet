# CharCet
 
 ![model](https://github.com/meiyichao/CharCet/blob/main/model.png)
 
 CharCet is a deep-learning framework that integrates DNA sequences and single-cell transcription data to predict chromatin accessibility.
 
 # Requirements
- torch(2.1.1)
- scikit-learn(1.3.1)
- pyfasta(0.5.2)
- pandas(1.3.4)
- numpy(1.23.5)

# Installation
CharCet can be downloaded by
```shell
git clone https://github.com/meiyichao/CharCet
```

# Instructions
We provide detailed step-by-step instructions for running CharCet model including data preprocessing, model training, and model testing.

## Data preprocessing
**Step 1**: Preprocess scRNA-seq data and generate expression matrix(C x N)

We provide 'GenerationExpression_matrix.R' for generating expression matrices of cell types.
```R
Rscript Generation_expression_matrix.R <input_data_directory> <ouput_data_directory>
<input_data_directory>:input data directory
<ouput_data_directory>:ouput data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz 
│   ├── matrix.mtx.gz
…   …
│   │
```
The directory structure has three files. These three files can be obtained by processing raw single-cell sequencing data using the cellRanger tool. Barcode contains cellular information; Features contain gene information; Matrix contains expression count information.

The ouput_data_directory structure is as follows:
```
├──ouput_data_path
│   ├── final_express_matrix.txt
…   …
│   │
```
The directory structure has one file, which is "final_express_matrix.txt". The size of the expression matrix is C x N, where C represents the number of cell types and N represents the number of highly variable genes(hvg).The format of the preprocessed expression matrix is as follows:
```
	hvg_1	        hvg_2	        hvg_3	        ...     hvg_N
1	0.006225735	-0.152664427	-0.254163005	...	-0.038108164
2	-0.192960961	-0.17115587	-0.192574967	...	-0.051457183
...	...     	...     	...     	...  	...
19	-0.352929977	0.171524705	0.532515698	...	-0.065238186
```

**Step 2**: Preprocess scATAC-seq data and obtain cell type-specific peaks

This step is a preliminary preprocessing of scATAC seq data. We provide "Generation_celltype_specific_packs.R" to obtain cell type-specific peaks. Considering that a cell type has many cells, for each peak of that cell type, if at least 1/5 of the cells have open signals on that peak, the peak is considered chromatin accessible and retained, otherwise inaccessible and filtered. The obtained cell type specific peaks are used for subsequent analysis.

**Step 3**: Map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval)

We use the intersect function of the bedtools tool to map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval), and the mapped region is marked as "1", indicating that it is open.The code we provide for this step is "bedtools_intersect.py"
```python
python bedtools_intersect.py
```

**Step 4**: Generating label matrix for Classification(L x C)

After obtaining the cell type-specific peaks mentioned above, we retained genomic loci with clear ATAC seq signals in at least one cell type for subsequent analysis. We provide "Generating_label_matrix_class.R" to generate the label matrix.The label matrix size is `L x C` where L is the number of candidate regulatory loci and C is the number of cell types.The format of the generated label matrix is as follows:
```
        	1       2       3       ...     C
region_1	0	1	0	...	0
region_2	1	1	0	...	1
...     	...    	...    	...    	...  	...
region_L	0	0	1	...	1
```
**Step 5**: Generating label matrix for Regression(L x C)

The generation of the label matrix for regression tasks is similar to that for classification tasks, and the label matrix can be generated through "bedtools_intersect.py" and "Generating_label_matrix_regress.R"

## Model training

We provide `Classification_model.py` and `Regression_model.py` for run CharCet in a classication and regression settings, respectively.
```python
python Classification_model.py <FOLD_ID> 
FOLD_ID: cross validation fold id, from 1-19
```
```python
python Regression_model.py <FOLD_ID> 
FOLD_ID: cross validation fold id, from 1-19
```

## Model testing

We provide `Classification_test.py` and `Regression_test.py` for testing CharCet.
```python
python Classification_test.py <FOLD_ID> 
FOLD_ID: cross validation fold id, from 1-19
```
```python
python Regression_test.py <FOLD_ID> 
FOLD_ID: cross validation fold id, from 1-19
```

# Contact us

**Yichao Mei**: meiyichao0121@163.com <br>
**Junxiang Gao**: gao200@mail.hzau.edu.cn <br>

























