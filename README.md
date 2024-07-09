# CharCet
 
 ![model](https://github.com/meiyichao/CharCet/blob/main/model.png)
 
 CharCet is a hybrid neural network that includes a convolutional module for sequence data and a feedforward module for predicting chromatin accessibility.
 
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
We provide detailed step-by-step instructions for running CharCet model including data preprocessing, model training, and model test.

## Data preprocessing
**Step 1**: Download scATAC-seq and scRNA-seq data annotated with cell types(h5ad files)

After downloading the data,We provide '1.Extracting_h5ad_file_information.py' to extract matrix(.csv) and metadata(.csv) information from the h5ad file,one can execute the following code to obtain the corresponding information.

```python
python 1.Extracting_h5ad_file_information.py
```
**Step 2**: Preprocess scRNA-seq data and generate expression matrix(C x N)

We provide "2. Generation_expression_matrix. R" to preprocess scRNA seq data, where the size of the expression matrix is C x N, where C represents the number of cell types and N represents the number of highly variable genes(hvg).The format of the preprocessed expression matrix is as follows:
```
	hvg_1	        hvg_2	        hvg_3	        ...     hvg_N
1	0.006225735	-0.152664427	-0.254163005	...	-0.038108164
2	-0.192960961	-0.17115587	-0.192574967	...	-0.051457183
...	...     	...     	...     	...  	...
19	-0.352929977	0.171524705	0.532515698	...	-0.065238186
```
**Step 3**: Preprocess scATAC-seq data and obtain cell type-specific peaks

This step is a preliminary preprocessing of scATAC seq data. We provide "3.Generation_celltype_specific_packs.R" to obtain cell type-specific peaks. Considering that a cell type has many cells, for each peak of that cell type, if at least 1/5 of the cells have open signals on that peak, the peak is considered chromatin accessible and retained, otherwise inaccessible and filtered. The obtained cell type specific peaks are used for subsequent analysis.

**Step 4**: Map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval)

We use the intersect function of the bedtools tool to map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval), and the mapped region is marked as "1", indicating that it is open.The code we provide for this step is "4.bedtools_intersect_class.py"
```python
python 4.bedtools_intersect_class.py
```

**Step 5**: Generating label matrix for Classification(L x C)

After obtaining the cell type-specific peaks mentioned above, we retained genomic loci with clear ATAC seq signals in at least one cell type for subsequent analysis. We provide "5.Generating_label_matrix_class.R" to generate the label matrix.The label matrix size is `L x C` where L is the number of candidate regulatory loci and C is the number of cell types.The format of the generated label matrix is as follows:
```
        	1       2       3       ...     C
region_1	0	1	0	...	0
region_2	1	1	0	...	1
...     	...    	...    	...    	...  	...
region_L	0	0	1	...	1
```
**Step 6**: Getfasta

For the remaining genomic loci mentioned above, each locus has a length of 200bp. For the prediction of each locus, we use a 1000bp DNA sequence around it. To facilitate the acquisition of this DNA sequence, we expand the upstream and downstream ranges of each locus by 400bp each (this operation only obtains a 1000bp DNA sequence, and the actual prediction range for each locus is still 200bp), and then use the getfasta function of the bedtools tool to obtain the base sequence. The provided code is "6.Getfasta.py"
```python
python 6.Getfasta.py
```

**Step 7**: Generating label matrix for Regression(L x C)

The generation of the label matrix for regression tasks is similar to that for classification tasks, and the label matrix can be generated through "7.bedtools_intersect_regress.py" and "8.Generating_label_matrix_regress.R"

## Model training and test

We provide `Classification_model.py` and `Regression_model.py` for run CharCet in a classication and regression settings, respectively.
```python
python Classification_model.py <FOLD_ID> <Sample_Ratio>
FOLD_ID: cross validation fold id, from 1-19
Sample_Ratio:Sample ratio for training,from 0-1
```
```python
python Regression_model.py <FOLD_ID> <Sample_Ratio>
FOLD_ID: cross validation fold id, from 1-19
Sample_Ratio:Sample ratio for training,from 0-1
```
```python
python Classification_test.py <FOLD_ID> <Sample_Ratio>
FOLD_ID: cross validation fold id, from 1-19
Sample_Ratio:Sample ratio for testing,from 0-1
```
```python
python Regression_test.py <FOLD_ID> <Sample_Ratio>
FOLD_ID: cross validation fold id, from 1-19
Sample_Ratio:Sample ratio for testing,from 0-1
```


























