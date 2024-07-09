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

**Step 5**: Generating label matrix for Classification(L x C)

After obtaining the cell type-specific peaks mentioned above, we retained genomic loci with clear ATAC seq signals in at least one cell type for subsequent analysis. We provide "5.Generating_label_matrix_class.R" to generate the label matrix.The label matrix size is `L x C` where L is the number of candidate regulatory loci and C is the number of cell types.The format of the generated label matrix is as follows:
```
        	1       2       3       ...     C
region_1	0	1	0	...	0
region_2	1	1	0	...	1
...     	...    	...    	...    	...  	...
region_L	0	0	1	...	1
```
**Step 6**: Generating label matrix (L x C)

We provide scripts for generating both binary label matrix (classification) and continuous label matrix (regression) here.

The label matrix size is `L X C` where L is the number of candidate regulatory loci and C is the number of cell lines.

Use the following two scripts for generating binary label matrix
```python
python 3.3.Generate_label.py <PEAK_FILE> <CELL_SET> <OUTPUT> / 3.4.Generate_label.py <PEAK_FILE> <CELL_SET> <OUTPUT>
PEAK_FILE: the generated union peak file in `Step 3` (e.g. `union.peaks.bed`)
CELL_SET: cell id set
OUTPUT: output label matrix file
```
**Step 7**: Normalizing reads count

For reads count across different cell line, we normalize it by log transformation.
```python
python 3.5.Normalize_readscount.py <CELL_SET> <OUTPUT>
CELL_SET: cell id set
OUTPUT: output normalized reads count matrix file
```
**NOTES**: If one need to run DeepCAGE with custom data, what he/she needs to do is to generate three matrices (`TF expression matrix`, `motif score matrix` and `label matrix`) by own. 

## Model training and test

We provide `4.classification.py` and `5.Regression.py` for run DeepCAGE in a classication and regression settings, respectively.
```python
python 4.classification.py <GPU_ID> <FOLD_ID>
GPU_ID: GPU card id, default: 0
FOLD_ID: cross validation fold id, from 0-4
```
```python
python 5.Regression.py <GPU_ID> <FOLD_ID>
GPU_ID: GPU card id, default: 0
FOLD_ID: cross validation fold id, from 0-4
```
Note that the deault setting will be multi-gpu model. The trained model will be saved in `data/models` folder and prediction outcome will be saved in `data` folder.


# License
This project is licensed under the MIT License - see the LICENSE.md file for details


























