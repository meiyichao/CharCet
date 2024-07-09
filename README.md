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
We provide "2. Generation_expression_matrix. R" to preprocess scRNA seq data, where the size of the expression matrix is C x N, where C represents the number of cell types and N represents the number of highly variable genes(hvg). The preprocessed expression matrix is as follows:
```
	hvg_1	        hvg_2	        hvg_3	        ...     hvg_N
1	0.006225735	-0.152664427	-0.254163005	...	-0.038108164
2	-0.192960961	-0.17115587	-0.192574967	...	-0.051457183
...	...     	...     	...     	...  	...
19	-0.352929977	0.171524705	0.532515698	...	-0.065238186
```
**Step 3**: Loci filtering and candidate regulatory regions selection

Please refer to `Supplementary Figure 1` for candidate regulatory regions selection strategy. Directly run `bash 3.0.Generate_peak_bin.sh` to generate candidate regulatory regions set (`union.peaks.bed` and `union.peaks.pad1k.bed`)

**Step 4**: Generating expression matrix (N x C)

The TF gene expression matrix size is `N x C` where N is the number of TFs and C is the number of cell lines. 

```python
python 3.1.Generate_tf_exp.py <CELL_SET> <OUTPUT>
CELL_SET: cell id set
OUTPUT: output expression matrix file
```
**Step 5**: Generating motif score matrix (L x N)

The motif score matrix size is `L x N` where L is the number of candidate regulatory loci and N is the number of the coresponding TFs.

```python
python 3.2.Generate_motif_score.py <PEAK_FILE> <MOTIF_FILE> <OUTPUT>
PEAK_FILE: the generated union peak file in `Step 3` (e.g. `union.peaks.bed`)
MOTIF_FILE: motif file in homer format
OUTPUT: output motif score matrix file
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


























