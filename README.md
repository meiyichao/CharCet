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
Rscript Generation_expression_matrix.R <input_data_directory> <output_data_directory>
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz 
│   ├── matrix.mtx.gz
```
The directory structure has three files. These three files can be obtained by processing raw single-cell sequencing data using the [cellRanger](https://www.10xgenomics.com/support/software/cell-ranger/latest) tool. Barcode contains cellular information; Features contain genes information; Matrix contains expression count information.

The output_data_directory structure is as follows:
```
├──output_data_path
│   ├── final_exp_matrix.txt
```
The directory structure has one file, which is "final_exp_matrix.txt". The size of the expression matrix is C x N, where C represents the number of cell types and N represents the number of highly variable genes(hvg).The format of the preprocessed expression matrix is as follows:
```
	hvg_1	        hvg_2	        hvg_3	        ...     hvg_N
1	0.006225735	-0.152664427	-0.254163005	...	-0.038108164
2	-0.192960961	-0.17115587	-0.192574967	...	-0.051457183
...	...     	...     	...     	...  	...
19	-0.352929977	0.171524705	0.532515698	...	-0.065238186
```

**Step 2**: Preprocess scATAC-seq data and obtain cell type-specific peaks

We provide 'Generation_celltype_specific_peaks.R' to preprocess scATAC-seq data and obtain cell type-specific peaks.
```R
Rscript Generation_celltype_specific_peaks.R <input_data_directory> <output_data_directory>
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz 
│   ├── matrix.mtx.gz
```
The directory structure has three files. These three files can be obtained by processing raw single-cell sequencing data using the cellRanger tool. Barcode contains cellular information; Features contain peaks information; Matrix contains count information.

The output_data_directory structure is as follows:
```
├──output_data_directory
│   ├── 1.bed
│   ├── 2.bed
│   ├── 3.bed
…   …
│   ├── n.bed
```
The directory structure has n files, "n" is the number of cell types. The format of the bed file is as follows:
```
chr1	 903617	   907386
chr1	 958518	   963388
...	 ...       ...    
chr22	 50625295  50629340
```
Considering that a cell type has many cells, for each peak of that cell type, if at least 1/5 of the cells have open signals on that peak, the peak is considered chromatin accessible and retained, otherwise inaccessible and filtered. The obtained cell type specific peaks are used for subsequent analysis.

**Step 3**: Map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval)

We provide 'bedtools_intersect.py' to Map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval).
```python
python bedtools_intersect.py <task> <input_data_directory> <output_data_directory>
<task>:classification or regression
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── 1.bed
│   ├── 2.bed
│   ├── 3.bed
…   …
│   ├── n.bed
```
The output_data_directory structure is as follows:
```
├──output_data_directory
│   ├── class
│   	├── 1_class.bed
│   	├── 2_class.bed
│   	├── 3_class.bed
…   	…
│   	├── n_class.bed
│   ├── regress
│   	├── 1_regress.bed
│   	├── 2_regress.bed
│   	├── 3_regress.bed
…   	…
│   	├── n_regress.bed
```
The 200bp non overlapping interval hg19 can be obtained using the "makewindows" function of the [bedtools](https://bedtools.readthedocs.io/en/latest/) tool. We use the intersect function of the bedtools tool to map cell type-specific peaks to the human reference genome of hg19 (200bp non overlapping interval), and the mapped region is marked as "1", indicating that it is open. In this step, the output directory's bed file contains an additional column of open information compared to the input directory's bed file, and all loci have a length of 200bp. The format of the output directory's bed file is as follows:
```
chr1	 0	   200	       0
chr1	 200	   400         0
...	 ...       ...         ...  
chr22	 35922000  35922200    1
...	 ...       ...         ...
```

**Step 4**: Generating label matrix for Classification(L x C)

We provide "Generating_label_matrix_class.R" to generate the label matrix of classification task.
```R
Rscript Generating_label_matrix_class.R <input_data_directory> <output_data_directory>
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── 1_class.bed
│   ├── 2_class.bed
│   ├── 3_class.bed
…   …
│   ├── n_class.bed
```
The output_data_directory structure is as follows:
```
├──output_data_directory
│   ├── union.peaks.bed
│   ├── union.peaks.labels.class.txt
```
After obtaining the bed file containing open information mentioned above, we retained genomic loci with clear ATAC-seq signals(locus is marked as "1") in at least one cell type for subsequent analysis. We provide "Generating_label_matrix_class.R" to generate the label matrix.The label matrix size is `L x C` where L is the number of candidate regulatory loci and C is the number of cell types.The format of the generated label matrix is as follows:
```
        		1       2       3       ...     C
chr1:816400-817400	0	1	0	...	0
chr1:816600-817600	1	1	0	...	1
...     		...    	...    	...    	...  	...
chr22:50783800:50784800	0	0	1	...	1
```
It should be noted that the row names of the label matrix are 1000bp long loci obtained by expanding 400bp upstream and downstream of 200bp long loci. For the extraction of DNA sequences of these 1000 bp long loci, we use the following code to generate:
```
bedtools getfasta -fi GRCh38.p14.genome.fa -bed union.peaks.bed -fo union.peaks.pad1k.fa
<-fi>:reference genome
<-bed>:the input bed file
<-fo>:fasta file containing all 1000bp Loci sequences
```
The output file of the above code is "union.peaks.pad1k.fa"

**Step 5**: Generating label matrix for Regression(L x C)

We provide "Generating_label_matrix_regress.R" to generate the label matrix of regression task.
```R
Rscript Generating_label_matrix_regress.R <input_data_directory> <output_data_directory>
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── 1_regress.bed
│   ├── 2_regress.bed
│   ├── 3_regress.bed
…   …
│   ├── n_regress.bed
```
The output_data_directory structure is as follows:
```
├──output_data_directory
│   ├── union.peaks.labels.regress.txt
```

## Model training

We provide `Classification_model.py` and `Regression_model.py` for run CharCet in a classication and regression settings, respectively.
```python
python Classification_model.py <FOLD_ID> <input_data_directory> <output_data_directory>
FOLD_ID: cross validation fold id, from 1-19
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
```python
python Regression_model.py <FOLD_ID> <input_data_directory> <output_data_directory>
FOLD_ID: cross validation fold id, from 1-19
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── final_exp_matrix.txt
│   ├── union.peaks.labels.class.txt
│   ├── union.peaks.labels.regress.txt
│   ├── union.peaks.pad1k.fa
│   ├── Leave_one_out_cross_validation.txt
```
The file 'Leave_one_out_crossnvalidation.txt' divides the cell types into training and testing sets based on leave one cross validation. The file has been provided in the train folder.

The output_data_directory structure is as follows:
```
│   ├── class
│   	├── model.pth
│   ├── regress
│   	├── model.pth
```

## Model testing

We provide `Classification_test.py` and `Regression_test.py` for testing CharCet.
```python
python Classification_test.py <FOLD_ID> <input_data_directory> <output_data_directory>
FOLD_ID: cross validation fold id, from 1-19
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
```python
python Regression_test.py <FOLD_ID> <input_data_directory> <output_data_directory>
FOLD_ID: cross validation fold id, from 1-19
<input_data_directory>:input data directory
<output_data_directory>:output data directory
```
The input_data_directory structure is as follows:
```
├──input_data_directory
│   ├── final_exp_matrix.txt
│   ├── union.peaks.labels.class.txt
│   ├── union.peaks.labels.regress.txt
│   ├── union.peaks.pad1k.fa
│   ├── Leave_one_out_cross_validation.txt
```
The file 'Leave_one_out_crossnvalidation.txt' divides the cell types into training and testing sets based on leave one cross validation. The file has been provided in the train folder.

The output_data_directory structure is as follows:
```
├──output_data_directory
│   ├── class
│   	├── 1_class.csv
│   	├── 2_class.csv
│   	├── 3_class.csv
…   	…
│   	├── n_class.csv
│   ├── regress
│   	├── 1_regress.csv
│   	├── 2_regress.csv
│   	├── 3_regress.csv
…   	…
│   	├── n_regress.csv
```
Save the prediction results in a CSV file.

The format of the classification results is as follows:
```
        		1
chr1:816400-817400	0
chr1:816600-817600	1
...     		...
chr22:50783800:50784800	0
```
The format of the regression results is as follows:
```
        		1
chr1:816400-817400	0.4860
chr1:816600-817600	1.1356
...     		...
chr22:50783800:50784800	3.2389
```

# Contact us

**Yichao Mei**: meiyichao0121@163.com <br>
**Junxiang Gao**: gao200@mail.hzau.edu.cn <br>

























