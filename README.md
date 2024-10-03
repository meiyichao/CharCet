 
 ![model](https://github.com/meiyichao/CharCet/blob/main/model.png)
 
 CharCet is a deep-learning framework that integrates DNA sequences and single-cell transcription data to predict chromatin accessibility for cell types.
 
 # Requirements
 
 ## Python
- torch(2.1.1)
- scikit-learn(1.3.1)
- pyfasta(0.5.2)
- pandas(1.3.4)
- numpy(1.23.5)

 ## R
- Seurat(4.4.0)
- tidyverse(2.0.0)
- patchwork(1.1.3)
- dplyr(1.1.3)
- gtools(3.9.4)
- Matrix(1.6.5)
- data.table(1.14.8)

# Installation
CharCet can be downloaded using the following command.
```shell
git clone https://github.com/meiyichao/CharCet
```

# Instructions
In order to standardize the input data format, CharCet requires single-cell transcriptome and chromatin accessibility data to be annotated with cell types and converted into R Data Serialization (RDS) files using Seurat's saveRDS() function. RDS is a widely used binary file format for storing R objects, which is particularly convenient for subsequent analysis and sharing. The user specifies the storage directory for RDS files and passes it to the R script using the "input_data_directory" parameter. The output file directory is passed to the R script using the "output_data_directory" parameter. The detailed instructions for running the CharCet model cover data preprocessing, model training, and model testing.

## Data preprocessing
**Step 1**: Preprocess scRNA-seq data and generate expression matrix
Users are required to execute the R script Generation_expression_matrix.R using the Rscript command. This script generates expression matrices for each cell type and stores them in the final_exp_matrix.txt file.
```R
Rscript Generation_expression_matrix.R <input_data_directory> <output_data_directory>
```
The size of the expression matrix is C x N, where C represents the number of cell types and N represents the number of highly variable genes(hvg).The format of the preprocessed expression matrix is as follows:
```
	hvg_1	        hvg_2	        hvg_3	        ...     hvg_N
1	0.006225735	-0.152664427	-0.254163005	...	-0.038108164
2	-0.192960961	-0.17115587	-0.192574967	...	-0.051457183
...	...     	...     	...     	...  	...
19	-0.352929977	0.171524705	0.532515698	...	-0.065238186
```

**Step 2**: Preprocess scATAC-seq data and obtain cell type-specific peaks
Similar to Step 1, users are required to execute the R script using the Rscript command to produce accessible files in *.bed format for each cell type.
```R
Rscript Generation_celltype_specific_peaks.R <input_data_directory> <output_data_directory>
```

**Step 3**: The CharCet model utilizes a Python command to execute the bedtools_intersect.py script, enabling the mapping of cell type-specific peaks to the human reference genome. The "task" parameter within the command can be configured as either "classification" or "regression" based on the user's requirements.
```python
python bedtools_intersect.py <task> <input_data_directory> <output_data_directory>
```

**Step 4**: The CharCet model utilizes the Rscript command to generate a label matrix. For classification tasks, users are required to execute the R script Generating_label_matrix_class.R, and the resulting output will be directed to the file union.peaks.labels.class.txt. In the case of regression tasks, users should execute Generating_label_matrix_regress.R, and the output will be saved to the file union.peaks.labels.regress.txt. When the model is employed for both classification and regression tasks, it is necessary to run both scripts accordingly.
```R
Rscript Generating_label_matrix_class.R <input_data_directory> <output_data_directory>
Rscript Generating_label_matrix_regress.R <input_data_directory> <output_data_directory>
```

## Model training
CharCet offers two programs: Classification_model.py and Regression_model.py. They are used for training models in classification and regression. Users can run either one or both of them using the provided commands. The model utilizes the leave-one-out cross-validation method, where one cell type is chosen as the test set and the remaining cell types serve as the training set each time the program is executed. The cell type number used as the test set is specified by the “cell_type_id”parameter.
```python
python Classification_model.py <cell_type_id> <input_data_directory> <output_data_directory>
python Regression_model.py <cell_type_id> <input_data_directory> <output_data_directory>
```

## Model testing
Similar to model training, CharCet offers two programs for model testing: Classification_test.py and Regression_test.py, designed for testing classification and regression models, respectively. Users can run these programs using the provided commands. The test set's cell type number is determined by the "cell_type_id" parameter.
```python
python Classification_test.py <FOLD_ID> <input_data_directory> <output_data_directory>
python Regression_test.py <FOLD_ID> <input_data_directory> <output_data_directory>
```
```python
python Regression_test.py <FOLD_ID> <input_data_directory> <output_data_directory>
```
CharCet stores the predicted results for each cell type as a *.CSV file in the following format.
```
chr1:816400-817400	0.4860
chr1:816600-817600	3.2389
...     		...
chr22:50783800:50784800	0.7624
```

# Contact us

**Yichao Mei**: meiyichao0121@163.com <br>
**Junxiang Gao**: gao200@mail.hzau.edu.cn <br>

























