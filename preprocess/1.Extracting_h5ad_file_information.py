# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
import pandas as pd

Pbmc10k_scRNA_data = sc.read_h5ad("/home/ycmei/model_demo/data/10x-Multiome-Pbmc10k-RNA.h5ad")                         
RNA=Pbmc10k_scRNA_data.X.toarray()
pd.DataFrame(data=RNA, index=Pbmc10k_scRNA_data.obs_names, columns=Pbmc10k_scRNA_data.var_names).to_csv('Pbmc10k_RNA_matrix.csv')
pd.DataFrame(Pbmc10k_scRNA_data.obs).to_csv("Pbmc10k_RNA_metadata.csv")


Pbmc10k_scATAC_data = sc.read_h5ad("/home/ycmei/10x-Multiome-Pbmc10k-ATAC.h5ad")
ATAC=Pbmc10k_scATAC_data.X.toarray()
pd.DataFrame(data=ATAC, index=Pbmc10k_scATAC_data.obs_names, columns=Pbmc10k_scATAC_data.var_names).to_csv('Pbmc10k_ATAC_matrix.csv')
pd.DataFrame(Pbmc10k_scATAC_data.obs).to_csv("Pbmc10k_ATAC_metadata.csv")

