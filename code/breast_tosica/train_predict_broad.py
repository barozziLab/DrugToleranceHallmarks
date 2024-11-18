import TOSICA
import scanpy as sc
import numpy as np
import warnings 
warnings.filterwarnings ("ignore")
import torch
import anndata

# read reference data
ref = sc.read('reference_data/epithelial_cells.h5ad')

# define groups of cell types
author_cell_type_groups = {
    "basal": "Basal",
    "LummHR-major": "LumHR",
    "LummHR-active": "LumHR",
    "LummHR-SCGB": "LumHR",
    "Lumsec-major": "LumSec",
    "Lumsec-basal": "LumSec",
    "Lumsec-HLA": "LumSec",
    "Lumsec-KIT": "LumSec",
    "Lumsec-lac": "LumSec",
    "Lumsec-myo": "LumSec",
    "Lumsec-prol": "LumSec"
}

ref.obs['author_cell_type_group'] = ref.obs['author_cell_type'].map(author_cell_type_groups)

# change var_names to gene symbols
ref.var.ensembl = ref.var_names
ref.var_names = ref.var.feature_name
# pre-process reference data
sc.pp.normalize_total(ref, target_sum=1e4)
sc.pp.log1p(ref)
sc.pp.highly_variable_genes(ref)
ref = ref[:, ref.var.highly_variable]

ref = ref[:,ref.var_names]
print(ref)
print(ref.obs.author_cell_type.value_counts())

# read query data
query = sc.read('MCF7_full.h5ad')

# pre-process query data
sc.pp.normalize_total(query, target_sum=1e4)
sc.pp.log1p(query)
sc.pp.highly_variable_genes(query)

# method 1: all features that are highly variable in ref and exist in query
overlap_hvg = query.var.features[query.var.features.isin(ref.var.feature_name)]

query_1 = query[:, query.var.features.isin(overlap_hvg)]
ref_1 = ref[:, ref.var.feature_name.isin(overlap_hvg)]
query_1 = query_1[:,ref_1.var_names]

#train
TOSICA.train(ref_1, gmt_path='human_gobp', label_name='author_cell_type_group', project='A100/hBreast_epithelial_hr_sec_method_1_full')
#predict


# method 2: all features that are highly variable in ref and highly variable in query
# filter query
query_2 = query[:, query.var.highly_variable]

overlap_hvg = query_2.var.features[query_2.var.features.isin(ref.var.feature_name)]

query_2 = query[:, query.var.features.isin(overlap_hvg)]
ref_2 = ref[:, ref.var.feature_name.isin(overlap_hvg)]
query_2 = query_2[:,ref_2.var_names]

#train
TOSICA.train(ref_2, gmt_path='human_gobp', label_name='author_cell_type_group', project='A100/hBreast_epithelial_hr_sec_method_2_full')

# method 3: all features that are highly variable in ref, fill in non-existing ones in query
missing_genes = ref.var.feature_name[~ref.var.feature_name.isin(query.var.features)]

# Create a new anndata object to store the missing genes with counts set to 0
missing_genes_data = anndata.AnnData(
    X=np.zeros((len(query.obs), len(missing_genes))),  # Initialize with zeros
    obs=query.obs,  # Use the same observations as the query object
    var=ref.var.loc[missing_genes],  # Extract rows corresponding to missing genes
)

query_3 = anndata.concat([query, missing_genes_data], axis = 1)
ref_3 = ref
query_3 = query_3[:,ref_3.var_names]

#train
TOSICA.train(ref_3, gmt_path='human_gobp', label_name='author_cell_type_group', project='A100/hBreast_epithelial_hr_sec_method_3_full')

# PREDICTIONS

# method 1:
model_weight_path = 'A100/hBreast_epithelial_hr_sec_method_1_full/model-9.pth'
predictions_1 = TOSICA.pre(query_1, model_weight_path = model_weight_path, project='A100/hBreast_epithelial_hr_sec_method_1_full')

# predictions_1.__dict__['_raw'].__dict__['_var'] = predictions_1.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'}) # bugfix to save the object
predictions_1.obs.to_csv('A100/predictions_hr_sec_1.txt', sep='\t', index=True)

# method 2:
model_weight_path = 'A100/hBreast_epithelial_hr_sec_method_2_full/model-9.pth'
predictions_2 = TOSICA.pre(query_2, model_weight_path = model_weight_path, project='A100/hBreast_epithelial_hr_sec_method_2_full')

# predictions_2.__dict__['_raw'].__dict__['_var'] = predictions_2.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'}) # bugfix to save the object
predictions_2.obs.to_csv('A100/predictions_hr_sec_2.txt', sep='\t', index=True)

# method 3:
model_weight_path = 'A100/hBreast_epithelial_hr_sec_method_3_full/model-9.pth'
predictions_3 = TOSICA.pre(query_3, model_weight_path = model_weight_path, project='A100/hBreast_epithelial_hr_sec_method_3_full')

# predictions_3.__dict__['_raw'].__dict__['_var'] = predictions_3.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'}) # bugfix to save the object
predictions_3.obs.to_csv('A100/predictions_hr_sec_3.txt', sep='\t', index=True)