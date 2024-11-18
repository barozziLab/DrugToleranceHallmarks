
#~~~~~~~~~~~~~
#Prepare Input
##############

# scTransform normalized data
mat <- as.matrix(so@assays$SCT@data)


#~~~~~~~~~~~~~~~~~~~~~
#Compute TF-activities
######################

net_tf <- get_collectri(organism='human', split_complexes=FALSE)

# Run ulm
net_tf_acts <- run_ulm(mat=mat, 
                       net=net_tf, 
                       .source='source', 
                       .target='target',
                       .mor='mor', 
                       minsize = 5)

# Extract ulm and store it in a new assay called tfsulm in the Seurat object
so[['tfsulm']] <- net_tf_acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Save net_tf_acts and clean memory
saveRDS(object = net_tf_acts, file = acts_tf_path)
rm(net_tf_acts)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Compute Signalling Activities
##############################

#net_sig <- get_progeny(organism = 'human', top = 500)

# Run mlm
#net_sig_acts <- run_mlm(mat=mat, 
#                        net=net_sig, 
#                        .source='source', 
#                        .target='target', 
#                        .mor='weight', 
#                        minsize = 5)

# Extract mlm and store it in a new assay called pathwaysmlm in the Seurat object
#so[['pathwaysmlm']] <- net_sig_acts %>%
#  pivot_wider(id_cols = 'source', names_from = 'condition',
#              values_from = 'score') %>%
#  column_to_rownames('source') %>%
#  Seurat::CreateAssayObject(.)

# Save net_tf_acts and clean memory
#saveRDS(object = net_sig_acts, file = acts_sig_path)
#rm(net_sig_acts)


#~~~~~~~~~~~
#Save Output
############

# Save the updated Seurat object
saveRDS(object = so, file = so_decoupleR_path)

