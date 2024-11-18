
data_filter <- function(d_sub_long, tf_expr_frac, sample_w, filt_tf_expr, filt_tf_act) {
  
  #filter for TFs expressed
  tf_w_expr <- tf_expr_frac %>% 
    dplyr::filter(cell_type == sample_w & fraction >= filt_tf_expr) %>% 
    pull(TF_gene)
  #filter for TFs being significantly active
  tf_w_act <- tf_stats %>% 
    dplyr::filter(cell_type == sample_w & n_cell_frac >= filt_tf_act) %>% 
    pull(source)
  tf_w <- intersect(tf_w_expr, tf_w_act)
  
  #sort the filtered TFs based on value in the sample/comparison of interest
  d_bp_w <- d_sub_long %>% 
    dplyr::filter(name == sample_w & gene %in% tf_w) %>% 
    arrange(value)
  d_sub_bp <- d_sub_long %>% 
    dplyr::filter(gene %in% d_bp_w$gene) %>% 
    mutate(gene = factor(gene, levels=d_bp_w$gene))
  
  #annotate TFs by up/down activity
  d_sub_bp_up <- d_sub_bp %>% 
    dplyr::filter(name == sample_w & value > 0.5) %>% 
    pull(gene) %>% sort()
  d_sub_bp_down <- d_sub_bp %>% 
    dplyr::filter(name == sample_w & value < -0.5) %>% 
    pull(gene) %>% sort()
  d_sub_bp <- d_sub_bp %>% 
    mutate(group = ifelse(gene %in% d_sub_bp_up, "up", ifelse(gene %in% d_sub_bp_down, "down", "else")))
  
  return(d_sub_bp)
  
}

customScatter <- function(data, x_id, x_lab, y_id, y_lab, lim) {
  
  pcc <- cor.test(data[[x_id]], data[[y_id]], method = "spearman")
  pcc_anno <- paste0("SCC = ", formatC(pcc$estimate, 2))
  
  data <- data %>% mutate(gene2 = ifelse(highlight, as.character(gene), ""))
  
  plt <- ggplot(data, aes(x=get(x_id), y=get(y_id))) + 
    geom_point() +
    xlim(lim) + 
    ylim(lim) + 
    geom_abline(slope=1, intercept=0, color = "grey") +
    geom_hline(yintercept=0, linetype=2, color = "grey") +
    geom_vline(xintercept=0, linetype=2, color = "grey") +
    ggtitle(pcc_anno) +
    xlab(x_lab) + 
    ylab(y_lab) +
    theme_bw() +
    theme(legend.position = "none")
  
  plt_anno <- ggplot(data, aes(x=get(x_id), y=get(y_id), color = highlight, label = gene2)) + 
    geom_point() +
    geom_text_repel(max.overlaps = 30) +
    scale_color_manual(values=c('darkgrey', 'red')) +
    xlim(lim) + 
    ylim(lim) + 
    geom_abline(slope=1, intercept=0, color = "grey") +
    geom_hline(yintercept=0, linetype=2, color = "grey") +
    geom_vline(xintercept=0, linetype=2, color = "grey") +
    ggtitle(pcc_anno) +
    xlab(x_lab) + 
    ylab(y_lab) +
    theme_bw() +
    theme(legend.position = "none")
  
  return(list(plt, plt_anno))
  
}
