DiagMultiKPlot_par = function(ks, res_by_k, count_vecs) {
  
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))
  suppressPackageStartupMessages(library(cowplot))
  
  # set a data frame for the data to plot
  tog <- as.data.frame(table(ks)[table(ks) > 1])
  
  # calculate rpac
  pacobj <- CalcPAC_par(x1 = 0.1, x2 = 0.9, k_vec = names(count_vecs), count_vecs = count_vecs)
  tog$rpac <- pacobj$rPAC
  tog$one_minus_rpac  <- 1-tog$rpac
  
  # Plot Freq of # of runs
  freqPlot <- ggplot(data=tog, aes(x=ks, y=Freq)) +
    geom_bar(stat="identity") +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    geom_text(aes(label=Freq),  vjust=-0.3, size=3.5) +
    scale_x_discrete("K") +
    scale_y_continuous("Number of clustering runs") +
    geom_hline(yintercept=100, linetype="dashed", color = "black")
  
  # Plot rPAC for each K
  rpacPlot <- ggplot(data=tog, aes(x=ks, y=rpac,group=1)) +
    geom_point(shape=21, color="black", fill="black", size=2) +
    geom_line() +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    scale_x_discrete("K") +
    scale_y_continuous("rPAC")
  
  # Plot (1-rPAC) Vs freq for each K
  # first find optimal K using the convex hull
  optK <- findOptK(tog)
  cat("Optimal K: ", optK)
  
  scatPlot <- ggplot(data=tog, aes(x=one_minus_rpac, y=Freq)) +
    geom_point(shape=21, color="black", fill="black", size=1.5) +
    geom_path(color="grey", alpha=0.75, linetype=2) +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    scale_x_continuous("1 - rPAC") +
    scale_y_continuous("Number of clustering runs") +
    geom_hline(yintercept=100, linetype="dashed", color = "black") +
    geom_label_repel(aes(label = ks), segment.color = 'grey50', size=3)  +
    geom_path(data=tog[match(findOptK(tog), tog$ks), ])
  
  list("plot" = plot_grid(freqPlot, rpacPlot, scatPlot, ncol=3), "tog" = tog)
  
}

mk_res_prostate <- readRDS("multik_res_prostate_no_2d_acute_0.8.rds")
res_list_prostate <- DiagMultiKPlot_par(ks = mk_res_prostate$ks, res_by_k = mk_res_prostate$res_by_k, count_vecs = mk_res_prostate$count_vecs)
write_tsv(res_list_prostate$tog, file = "multik_res_prostate_no_2d_acute.tsv")