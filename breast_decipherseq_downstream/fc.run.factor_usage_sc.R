
metadata_ordered <- read.xls("metadata_filtered_ordered.xlsx") %>% 
  as_tibble() %>%
  mutate(sample = paste0(description, "_batch-", batch)) %>%
  mutate(class = ifelse(ESR1_status == "mut", "ESR1_mutant", class)) %>%
  mutate(class = factor(class, levels = c("treatment-naive", "ESR1_mutant", "first_line", "resistant", "second_line")))

H.data.long <- H.data.long %>%
  mutate(gene_program = as.numeric(str_remove(gene_program, "R33_Program")))

#re-order description based on metadata_ordered
H.data.long <- H.data.long %>%
  mutate(description = factor(description, levels = unique(metadata_ordered$description)))

#ridge plots

gene_programs <- H.data.long %>% pull(gene_program) %>% unique()

plt_ridge <- list()

for (gene_program_w in gene_programs) {
  
  d_in <- H.data.long %>% 
    dplyr::filter(gene_program == gene_program_w) %>% 
    mutate(value_sat = ifelse(value >= 1, 1, value))
  
  plt_ridge[[gene_program_w]] <- ggplot(d_in, aes(x = value_sat, y = description)) + 
    geom_density_ridges() +
    xlab("Program Activity (log10)") +
    ylab("Gene Program") +
    ggtitle(paste0("Gene Program ", gene_program_w)) +
    theme_bw()
  
}

out_path <- paste0(out_folder, "/H.genes.programs.ridgeplot.pdf")
pdf(out_path, width = 5, height = 6)
for (gene_program in gene_programs) {
  plot(plt_ridge[[gene_program]])
}
dev.off()

#NOTE: it seems the programs to have values of either <0.01 or >0.01 - classification possible!

#classify
H.data.long <- H.data.long %>%
  mutate(value_bin = value >= 0.01)

H.stats.cls <- H.data.long %>% 
  group_by(description, gene_program) %>% 
  summarise(gene_program_pos = sum(value_bin), n = n()) %>% 
  mutate(gene_program_frac = gene_program_pos / n) %>% 
  ungroup()

H.stats.cls <- H.stats.cls %>%
  left_join(metadata_ordered %>% select(description, class, treatment) %>% unique(), by = "description")

#ggplot(data = H.stats.cls, aes(x=description, y=gene_program_frac)) +
#  geom_boxplot() +
#  coord_flip() +
#  facet_wrap(~gene_program) +
#  theme_bw()

#bubble plot

H.stats.cls <- H.stats.cls %>% 
  mutate(gene_program = factor(gene_program, levels = programs_order))

plt <- ggplot(H.stats.cls, aes(x = description, y = gene_program)) + 
  geom_point(aes(size = gene_program_frac, fill = gene_program_frac), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="Fraction") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.genes.programs.cls.bubbleplot.pdf")
pdf(out_path, width = 7.5, height = 7)
plot(plt)
dev.off()

#general box plot

plt <- ggplot(data = H.stats.cls, aes(x=gene_program, y=gene_program_frac)) +
  geom_boxplot() +
  theme_bw()

out_path <- paste0(out_folder, "/H.genes.programs.cls.boxplot.pdf")
pdf(out_path, width = 7.5, height = 3)
plot(plt)
dev.off()

#ggplot(data = H.stats.cls %>% group_by(gene_program, class) %>% summarise(n = n(), frac = mean(gene_program_frac)), aes(fill=class, y=frac, x=gene_program)) + 
#  geom_bar(position="dodge", stat="identity")

#generalists vs specialist ~ how many programs used by each cell
#not great results, unfortunately 

H.stats.cells <- H.data.long %>% 
  rowwise() %>%
  mutate(rownames = paste(sample, barcode, sep = "_")) %>%
  group_by(rownames) %>% 
  summarise(gene_program_active = sum(value_bin), n = n()) %>% 
  mutate(gene_program_active_frac = gene_program_active / n) %>% 
  ungroup()

H.stats.cells <- H.stats.cells %>% 
  left_join(name_tibble, by = "rownames") %>%
  mutate(description = factor(description, levels = unique(metadata_ordered$description)))

plt <- ggplot(data = H.stats.cells, aes(x=description, y=gene_program_active)) +
  geom_boxplot() +
  coord_flip() +
  ylab("Number of active gene programs") +
  theme_bw()

out_path <- paste0(out_folder, "/H.genes.programs_active.pdf")
pdf(out_path, width = 5, height = 6)
plot(plt)
dev.off()

