# alorenzetti 20210214

# description ####
# this script will take 
# the list of genes
# interacting with lsm
# and compare that to
# the list of diff
# express genes in
# lsm knockout strain

# importing degenes from
# lsm knockout experiment
deres = list()
deres[["dfsig"]] = read.xlsx(xlsxFile = "~/gdrive/dlsm/dlsm_de_analysis_kallisto/results/DESeq2.xlsx") %>% 
  as_tibble() %>% 
  filter(pvalue < 0.01)

deres[["dfsig_locus_tags"]] = deres$dfsig %>% 
  select(locus_tag) %>% 
  unlist(use.names = F)

# finding intersection
# and plotting venn diagram
intersec = df$locus_tags_representative[df$locus_tags_representative %in% deres$dfsig_locus_tags]

# view intersection results in DE table
# deres$dfsig %>%
#   filter(deres$dfsig$locus_tag %in% intersec) %>% 
#   view()

plots$vennDiag$obj = venn(combinations = list(`Interação SmAP1` = df$locus_tags_representative,
                                              `Genes diferencialmente expressos` = deres$dfsig_locus_tags)) 
plots$vennDiag$plot = plot(plots$vennDiag$obj,
                           edges = list(col = "white"),
                           fills =  c("#E15759",
                                      "#4E79A7",
                                      "#B07AA1"),
                           legend = list(side = "bottom"))

ggsave(filename = "plots/vennDiagRIP_DE.png",
       plot = plots$vennDiag$plot,
       height = 2,
       width = 3,
       unit = "in",
       dpi = 300)
