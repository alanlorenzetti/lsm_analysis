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

venn(combinations = list(int = df$locus_tags_representative,
                         de = deres$dfsig_locus_tags)) %>% 
  plot()
