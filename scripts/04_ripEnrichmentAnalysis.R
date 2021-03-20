# 20210307 alorenzetti

# description ####
# this script will download and inspect
# the COG functional classes of transcripts
# binding to SmAP1 and check for 
# enrichment

# starting processing ####
# downloading COG data for halo_nr_tx
cog = read_tsv(file = "https://alanlorenzetti.github.io/halo_nr_tx/data/cog.tsv")

# merging tibs
# and adjusting COG categories
# for unknown instances
nrtxcog = nrtx %>% 
  left_join(x = ., y = cog,
            by = "representative") %>% 
  mutate(cog_category = case_when(is.na(cog_category) ~ "Undefined",
                                  TRUE ~ as.character(cog_category)))

# adjusting COG functional categories
# there are a few genes having more than
# one functional class; for those I am
# going to keep the first one in order
# to reduce complexity
# I will also give "Mobilome..." class
# to all transposases and ISH elements
nrtxcog = nrtxcog %>% 
  mutate(cog_category = str_replace(cog_category, "\\|.*$", "")) %>% 
  mutate(cog_category = case_when(geneClass == "Transposase" ~ "Mobilome: prophages, transposons",
                                  TRUE ~ as.character(cog_category)))

# testing for enrichment ####
# lsm interaction as primary clusters
# hypergeometric enrichment test of
# COG var inside lsm interaction status
vars = "cog_category"

# setting qval threshold
qthr = 0.01

# creating list to store results
enrich = list()
inputobj = nrtxcog
funcatobj = nrtxcog
for(j in inputobj$LSmInteraction %>% unique()){
  curRegGroup = inputobj %>%
    filter(LSmInteraction == j)
  for(k in vars){
    curRegGroupVec = curRegGroup %>% 
      dplyr::select(c(k)) %>% 
      unlist(use.names = F)
    
    curRegGroupLvs = curRegGroupVec %>% 
      unique()
    
    for(l in curRegGroupLvs){
      wb = sum(curRegGroupVec == l)
      vecu = funcatobj %>% 
        dplyr::select(c(k)) %>% 
        unlist(use.names = F)
      wu = sum(vecu == l)
      bu = sum(vecu != l)
      drawn = curRegGroupVec %>% length()
      
      pval = phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
      
      tib = tibble(regRule = j,
                   criteria = k,
                   level = l,
                   pval = pval)
      enrich = bind_rows(enrich, tib)
    }
  }
}

# correcting pvalues using BH method
# filtering by pval
enrich$qval = p.adjust(enrich$pval, method = "BH")
enrich = enrich[enrich$qval < qthr,]

# plotting counts of categories with enrichment score
cogEnrichPlot = nrtxcog %>%
  filter(LSmInteraction == "Sim") %>%
  group_by(cog_category) %>% 
  summarise(count = n()) %>%
  left_join(x = .,
            y = enrich %>%
              filter(regRule == "Sim"),
            by = c("cog_category" = "level")) %>%
  mutate(enrichStatus = case_when(is.na(qval) ~ NA_character_,
                                  TRUE ~ "*"),
         cog_category = factor(cog_category, levels = rev(cog_category)),
         size = case_when(count <= 30 ~ "Menos que 30 obsevações",
                          TRUE ~ "Mais que 30 observações")) %>% 
  mutate(size = factor(size, levels = size %>% unique())) %>% 
  ggplot(aes(x = count,
             y = cog_category)) +
  geom_col(fill = "white", col = "black") +
  geom_text(aes(label = enrichStatus),
            vjust = 0.75) +
  facet_wrap(~ size,
             scales = "free_x",
  ) + 
  xlab("Observações") + 
  ylab("Categoria do COG")  +
  theme(text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "bottom")

# plotting another version of contingency table
# to place besides the enrichement chart
tnpVsBackPlot = contTbl %>% 
  ggplot(aes(x = count,
             y = geneClass,
             fill = LSmInteraction,
             label = count)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_tableau(name = "Interação SmAP1") +
  ylab("Classe") +
  xlab("Frequência Relativa") +
  theme(legend.position = "bottom")

# arranging plots
panelEnrichRelFreq = ggarrange(plotlist = list(cogEnrichPlot,
                                               tnpVsBackPlot),
                               ncol = 1,
                               nrow = 2,
                               labels = "AUTO",
                               heights = c(1, 0.5))

# saving
ggsave(filename = "plots/panelEnrichAndRelFreq.png",
       plot = panelEnrichRelFreq,
       dpi = 300,
       unit = "in",
       height = 6,
       width = 7)
