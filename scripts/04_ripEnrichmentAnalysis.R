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
# also adjusting geneClass to include
# the class of undefined by COG
nrtxcog = nrtxcog %>% 
  mutate(cog_category = str_replace(cog_category, "\\|.*$", "")) %>% 
  mutate(cog_category = case_when(geneClass == "Transposase" ~ "Mobilome: prophages, transposons",
                                  TRUE ~ as.character(cog_category))) %>% 
  mutate(geneClass = case_when(cog_category == "Undefined" & 
                                 geneClass != "Transposase" ~ "Não definida",
                               TRUE ~ as.character(geneClass)),
         geneClass = factor(x = geneClass,
                            levels = c("Não definida",
                                       "Transposase",
                                       "Outra") %>%
                              rev()))

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
         size = case_when(count <= 30 ~ "< 30 Obs.",
                          TRUE ~ "> 30 Obs.")) %>% 
  mutate(size = factor(size, levels = size %>% unique())) %>% 
  ggplot(aes(x = count,
             y = cog_category)) +
  geom_col(fill = "white", col = "black", size = 0.25) +
  geom_text(aes(label = enrichStatus),
            vjust = 0.75,
            hjust = 1.25) +
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
tnpVsBackPlot = nrtxcog %>% 
  group_by(geneClass, LSmInteraction) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(x = count,
             y = geneClass,
             fill = LSmInteraction,
             label = count)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_tableau(name = "Interação com SmAP1") +
  ylab("Classe") +
  xlab("Frequência Relativa") +
  theme(legend.position = "bottom")

# arranging plots
panelEnrichRelFreq = ggarrange(plotlist = list(cogEnrichPlot,
                                               tnpVsBackPlot),
                               ncol = 1,
                               nrow = 2,
                               labels = "AUTO",
                               heights = c(1, 0.7))

# saving
ggsave(filename = "plots/panelEnrichAndRelFreq.png",
       plot = panelEnrichRelFreq,
       dpi = 300,
       unit = "in",
       height = 6.5,
       width = 7)

# I feel like exploring the Undefined category ####
# for that I will create a wordcloud of products
prodsTxt = nrtxcog %>% 
  filter(cog_category == "Undefined") %>% 
  dplyr::select(product) %>% 
  unlist(use.names = F) %>% 
  str_replace_all("protein", "") %>% 
  str_replace_all("ORF", "") %>% 
  str_replace_all("domain", "domain") %>% 
  str_replace_all("DUF\\d+", "DUF")

# creating corpus objct
prodsTxt = Corpus(VectorSource(prodsTxt))

# creating freq matrix
termM = TermDocumentMatrix(x = prodsTxt)
termM = termM %>% as.matrix()
termdf = tibble(word = names(rowSums(termM)),
                freq = rowSums(termM)) %>% 
  filter(freq >= 10) %>% 
  arrange(freq) %>% 
  mutate(word = factor(word, levels = word))

# plotting wordcloud
# col = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`
# wc = termdf %>% 
#   filter(freq >= 10) %>% 
#   ggplot(aes(label = word,
#              size = freq)) +
#   geom_text_wordcloud() +
#   scale_size_area(max_size = 12) +
#   theme_minimal()

wc = termdf %>% 
  filter(freq >= 20) %>% 
  arrange() %>% 
  ggplot(aes(y = word,
             x = freq)) +
  geom_col(fill = "white",
           color = "black",
           size = 0.25) +
  xlab("Frequência") +
  ylab("Termo")

# comparing means
compar = list(c("Não definida", "Outra"),
              c("Não definida", "Transposase"),
              c("Outra", "Transposase"))

comparGCclasses = nrtxcog %>% 
  mutate(geneClass = factor(x = geneClass,
                            levels = c("Não definida",
                                       "Transposase",
                                       "Outra"))) %>% 
  ggplot(aes(x = geneClass,
             y = GC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = 0.5,
                   alpha = 0.2) +
  stat_compare_means(comparisons = compar,
                     method = "wilcox.test",
                     mapping = aes(label = ..p.signif..)) +
  xlab("Classe")

# arranging plots
wcGCclassespanel = ggarrange(plotlist = list(wc,
                                             comparGCclasses),
                             labels = "AUTO",
                             nrow = 1,
                             ncol = 2,
                             widths = c(1, 0.8))

# saving
ggsave(filename = "plots/termFreq_GCclasses_panel.png",
       plot = wcGCclassespanel,
       dpi = 300,
       unit = "in",
       width = 7,
       height = 4)

# contingency table for lpi
# nrtxcog %>%
#   mutate(lpiStatus = case_when(norm_lpi >= 0.5 ~ "lpi >= 0.5",
#                                norm_lpi < 0.5 ~ "lpi < 0.5",
#                                TRUE ~ NA_character_)) %>% 
#   group_by(LSmInteraction, lpiStatus) %>% 
#   summarise(count = n()) %>% 
#   drop_na()
# 
# phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
# phyper(q= 16, m= 82, n= 1724, k = 227, lower.tail = F)
