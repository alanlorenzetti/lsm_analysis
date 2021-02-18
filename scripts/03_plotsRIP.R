# alorenzetti 20210213

# description ####
# this script will plot
# RIP results and
# perform enrichment
# analysis by type of feature

# starts here ####
# GC comparison ####
# getting GC distributions for those
# features binding to LSm vs non binding
# wilcoxon unpaired two sample test
# alt name: Mann–Whitney U test
nrtx %>% 
  ggplot(aes(y = GC,
             x = LSmInteraction)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(size = 0.5,
                alpha = 0.2) +
  stat_compare_means(method = "wilcox.test",
                     mapping = aes(label = ..p.signif..),
                     label.x = 1.5) +
  xlab(label = "Interação com SmAP1")

# contingency table and phyper test ####
# graphical representation of a contingency table
contTbl = nrtx %>% 
  group_by(LSmInteraction,
           geneClass) %>% 
  summarize(count = sum(n()))

contTbl %>% 
  ggplot(aes(x = LSmInteraction,
             y = geneClass,
             label = count
             )) +
  geom_point(size = 20,
             colour = "black",
             shape = 21,
             fill = "white") +
  geom_text(size = 5) +
  xlab("Interação com SmAP1") +
  ylab("Classe")

# hypergeometric test
# rationale from
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
blackdrawn = contTbl %>%
  ungroup() %>% 
  filter(LSmInteraction == "Sim" & geneClass == "Transposase") %>% 
  select(count) %>% 
  unlist(use.names = F)

whiteurn = contTbl %>%
  ungroup() %>%
  filter(LSmInteraction == "Não") %>% 
  select(count) %>% 
  unlist(use.names = F) %>% 
  sum()

blackurn = contTbl %>%
  ungroup() %>%
  filter(LSmInteraction == "Sim") %>% 
  select(count) %>% 
  unlist(use.names = F) %>% 
  sum()

drawn = contTbl %>%
  ungroup() %>%
  filter(geneClass == "Transposase") %>% 
  select(count) %>% 
  unlist(use.names = F) %>% 
  sum()
  
# the test confirms transposases are enriched for LSm binding
phyper(q = blackdrawn,
       m = whiteurn,
       n = blackurn,
       k = drawn)
