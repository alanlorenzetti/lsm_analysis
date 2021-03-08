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
enrich$qval = p.adjust(enrich$pval)
enrich = enrich[enrich$qval < qthr,]
