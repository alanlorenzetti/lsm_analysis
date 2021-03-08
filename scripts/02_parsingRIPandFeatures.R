# alorenzetti 20210211

# description ####
# this script will load files
# cross lsm results and 
# and plot charts

# start analysis ####
# loading interaction tracks
int = list()

# creating results dir
if(!dir.exists("results")){dir.create("results")}

# downloading
# halo non redundant transcriptome
nrtx = read_tsv("https://alanlorenzetti.github.io/halo_nr_tx/data/dictionary.tsv")

# loading and parsing halo annotation
annot = rtracklayer::import("data/Hsalinarum-gene-annotation-pfeiffer2019.gff3")
genes = subset(annot, type == "gene" | type == "pseudogene")
names(genes) = genes$locus_tag

# loading interaction regions generated
# using biorep1 vs control
int[["biorep1"]] = c(
  rtracklayer::import("data/biorep1-interaction-regions-entire-genome-fwd.gff3"),
  rtracklayer::import("data/biorep1-interaction-regions-entire-genome-rev.gff3")
)

# using biorep2 vs control
int[["biorep2"]] = c(
  rtracklayer::import("data/biorep2-interaction-regions-entire-genome-fwd.gff3"),
  rtracklayer::import("data/biorep2-interaction-regions-entire-genome-rev.gff3")
) 

# getting what genes intersect with interaction regions
overlaps = list()

# interaction regions of biorep1
overlaps[["biorep1"]] = genes[
  findOverlaps(query = genes, subject = int[["biorep1"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# interaction regions of biorep2
overlaps[["biorep2"]] = genes[
  findOverlaps(query = genes, subject = int[["biorep2"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# finding the union list of locus_tags
overlaps[["union"]] = c(overlaps[["biorep1"]],
                        overlaps[["biorep2"]]) %>% 
  unique()

# doing the same but trying to find
# interaction regions on the opposite strand
# interaction regions of biorep1 antisense
overlaps[["biorep1_as"]] = genes[
  findOverlaps(query = genes %>% invertStrand(), subject = int[["biorep1"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# interaction regions of biorep2 antisense
overlaps[["biorep2_as"]] = genes[
  findOverlaps(query = genes %>% invertStrand(), subject = int[["biorep2"]]) %>%
    as_tibble() %>% 
    select(queryHits) %>% 
    unique() %>% 
    unlist(use.names = F)
] %>% names()

# finding the union list of locus_tags antisense
overlaps[["union_as"]] = c(overlaps[["biorep1_as"]],
                        overlaps[["biorep2_as"]]) %>% 
  unique()

# creating an interaction tibble and
# getting the representative locus_tag
# for each one of the locus_tags within the union list
df = list()
df[["union"]] = tibble(locus_tag = overlaps[["union"]],
                       LSmInteraction = "Sim")

df[["redundant_df"]] = nrtx %>%
  separate_rows(locus_tag, sep = ",") %>% 
  left_join(x = ., y = df[["union"]], by = "locus_tag") %>% 
  drop_na()

df[["only_representative"]] = df[["redundant_df"]] %>% 
  dplyr::select(-locus_tag) %>% 
  unique()

df[["locus_tags_representative"]] = df[["only_representative"]] %>%
  select(representative) %>%
  unlist(use.names = F)

# doing the same for antisense
df[["union_as"]] = tibble(locus_tag = overlaps[["union_as"]],
                       LSmInteractionAS = "Sim")

df[["redundant_df_as"]] = nrtx %>%
  separate_rows(locus_tag, sep = ",") %>% 
  left_join(x = ., y = df[["union_as"]], by = "locus_tag") %>% 
  drop_na()

df[["only_representative_as"]] = df[["redundant_df_as"]] %>% 
  dplyr::select(-locus_tag) %>% 
  unique()

df[["locus_tags_representative_as"]] = df[["only_representative_as"]] %>%
  select(representative) %>%
  unlist(use.names = F)

# nonredundant list per biorep
overlaps[["biorep1_nr"]] = overlaps[["biorep1"]][overlaps[["biorep1"]] %in% df[["locus_tags_representative"]]]
overlaps[["biorep2_nr"]] = overlaps[["biorep2"]][overlaps[["biorep2"]] %in% df[["locus_tags_representative"]]]

# nonredundant list per biorep
overlaps[["biorep1_nr_as"]] = overlaps[["biorep1_as"]][overlaps[["biorep1_as"]] %in% df[["locus_tags_representative_as"]]]
overlaps[["biorep2_nr_as"]] = overlaps[["biorep2_as"]][overlaps[["biorep2_as"]] %in% df[["locus_tags_representative_as"]]]

# adding info to nrtx object
nrtx = nrtx %>% 
  left_join(x = .,
            y = df[["only_representative"]] %>% select(-product),
            by = "representative") %>% 
  mutate(LSmInteraction = case_when(is.na(LSmInteraction) ~ "Não",
                                    TRUE ~ as.character(LSmInteraction))) %>% 
  left_join(x = .,
            y = df[["only_representative_as"]] %>% select(-product),
            by = "representative") %>% 
  mutate(LSmInteractionAS = case_when(is.na(LSmInteractionAS) ~ "Não",
                                    TRUE ~ as.character(LSmInteractionAS)),
         geneClass = case_when(str_detect(product, "ISH|transposase") ~ "Transposase",
                               TRUE ~ "Outra"))
  
# adding GC of transcriptional units ####
# to nrtx object
# reading genome file
genome = readDNAStringSet("data/Hsalinarum.fa")
names(genome) = str_replace(names(genome), " .*$", "")  

# extracting sequences
seqs = getSeq(genome, genes)

# computing GC
# and creating GC tibble
gc = tibble(locus_tag = names(seqs),
            GC = letterFrequency(seqs, letters = "GC", as.prob = T) %>%
              as.numeric())
meangc = gc$GC %>%
  mean()

# merging dfs
nrtx = nrtx %>% 
  left_join(x = ., y = gc, by = c("representative" = "locus_tag")) %>% 
  mutate(GCdev = GC-meangc)

# saving object
write_tsv(x = nrtx,
          file = "results/interactionListNRTX.tsv")

write.xlsx(x = nrtx %>% 
             select(representative,
                    product,
                    locus_tag,
                    LSmInteraction,
                    LSmInteractionAS),
           file = "results/interactionListNRTX_basic.xlsx")

# parsing IS dataframe ####
# in order to classify
# each non redundant ORF according
# to its respective IS
isannot = rtracklayer::import("data/Hsalinarum-pfeiffer2019-mobileElements.gff3")

# finding locus_tags of transposases
# interacting with LSm
ltlsm = list()
ltlsm[["locus_tags"]] = nrtx %>%
  filter(LSmInteraction == "Sim" & geneClass == "Transposase") %>%
  select(representative) %>%
  unlist(use.names = F)

# subsetting from gene annot dataset
ltlsm[["annot"]] = genes[genes$locus_tag %in% ltlsm$locus_tags,]

ishits = findOverlaps(query = ltlsm$annot,
                      subject = isannot,
                      ignore.strand = T,
                      minoverlap = 50) %>% 
  as_tibble()

results = tibble(IS = isannot$mobile_element_type[ishits$subjectHits] %>% 
                   str_replace(., ".*:", ""),
                 locus_tag = ltlsm$locus_tags[ishits$queryHits])
results = results[mixedorder(results$IS),]
