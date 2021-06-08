# alorenzetti 20210420

# description ####
# this script will load files and
# get hfq interaction results for
# salmonella enterica

# start analysis ####
# loading interaction tracks
int = list()

# loading and parsing annotation
annot = rtracklayer::import("data/Senterica.gff")
genes = subset(annot, type == "gene" | type == "pseudogene")
names(genes) = genes$old_locus_tag

# loading interaction regions generated
int[["biorep1"]] = c(
  rtracklayer::import("data/senterica-interaction-regions-entire-genome-fwd.gff3"),
  rtracklayer::import("data/senterica-interaction-regions-entire-genome-rev.gff3")
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

# removing genes that don't have
# an old_locus_tag
overlaps$biorep1 = overlaps$biorep1[!is.na(overlaps$biorep1)]

# parsing orthologues from OrtholugeDB
orts = read_tsv(file = "http://www.pathogenomics.sfu.ca/ortholugedb/paired/download/plain?strain1=Halobacterium%20sp.%20NRC-1&strain2=Salmonella%20enterica%20subsp.%20enterica%20serovar%20Typhimurium%20str.%20SL1344")
orts = orts %>%
  select(Hsalinarum = `Locus Tag (Strain 1)`,
         Senterica = `Locus Tag (Strain 2)`)

# adding interacting genes to orts list
orts = left_join(x = orts,
                 y = tibble(locus_tag = overlaps$biorep1,
                            hfq_int = "Sim"),
                 by = c("Senterica" = "locus_tag")) %>% 
  mutate(hfq_int = case_when(is.na(hfq_int) ~ "Não",
                             TRUE ~ hfq_int))

# preparing to merge with nrtx obj
res = left_join(x = nrtxExp,
                y = orts,
                by = c("locus_tag" = "Hsalinarum")) %>% 
  drop_na() %>% 
  select(-locus_tag)

res = res[!duplicated(res$representative),]

# merging with nrtx obj but saving to other obj
nrtxsen = left_join(x = nrtx,
                    y = res,
                    by = "representative") %>% 
  mutate(hfq_int = case_when(is.na(hfq_int) ~ "Não",
                             TRUE ~ hfq_int))

# copying to clipboard 
nrtxsen %>% 
  clipr::write_clip()

