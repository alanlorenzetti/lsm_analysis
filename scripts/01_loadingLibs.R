# alorenzetti 20210211

# description ####
# this script will load the libs
# necessary to run this collection
# of scripts

# loading pacman
if(!require(pacman)){install.packages("pacman"); library(pacman)}

# packages
packs = c(
  "tidyverse",
  "rtracklayer",
  "GenomicRanges",
  "BSgenome",
  "Biostrings",
  "ggpubr",
  "ggbeeswarm",
  "gtools",
  "eulerr",
  "ggthemes",
  "tm",
  "ggwordcloud"
)

# loading
p_load(char = packs)

# setting ggplot2 theme
theme_set(theme_bw())
