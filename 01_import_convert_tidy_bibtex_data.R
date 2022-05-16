# Publication authorship
# 
# Import, convert, and tidy bibtex data downloaded from SCOPUS
#
# Copyright (C) 2022 by Dr. Steffen Blaschke (sbl.msc@cbs.dk).
# This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY).
# https://creativecommons.org/licenses/by/4.0/

# load libraries
library(bibliometrix) # bibliometrics
library(tidyverse) # tidy everything
library(tidytext) # tidy text
library(tm) # text mining

# set working directory
# setwd("YourDirectoryHere")

# tidy bib #
tidy_bib <- function(bib) {
  bib <- bib %>%
  # remove documents that are not articles (e.g., editorials, book reviews)
  filter(DT == "ARTICLE") %>%
  # remove articles from 2020 (i.e., keep only articles from 2010 to 2019)
  filter(PY != 2020) %>%
  # remove double DOIs
  filter(!duplicated(DI, fromLast = TRUE)) %>%
  # remove articles without cited references
  filter(!is.na(CR)) %>%
  # remove articles without abstracts
  filter(!is.na(AB)) %>%
  # remove anonymous authors
  filter(AU != "[ANONYMOUS] A")
  # return bib
  return(bib)
}

# remove copyright notice from abstracts
remove_copyright <- function(bib) {
  bib$AB <- bib$AB %>%
  # generic
  # dot, space, 20, two digits for year, characters, end of string
  str_remove("\\. 20\\d+.*$") %>%
  # space, COPYRIGHT, 20, two digits for year, characters, end of string
  str_remove("\\. COPYRIGHT .*$") %>%
  # dot, space, THE AUTHOR(S) 20, two digits for year, characters, end of string
  str_remove("\\. THE AUTHOR 20\\d+.*$") %>%
  # dot, space, THE AUTHOR(S) 20, two digits for year, characters, end of string
  str_remove("\\. THE AUTHOR\\(S\\) 20\\d+.*$") %>%
  # accounting
  # dot, space, SPRINGER SCIENCE, characters, end of string
  str_remove(" SPRINGER SCIENCE.*$") %>%
  # dot, space, CROWN COPYRIGHT, characters, end of string
  str_remove(" CROWN COPYRIGHT.*$") %>%
  # astronomy
  # dot, space, 2 017. THE AMERICAN ASTRONOMICAL SOCIETY, characters, end of string
  str_remove("&COPY;2010\\. THE AMERICAN ASTRONOMICAL SOCIETY.*$") %>%
  # dot, space, 2 017. THE AMERICAN ASTRONOMICAL SOCIETY, characters, end of string
  str_remove("2 017. THE AMERICAN ASTRONOMICAL SOCIETY.*$") %>%
  # dot, space, YEAR. THE AMERICAN ASTRONOMICAL SOCIETY, characters, end of string  
  str_remove("YEAR\\. THE AMERICAN ASTRONOMICAL SOCIETY.*$") %>%
  # dot, space, THE AMERICAN ASTRONOMICAL SOCIETY, characters, end of string
  str_remove("\\. THE AMERICAN ASTRONOMICAL SOCIETY.*$") %>%
  # dot, space, 20, two digits for year, THE AMERICAN ASTRONOMICAL SOCIETY, characters, end of string
  str_remove("20\\d+ THE AMERICAN ASTRONOMICAL SOCIETY.*$") %>%
  # dot, space, ASTRRONOMICAL SOCIETY OF AUSTRALIA, characters, end of string
  str_remove("\\. ASTRONOMICAL SOCIETY OF AUSTRALIA.*$") %>%
  # dot, space, JOURNAL COMPILATION AUSTRALIAN MAMMAL SOCIETY, characters, end of string
  str_remove("\\. JOURNAL COMPILATION AUSTRALIAN MAMMAL SOCIETY.*$") %>%
  # gastroenterology
  # dot, space, SPRINGER, characters, end of string
  str_remove("\\. SPRINGER.*$") %>%
  # space, SPRINGER 20, two digits for year, characters, end of string
  str_remove(" SPRINGER 20\\d+.*$") %>%
  # 20, two digits for yeaar, space, SPRINGER, characters, end of string
  str_remove("20\\d+ SPRINGER.*$")
  # return bib
  return(bib)
}

# tidy abstract
tidy_abstract <- function(bib) {
  bib$AB <- bib$AB %>%
  # lower case
  tolower() %>%
  # remove punctuation
  removePunctuation(preserve_intra_word_contractions = TRUE, preserve_intra_word_dashes = TRUE) %>%
  # remove stopwords
  removeWords(stopwords(kind = "en")) %>%
  # remove numbers
  removeNumbers() %>%
  # strip whitespace
  stripWhitespace() %>%
  # stem
  stemDocument()
  # return bib
  return(bib)
}

# import bibtex data
# accounting
acc.bib.raw <- convert2df("scopus_accounting_top_10_2010-2019.bib", dbsource = "scopus", format = "bibtex")
# astronomy
ast.bib.raw <- convert2df("scopus_astronomy_and_astrophysics_top_10_2010-2019.bib", dbsource = "scopus", format = "bibtex")
# gastroenterology
gas.bib.raw <- convert2df("scopus_gastroenterology_top_10_2010-2019.bib", dbsource = "scopus", format = "bibtex") %>%
  # remove "JOURNAL OF GASTRIC CANCER" (which is not on top 10 list)
  filter(SO != "JOURNAL OF GASTRIC CANCER")

# tidy bib, remove copyright notice, tidy abstract
acc.bib <- acc.bib.raw %>%
  tidy_bib() %>%
  remove_copyright() %>%
  tidy_abstract()

ast.bib <- ast.bib.raw %>%
  tidy_bib() %>%
  remove_copyright() %>%
  tidy_abstract()

gas.bib <- gas.bib.raw %>%
  tidy_bib() %>%
  remove_copyright() %>%
  tidy_abstract()

# clean up
detach("package:tm", unload = TRUE)
rm(tidy_bib, remove_copyright, tidy_abstract)
