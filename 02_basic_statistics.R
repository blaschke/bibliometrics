# Publication authorship
# 
# Basic statistics (articles per journal)
#
# Copyright (C) 2022 by Dr. Steffen Blaschke (sbl.msc@cbs.dk).
# This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY).
# https://creativecommons.org/licenses/by/4.0/

# articles per journal
# accounting
xtable::xtable(
  acc.bib %>%
  count(SO),
  type = "latex")
# astronomy
xtable::xtable(
  ast.bib %>%
    count(SO),
  type = "latex")
# gastroenterology
xtable::xtable(
  gas.bib %>%
    count(SO),
  type = "latex")

# check number of authors versus number of references
acc.bib %>%
  filter(length(str_split(AU, ";", simplify = TRUE)) > length(str_split(CR, ";", simplify = TRUE)))
ast.bib %>%
  filter(length(str_split(AU, ";", simplify = TRUE)) > length(str_split(CR, ",", simplify = TRUE))) %>%
acc.bib %>%
  filter(length(str_split(AU, ";", simplify = TRUE)) > length(str_split(CR, ";", simplify = TRUE)))