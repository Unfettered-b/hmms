# HMM Search Pipeline

This repository contains scripts for large-scale detection of protein superfamily
domains using HMMER.

## Requirements
- HMMER >= 3.3
- seqkit
- Python >= 3.9
- pandas

## Included Pfam models - BAR
- PF02147 – BAR
- PF06957 – F-BAR
- PF09104 – I-BAR
- PF03008 – DRF
- PF08515 – BAR_IMC
- PF17641 – BAR_dom_N
- PF17642 – BAR_dom_C


## ESCRT AsCOGs
- https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/asgard20/ 
- Required files:
    - asCOG.ali.tgz
    - asCOG.2020-10.def.tab

## Source
Pfam-A (current release at build time: Jan 2026)


hmm_search_bar.sh - compiles hmms for BAR domain and parses result to csv


hmm_search_bar2.sh - compiles hmms for BAR domain and parses result to csv - fixed bugs version


domtblout_parser.py - parses domtblout files to a single csv file


hmm_search_escrt.sh - compiles hmm for ESCRT proteins from AsCOGs alignments


get_ascog_ids.py - returns asCOG ids for given search term for the protein


