# BAR Superfamily HMM Search Pipeline

This repository contains scripts for large-scale detection of BAR superfamily
domains using HMMER.

## Requirements
- HMMER >= 3.3
- seqkit
- Python >= 3.9
- pandas

## Included Pfam models
- PF02147 – BAR
- PF06957 – F-BAR
- PF09104 – I-BAR
- PF03008 – DRF
- PF08515 – BAR_IMC
- PF17641 – BAR_dom_N
- PF17642 – BAR_dom_C

## Source
Pfam-A (current release at build time: Jan 2026)


hmm_search_bar.sh - compiles hmms for BAR domain and parses result to csv
hmm_search_bar2.sh - compiles hmms for BAR domain and parses result to csv - fixed bugs version 
domtblout_parser.py - parses domtblout files to a single csv file


