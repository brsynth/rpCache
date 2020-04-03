#!/bin/sh

# TERM=xterm
# bold=$(tput bold)
# normal=$(tput sgr0)

# # input_cache
# for file in chem_xref.tsv reac_xref.tsv chem_prop.tsv comp_xref.tsv
# do
#   echo ${bold}Downloading $file...${normal}
#   wget --quiet https://www.metanetx.org/cgi-bin/mnxget/mnxref/$file -P /home/input_cache/
# done
#
# echo ${bold}Downloading rules_rall_rp3.tar.gz...${normal}
# wget --quiet https://retrorules.org/dl/preparsed/rr02/rp3/hs -O /home/rules_rall_rp3.tar.gz
# tar xf /home/rules_rall_rp3.tar.gz -C /home/
# mv /home/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv /home/input_cache/rules_rall.tsv
#
# echo ${bold}Downloading rr02_more_data.tar.gz...${normal}
# wget --quiet https://retrorules.org/dl/this/is/not/a/secret/path/rr02 -O /home/rr02_more_data.tar.gz
# tar xf /home/rr02_more_data.tar.gz -C /home/
# mv /home/rr02_more_data/compounds.tsv /home/input_cache/rr_compounds.tsv
# mv /home/rr02_more_data/rxn_recipes.tsv /home/input_cache/

echo ${bold}Generating pickles...${normal}
python3 /home/rpCache.py
