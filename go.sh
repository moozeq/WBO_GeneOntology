#!/usr/bin/env zsh

alias activate=". ../venv/bin/activate"

activate

python3 splitter.py
python3 blaster.py
python3 hammer.py
python3 godb.py init

python3 ../venv/bin/fetch_associations.py --taxon_id 83333 -o data/ecoli.assocs
python3 ../venv/bin/find_enrichment.py --pval=0.05 --indent --min_overlap 0.5 --obo godb/go.obo --outfile data/enrichmentPfam.tsv data/studyPfam.txt data/population.txt data/ecoli.assocs
python3 ../venv/bin/find_enrichment.py --pval=0.05 --indent --min_overlap 0.5 --obo godb/go.obo --outfile data/enrichmentA.tsv data/studyA.txt data/population.txt data/ecoli.assocs
python3 ../venv/bin/find_enrichment.py --pval=0.05 --indent --min_overlap 0.5 --obo godb/go.obo --outfile data/enrichmentB.tsv data/studyB.txt data/population.txt data/ecoli.assocs

python3 godb.py annotate