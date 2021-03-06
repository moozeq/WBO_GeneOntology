#!/usr/bin/env zsh

alias activate=". ./venv/bin/activate"

activate

mkdir results

python3 splitter.py
python3 blaster.py
python3 hammer.py
python3 godb.py init

python3 ./venv/bin/fetch_associations.py --taxon_id 83333 -o results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentPfam.tsv results/studyPfam.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentA.tsv results/studyA.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentB.tsv results/studyB.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentAB.tsv --compare results/studyA.txt results/studyB.txt results/ecoli.assocs

python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentPfam.tsv results/studyPfam.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentA.tsv results/studyA.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentB.tsv results/studyB.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentAB.tsv --compare results/studyA.txt results/studyB.txt results/ecoli.assocs

python3 godb.py annotate
python3 godb.py analyze