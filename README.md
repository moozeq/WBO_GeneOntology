# Analysis

_E. coli_ proteins fragments and ontology analysis. You can read report in Polish [here](doc/report_PL.md).

## TL;JustRun

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

./go.sh
```

## Searching the most similar sequences

First step in the pipeline is to search [genes](data/genes_e_coli.fasta) based on [protein fragments](data/protein_fragments.fasta).

```bash
python3 blaster.py
```

Command above builds `BLAST` database and then search it using `hmmscan`.

## Searching PFAM domains

Next step in the pipeline searches for PFAM domains in found genes.

```bash
python3 hammer.py
```

It downloads `.hmm` database and based on it builds indexes which are used for `HMMER` software.

## Gene Ontology annotate

Then script uses results for annotating filtered genes.

```bash
python3 godb.py init
python3 ../venv/bin/fetch_associations.py --taxon_id 83333 -o results/ecoli.assocs
python3 godb.py annotate
```

It uses _E. coli_ `.gaf` file for finding proper `GO` IDs and then uses `go.obo` database to annotate
obtained genes. It also uses `goatools` to create `results/ecoli.assocs` file with easy mapping `gene -> GO IDs` for
_E. coli_.

## Overrepresented analysis

The last step is separated from others in pipeline and uses `goatools` software for analysis overrepresented terms in _E. coli_ genes from `A` and `B` groups.

```bash
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentPfam.tsv results/studyPfam.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentA.tsv results/studyA.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentB.tsv results/studyB.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --min_overlap 0.5 --obo godb/go.obo --outfile results/enrichmentAB.tsv --compare results/studyA.txt results/studyB.txt results/ecoli.assocs

python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentPfam.tsv results/studyPfam.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentA.tsv results/studyA.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentB.tsv results/studyB.txt results/population.txt results/ecoli.assocs
python3 ./venv/bin/find_enrichment.py --method bonferroni --pval 0.5 --min_overlap 0.5 --obo godb/go.obo --outfile results/bon_enrichmentAB.tsv --compare results/studyA.txt results/studyB.txt results/ecoli.assocs
python3 godb.py analyze
```