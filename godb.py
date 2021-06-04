import argparse
import csv
import itertools
import json
from dataclasses import dataclass
from pathlib import Path

from pronto import Ontology

from utils import Term, download_and_unpack

ENRICH_ROW_LEN = 14


@dataclass
class Enrich:
    go: str
    NS: str
    enrichment: str
    name: str
    ratio_in_study: str
    ratio_in_pop: str
    p_uncorrected: float
    depth: int
    study_count: int
    p_bonferroni: float
    p_sidak: float
    p_holm: float
    p_fdr_bh: float
    study_items: list

    def __post_init__(self):
        self.p_uncorrected = float(self.p_uncorrected)
        self.depth = int(self.depth)
        self.study_count = int(self.study_count)
        self.p_bonferroni = float(self.p_bonferroni)
        self.p_sidak = float(self.p_sidak)
        self.p_holm = float(self.p_holm)
        self.p_fdr_bh = float(self.p_fdr_bh)
        self.study_items = list(self.study_items.split(', '))

        i, o = self.ratio_in_study.split('/')
        self.ration_in_study_perc = 100 * int(i) / int(o)
        i, o = self.ratio_in_pop.split('/')
        self.ratio_in_pop_perc = 100 * int(i) / int(o)

    @staticmethod
    def to_tsv_row(enrich_dict: dict):
        return f'{enrich_dict["go"]}\t{enrich_dict["NS"]}\t{enrich_dict["name"]}\t{enrich_dict["ratio_in_study"]}\t{enrich_dict["ratio_in_pop"]}\t{enrich_dict["p_uncorrected"]}\t{enrich_dict["p_bonferroni"]}'


class GO:
    def __init__(self, db_name: str):
        self.db_name = db_name
        self.ontology = Ontology(db_name)

    def term_from_gos(self, gos: list) -> dict:
        terms = {}
        for go in gos:
            term = self.ontology[go]
            terms[go] = {'name': term.name, 'definition': term.definition}
        return terms

    def annotate(self, input_file: str, assocs_file: str, out: str):
        with open(input_file) as f:
            study = json.load(f)
        with open(assocs_file) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                gene, gos = row
                if gene in study:
                    study[gene]['go'] = self.term_from_gos(gos.split(';'))

        with open(out, 'w') as f:
            json.dump(study, f, indent=4)

    @classmethod
    def go_db_from_url(cls,
                       godb_url: str = 'http://purl.obolibrary.org/obo/go.obo'):
        godb = 'godb/go.obo'
        if not Path(godb).exists():
            Term.info(f'GO database not found, downloading...')
            if download_and_unpack(godb_url, godb):
                Term.ok(f'GO database stored at = {godb}')

        return cls(godb)

    @staticmethod
    def analyze_enrichment(input_file: str):
        with open(input_file) as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)
            enrichments = [
                Enrich(*row)
                for row in reader
                if len(row) == ENRICH_ROW_LEN
            ]
            enrichments_cat = {
                term: sorted([e.__dict__ for e in enrichments if e.NS == term], key=lambda e: e['p_uncorrected'])[:5]
                for term in ['BP', 'MF', 'CC']
            }
            return enrichments_cat


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', type=str, choices=['init', 'annotate', 'analyze'],
                        help='GO mode, "init" for databases download, "annotate" for annotation, "analyze" for enrichment analysis')
    args = parser.parse_args()
    if args.mode == 'init':
        go = GO.go_db_from_url()
    elif args.mode == 'annotate':
        go = GO.go_db_from_url()
        go.annotate('results/study_annotated.json', 'results/ecoli.assocs', 'results/study_annotated.json')
    elif args.mode == 'analyze':
        a = GO.analyze_enrichment('results/enrichmentA.tsv')
        b = GO.analyze_enrichment('results/enrichmentB.tsv')
        ab = GO.analyze_enrichment('results/enrichmentAB.tsv')
        pfam = GO.analyze_enrichment('results/enrichmentPfam.tsv')
        results = {
            'a': a,
            'b': b,
            'ab': ab,
            'pfam': pfam
        }
        with open('results/FINAL_ENRICHMENT.json', 'w') as f:
            json.dump(results, f, indent=4)

        for res_type, res in results.items():
            with open(f'results/FINAL_ENRICHMENT_{res_type.upper()}.tsv', 'w') as f:
                all_res = [r for r in res.values()]
                gathered_res = list(itertools.chain(*all_res))
                gathered_res = [Enrich.to_tsv_row(e) for e in gathered_res]
                f.writelines('\n'.join(gathered_res))
    else:
        raise Exception('Wrong mode')

