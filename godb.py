import argparse
import csv
import json
from pathlib import Path

from pronto import Ontology

from utils import Term, download_and_unpack


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
                       tax_id: int,
                       godb_url: str = 'http://purl.obolibrary.org/obo/go.obo',
                       gafdb_url: str = 'http://current.geneontology.org/annotations/ecocyc.gaf.gz'):
        godb = 'godb/go.obo'
        ecocycdb = 'godb/ecocyc.gaf'
        if not Path(godb).exists():
            Term.info(f'GO database not found, downloading...')
            if download_and_unpack(godb_url, godb):
                Term.ok(f'GO database stored at = {godb}')
        if not Path(ecocycdb).exists():
            Term.info(f'GAF database not found, downloading...')
            if download_and_unpack(gafdb_url, ecocycdb):
                Term.ok(f'GAF database stored at = {ecocycdb}')

        return cls(godb)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', type=str, choices=['init', 'annotate'],
                        help='GO mode, "init" for databases download, "annotate" for analysis')
    args = parser.parse_args()
    if args.mode == 'init':
        go = GO.go_db_from_url(83333)
    elif args.mode == 'annotate':
        go = GO.go_db_from_url(83333)
        go.annotate('data/study_annotated.json', 'data/ecoli.assocs', 'data/study_annotated.json')
    else:
        raise Exception('Wrong mode')

