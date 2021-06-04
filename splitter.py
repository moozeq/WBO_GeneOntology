from utils import read_seqs

if __name__ == '__main__':
    genes = read_seqs('data/genes_e_coli.fasta')
    genes = [gene.id for gene in genes]
    with open('results/population.txt', 'w') as f:
        f.write('\n'.join(genes))
    with open('results/studyA.txt', 'w') as f:
        f.write('\n'.join(gene for gene in genes if gene[-1] == 'A'))
    with open('results/studyB.txt', 'w') as f:
        f.write('\n'.join(gene for gene in genes if gene[-1] == 'B'))
