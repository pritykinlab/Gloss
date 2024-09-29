## basically a class that holds the pathway acquisition logic
## im gonna supply the pathways myself, but will give the options to refresh kegg
import pickle

class PathwaysModule():
    def __init__(self, pathway_string, adata_genes):
        self.adata_genes = adata_genes
        self.pathways = {}
        self.accepted_pathways = set('hallmark', 'kegg')
        if pathway_string not in self.accepted_pathways:
            raise ValueError('pathways string not in {}'.format(", ".join(self.accepted_pathways)))
        elif pathway_string == 'hallmark':
            self.pathways = self.get_hallmark_pathways()
        elif pathway_string == 'kegg':
            self.pathways = self.get_kegg_pathways()

    def get_hallmark_pathways(self):
        pathways = {}
        hallmark_pathways_path = './data/mouse_hallmark_genes.gmt.txt'
        with open(hallmark_pathways_path, "r") as f:
            for line in f:
                line = line.strip()
                words = line.split("\t")
                pathway = words[2:]
                pathway = list(set(pathway).intersection(self.adata_genes))
                pathways[words[0]] = pathway
        return pathways

    def get_kegg_pathways(self):
        with open('./data/kegg_pathways_sep18_2024.pkl', 'rb') as handle:
            kegg_gene_sets = pickle.load(handle)
        for key, value in kegg_gene_sets.items():
            kegg_gene_sets[key] = list(set([gene for gene in value if gene in self.adata_genes]))
            
        empty = []
        for key in kegg_gene_sets:
            if not kegg_gene_sets[key]:
                empty.append(key)
        for pathway in empty:
            kegg_gene_sets.pop(pathway)
            
        return kegg_gene_sets