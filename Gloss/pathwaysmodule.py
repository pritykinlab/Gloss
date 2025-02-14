## basically a class that holds the pathway acquisition logic
## im gonna supply the pathways myself, but will give the options to refresh kegg
import pickle
import os 

class PathwaysModule():
    def __init__(self, pathway_string, adata_genes):
        self.package_dir = os.path.dirname(os.path.abspath(__file__))
        self.adata_genes = adata_genes
        self.pathways = {}
        self.accepted_pathways = set(['hallmark', 'kegg'])
        if type(pathway_string) == str:
            if pathway_string not in self.accepted_pathways:
                raise ValueError('pathways string not in {}'.format(", ".join(self.accepted_pathways)))
            elif pathway_string == 'hallmark':
                self.pathways = self._get_hallmark_pathways()
            elif pathway_string == 'kegg':
                self.pathways = self._get_kegg_pathways()
        elif type(pathway_string) == dict:
            self.pathways = self._filter_pathways(pathway_string)

        self.pathway_genes = self._get_pathway_genes()

    def _get_hallmark_pathways(self):
        pathways = {}
        hallmark_pathways_path = './data/mouse_hallmark_genes.gmt.txt'
        data_file_path = os.path.join(self.package_dir, hallmark_pathways_path)
        with open(data_file_path, "r") as f:
            for line in f:
                line = line.strip()
                words = line.split("\t")
                pathway = words[2:]
                pathway = list(set(pathway).intersection(self.adata_genes))
                pathways[words[0]] = pathway

        return self._non_empty_pathways(pathways)

    def _get_kegg_pathways(self):
        kegg_pathways_path = './data/kegg_pathways_sep18_2024.pkl'
        data_file_path = os.path.join(self.package_dir, kegg_pathways_path)
        with open(data_file_path, 'rb') as handle:
            kegg_gene_sets = pickle.load(handle)
        for key, value in kegg_gene_sets.items():
            kegg_gene_sets[key] = list(set([gene for gene in value if gene in self.adata_genes]))
            
        return self._non_empty_pathways(kegg_gene_sets)
    
    def _non_empty_pathways(self, pathway_gene_sets):
        empty = []
        for key in pathway_gene_sets:
            if not pathway_gene_sets[key]:
                empty.append(key)
        for pathway in empty:
            pathway_gene_sets.pop(pathway)
            
        return pathway_gene_sets
    
    def _get_pathway_genes(self):
        return set([ x for mylist in self.pathways.values() for x in mylist ])
    
    def _filter_pathways(self, pathway_gene_sets):
        for key, value in pathway_gene_sets.items():
            pathway_gene_sets[key] = list(set([gene for gene in value if gene in self.adata_genes]))
        return self._non_empty_pathways(pathway_gene_sets)