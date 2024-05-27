import requests
from requests.adapters import HTTPAdapter, Retry
from requests.exceptions import JSONDecodeError
from tqdm import tqdm
import csv
from biocypher_metta.adapters import Adapter


# Example pathway input file:
# R-GGA-199992	trans-Golgi Network Vesicle Budding	Gallus gallus
# R-HSA-164843	2-LTR circle formation	Homo sapiens
# R-HSA-73843	5-Phosphoribose 1-diphosphate biosynthesis	Homo sapiens
# R-HSA-1971475	A tetrasaccharide linker sequence is required for GAG synthesis	Homo sapiens
# R-HSA-5619084	ABC transporter disorders	Homo sapiens


class ReactomePathwayAdapter(Adapter):

    def __init__(self, filepath, pubmed_map_path, write_properties, add_provenance):

        self.filepath = filepath
        self.pubmed_map_path = pubmed_map_path
        self.load_pubmed_map()
        self.label = 'pathway'
        self.dataset = 'pathway'
        self.source = "REACTOME"
        self.source_url = "https://reactome.org"

        super(ReactomePathwayAdapter, self).__init__(write_properties, add_provenance)
    
    def load_pubmed_map(self):
        self.pubmed_map = {}
        with open(self.pubmed_map_path, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                pathway_id, pubmed_id = row[0], row[0]
                self.pubmed_map[pathway_id] = pubmed_id

    def get_nodes(self):
        with open(self.filepath) as input:
            for line in input:
                id, name, species = line.strip().split('\t')
                if species == 'Homo sapiens':
                    props = {}
                    if self.write_properties:
                        props['pathway_name'] = name
                    
                        pubmed_id = self.pubmed_map.get(id, None)
                        if pubmed_id is not None:
                            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{self.pubmed_map[id]}"
                            props['evidence'] = pubmed_url,
                        
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield id, self.label, props