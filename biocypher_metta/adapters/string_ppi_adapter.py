# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher_metta.adapters import Adapter
import pickle
import csv
import gzip

# Imports STRING Protein-Protein interactions

# protein1 protein2 combined_score
# 9606.ENSP00000000233 9606.ENSP00000356607 173
# 9606.ENSP00000000233 9606.ENSP00000427567 154
# 9606.ENSP00000000233 9606.ENSP00000253413 151
# 9606.ENSP00000000233 9606.ENSP00000493357 471
# 9606.ENSP00000000233 9606.ENSP00000324127 201
# 9606.ENSP00000000233 9606.ENSP00000325266 180
# 9606.ENSP00000000233 9606.ENSP00000320935 181

class StringPPIAdapter(Adapter):
    def __init__(self, filepath, ensembl_to_uniprot_map):
        """
        Constructs StringPPI adapter that returns edges between proteins
        :param filepath: Path to the TSV file downloaded from String
        :param ensembl_to_uniprot_map: file containing pickled dictionary mapping Ensemble Protein IDs to Uniprot IDs
        """
        self.filepath = filepath

        with open(ensembl_to_uniprot_map, "rb") as f:
            self.ensembl2uniprot = pickle.load(f)

        self.label = "interacts_with"
        self.source = "STRING"
        self.source_url = "https://string-db.org/"
        self.version = "v12.0"

    def get_edges(self):
        with gzip.open(self.filepath, "rt") as fp:
            table = csv.reader(fp, delimiter=" ", quotechar='"')
            table.__next__() # skip header
            for row in table:
                protein1 = row[0].split(".")[1]
                protein2 = row[1].split(".")[1]
                if protein1 in self.ensembl2uniprot and protein2 in self.ensembl2uniprot:
                    protein1_uniprot = self.ensembl2uniprot[protein1]
                    protein2_uniprot = self.ensembl2uniprot[protein2]
                    _id = protein1_uniprot + "_" + protein2_uniprot + "_" + self.label
                    _source = protein1_uniprot
                    _target = protein2_uniprot
                    _props = {
                        "score": float(row[2]) / 1000, # divide by 1000 to normalize score
                    }
                    yield _id, _source, _target, self.label, _props