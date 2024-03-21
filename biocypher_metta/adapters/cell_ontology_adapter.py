from biocypher_metta.adapters import Adapter
import gzip
import csv

class CellOntologyAdapter(Adapter):
    """
    Adapter for Cell Ontology dataset
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.label = "cell"
        self.source = "Cell Ontology"
        self.source_url = "https://bioportal.bioontology.org/ontologies/CL" 
        super(CellOntologyAdapter, self).__init__()

    def get_nodes(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.reader(f)
            next(reader)  
            for row in reader:
                cell_id = row[0]
                cell_label = row[1]
                parent_id = row[2] if row[2] != "None" else None
                
                properties = {
                    "label": cell_label,
                    "parent_id": parent_id
                }

                yield cell_id, self.label, properties

    def get_edges(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.reader(f)
            next(reader)  
            for row in reader:
                cell_id = row[0]
                parent_id = row[2] if row[2] != "None" else None

                if parent_id:
                    yield parent_id, cell_id, "is_parent_of", {}
                else:
                    yield "ROOT", cell_id, "is_root_of", {}

