from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter
class CellOntologyAdapter(OntologyAdapter):
    def __init__(self, write_properties, add_provenance, ontology, type, label="cell", dry_run=False):
        super().__init__(write_properties, add_provenance, ontology, type, label, dry_run)
        self.source = "Cell Ontology"
        self.source_url = "http://purl.obolibrary.org/obo/clo.owl"
        self.ONTOLOGIES = {
            'clo': 'http://purl.obolibrary.org/obo/clo.owl'
        }
