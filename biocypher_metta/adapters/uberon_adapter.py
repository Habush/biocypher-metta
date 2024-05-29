from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter

class UberonAdapter(OntologyAdapter):
    def __init__(self, write_properties, add_provenance, ontology, type, label="ontology term", dry_run=False):
        super().__init__(write_properties, add_provenance, ontology, type, label, dry_run)
        self.source = "UBERON"
        self.source_url = "http://purl.obolibrary.org/obo/uberon.owl"
        self.ONTOLOGIES = {
            'uberon': 'http://purl.obolibrary.org/obo/uberon.owl'
        }
