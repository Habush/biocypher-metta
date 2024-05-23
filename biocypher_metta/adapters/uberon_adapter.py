from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter

class UberonAdapter(OntologyAdapter):
    def __init__(self, write_properties, add_provenance, type, label='ontology term', dry_run=False):
        super().__init__(write_properties, add_provenance, type, label, dry_run)
        self.source = "UBERON"
        self.source_url = "http://purl.obolibrary.org/obo/uberon.owl"
        self.ONTOLOGIES = {
            'uberon': 'http://purl.obolibrary.org/obo/uberon.owl'
        }

    def get_graph(self, ontology='uberon'):
        return super().get_graph(ontology)
