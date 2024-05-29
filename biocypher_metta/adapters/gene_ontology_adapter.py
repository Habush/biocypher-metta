from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter

class GeneOntologyAdapter(OntologyAdapter):
    def __init__(self, write_properties, add_provenance, ontology, type, label="ontology term", dry_run=False):
        super().__init__(write_properties, add_provenance, ontology, type, label, dry_run)
        self.source = "Gene Ontology"
        self.source_url = "http://purl.obolibrary.org/obo/go.owl"
        self.ontology = ontology
        self.ONTOLOGIES = {
            'go': 'http://purl.obolibrary.org/obo/go.owl'
        }

    def find_go_nodes(self, graph):
        # subontologies are defined as `namespaces
        nodes_in_namespaces = list(graph.subject_objects(predicate=OntologyAdapter.NAMESPACE))

        node_namespace_lookup = {}
        for n in nodes_in_namespaces:
            node = n[0]
            namespace = n[1]

            node_key = OntologyAdapter.to_key(node)
            node_namespace_lookup[node_key] = str(namespace)
        return node_namespace_lookup

    def get_nodes(self):
        nodes = super().get_nodes()
        self.graph = self.get_graph()
        # Find GO nodes and their namespaces
        nodes_in_go_namespaces = self.find_go_nodes(self.graph)
        for node_id, label, props in nodes:
            props['subontology'] = nodes_in_go_namespaces.get(node_id, None)
            yield node_id, label, props

