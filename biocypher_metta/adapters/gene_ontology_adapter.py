import rdflib
from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter

class GeneOntologyAdapter(OntologyAdapter):
    ONTOLOGIES = {
        'go': 'http://purl.obolibrary.org/obo/go.owl'
    }

    NAMESPACE = rdflib.term.URIRef('http://www.geneontology.org/formats/oboInOwl#hasOBONamespace')
    EXACT_SYNONYM = rdflib.term.URIRef('http://www.geneontology.org/formats/oboInOwl#hasExactSynonym')
    RELATED_SYNONYM = rdflib.term.URIRef('http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym')

    def __init__(self, write_properties, add_provenance, ontology, type, label='go', dry_run=False):        
        super(GeneOntologyAdapter, self).__init__(write_properties, add_provenance, ontology, type, label, dry_run)

    def get_ontology_source(self, ontology):
        """
        Returns the source and source URL for the Gene Ontology.
        """
        if ontology == 'go':
            return 'Gene Ontology', 'http://purl.obolibrary.org/obo/go.owl'
        else:
            return None, None

    def find_go_nodes(self, graph):
        # subontologies are defined as `namespaces
        nodes_in_namespaces = list(graph.subject_objects(predicate=GeneOntologyAdapter.NAMESPACE))

        node_namespace_lookup = {}
        for n in nodes_in_namespaces:
            node = n[0]
            namespace = n[1]

            node_key = OntologyAdapter.to_key(node)
            node_namespace_lookup[node_key] = str(namespace)
        return node_namespace_lookup

    def get_nodes(self):
        nodes = super().get_nodes()

        # Find GO nodes and their namespaces
        if self.graph is not None:
            nodes_in_go_namespaces = self.find_go_nodes(self.graph)
            for node_id, label, props in nodes:
                props['subontology'] = nodes_in_go_namespaces.get(node_id, None)
                yield node_id, label, props
        else:
            for node_id, label, props in nodes:
                yield node_id, label, props