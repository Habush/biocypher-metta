import rdflib
from owlready2 import *
from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter

class CellOntologyAdapter(OntologyAdapter):
    def __init__(self, write_properties, add_provenance, type, label='cell', dry_run=False):
        self.label = label
        self.source = "Cell Ontology"
        self.source_url = "http://purl.obolibrary.org/obo/clo.owl"
        super(CellOntologyAdapter, self).__init__(write_properties, add_provenance, label, type, dry_run)


        self.ONTOLOGIES = {
            'clo': 'http://purl.obolibrary.org/obo/clo.owl'
        }


    def get_graph(self, ontology='clo'):
        if ontology not in self.ONTOLOGIES:
            raise ValueError(f"Ontology '{ontology}' is not defined in this adapter.")
        
        onto = get_ontology(self.ONTOLOGIES[ontology]).load()
        self.graph = default_world.as_rdflib_graph()
        self.clear_cache()
        return self.graph


    def get_nodes(self):
        self.graph = self.get_graph()
        self.cache_node_properties()

        nodes = self.graph.all_nodes()

        i = 0  # dry run is set to true just output the first 1000 nodes
        for node in nodes:
            if i > 100 and self.dry_run:
                break
            if not isinstance(node, rdflib.term.URIRef):
                continue

            term_id = OntologyAdapter.to_key(node)
            term_name = ', '.join(self.get_all_property_values_from_node(node, 'term_names'))
            description = ' '.join(self.get_all_property_values_from_node(node, 'descriptions'))
            synonyms = self.get_all_property_values_from_node(node, 'related_synonyms') + self.get_all_property_values_from_node(node, 'exact_synonyms')

            # Check for problematic characters or formatting
            if '"' in description:
                print(f"Skipping node {term_id} due to problematic characters or formatting.")
                continue

            props = {}
            if self.write_properties:
                props['term_name'] = term_name
                props['description'] = description
                props['synonyms'] = synonyms

                if self.add_provenance:
                    props['source'] = self.source
                    props['source_url'] = self.source_url
            i += 1
            yield term_id, self.label, props
    def get_edges(self):
        self.graph = self.get_graph()
        self.cache_edge_properties()
        for predicate in OntologyAdapter.PREDICATES:
            edges = list(self.graph.subject_objects(predicate=predicate, unique=True))
            i = 0  # dry run is set to true just output the first 100 relationships
            for edge in edges:
                if i > 100 and self.dry_run:
                    break
                from_node, to_node = edge

                if self.is_blank(from_node):
                    continue

                if self.is_blank(to_node) and self.is_a_restriction_block(to_node):
                    restriction_predicate, restriction_node = self.read_restriction_block(to_node)
                    if restriction_predicate is None or restriction_node is None:
                        continue

                    predicate = restriction_predicate
                    to_node = restriction_node

                if self.type == 'edge':
                    from_node_key = OntologyAdapter.to_key(from_node)
                    predicate_key = OntologyAdapter.to_key(predicate)
                    to_node_key = OntologyAdapter.to_key(to_node)

                    if predicate == OntologyAdapter.DB_XREF:
                        if to_node.__class__ == rdflib.term.Literal:
                            if str(to_node) == str(from_node):
                                print('Skipping self xref for: ' + from_node_key)
                                continue

                            if len(str(to_node).split(':')) != 2:
                                print('Unsupported format for xref: ' + str(to_node))
                                continue

                            to_node_key = str(to_node).replace(':', '_')

                            if from_node_key == to_node_key:
                                print('Skipping self xref for: ' + from_node_key)
                                continue
                        else:
                            print('Ignoring non literal xref: {}'.format(str(to_node)))
                            continue

                    predicate_name = self.predicate_name(predicate)
                    if predicate_name == 'dbxref':
                        continue  
                    props = {}
                    if self.write_properties:
                        props['rel_type'] = predicate_name
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield from_node_key, to_node_key, self.label, props
                    i += 1