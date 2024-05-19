import os
import json
import rdflib
from owlready2 import *
from biocypher_metta.adapters import Adapter

class OntologyAdapter(Adapter):
    HAS_PART = rdflib.term.URIRef('http://purl.obolibrary.org/obo/BFO_0000051')
    PART_OF = rdflib.term.URIRef('http://purl.obolibrary.org/obo/BFO_0000050')
    SUBCLASS = rdflib.term.URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')
    DB_XREF = rdflib.term.URIRef('http://www.geneontology.org/formats/oboInOwl#hasDbXref')

    LABEL = rdflib.term.URIRef('http://www.w3.org/2000/01/rdf-schema#label')
    RESTRICTION = rdflib.term.URIRef('http://www.w3.org/2002/07/owl#Restriction')
    TYPE = rdflib.term.URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
    ON_PROPERTY = rdflib.term.URIRef('http://www.w3.org/2002/07/owl#onProperty')
    SOME_VALUES_FROM = rdflib.term.URIRef('http://www.w3.org/2002/07/owl#someValuesFrom')
    ALL_VALUES_FROM = rdflib.term.URIRef('http://www.w3.org/2002/07/owl#allValuesFrom')
    EXACT_SYNONYM = rdflib.term.URIRef('http://www.geneontology.org/formats/oboInOwl#hasExactSynonym')
    RELATED_SYNONYM = rdflib.term.URIRef('http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym')
    DESCRIPTION = rdflib.term.URIRef('http://purl.obolibrary.org/obo/IAO_0000115')

    PREDICATES = [SUBCLASS, DB_XREF]
    RESTRICTION_PREDICATES = [HAS_PART, PART_OF]

    def __init__(self, write_properties, add_provenance, label, type='node', dry_run=False, subontologies=None):
        super().__init__(write_properties, add_provenance)
        self.type = type
        self.dry_run = dry_run
        self.subontologies = subontologies if subontologies else []
        self.label = label
        self.ONTOLOGIES = {}

    def get_graph(self, ontology):
        if ontology not in self.ONTOLOGIES:
            raise ValueError(f"Ontology '{ontology}' is not defined in this adapter.")
        
        onto = get_ontology(self.ONTOLOGIES[ontology]).load()
        self.graph = default_world.as_rdflib_graph()
        self.clear_cache()
        return self.graph

    def get_nodes(self):
        for ontology in self.ONTOLOGIES.keys():
            self.graph = self.get_graph(ontology)
            self.cache_node_properties()

            nodes = list(self.graph.subjects())
            if self.subontologies:
                nodes = self.filter_nodes_by_subontology(nodes, self.graph)

            i = 0  # dry run is set to true just output the first 1000 nodes
            for node in nodes:
                if i > 100 and self.dry_run:
                    break
                if not isinstance(node, rdflib.term.URIRef):
                    continue

                term_id = OntologyAdapter.to_key(node)
                props = {
                    'term_name': ', '.join(self.get_all_property_values_from_node(node, 'term_names')),
                    'description': ' '.join(self.get_all_property_values_from_node(node, 'descriptions')),
                    'synonyms': self.get_all_property_values_from_node(node, 'related_synonyms') +
                                self.get_all_property_values_from_node(node, 'exact_synonyms'),
                }
                if self.subontologies:
                    props['subontology'] = self.get_subontology(node, self.graph)
                
                i += 1
                yield term_id, self.label, props

    def get_edges(self):
        for ontology in self.ONTOLOGIES.keys():
            self.graph = self.get_graph(ontology)
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
                            continue  # TODO should we skip dbxref edges?
                        props = {
                            'rel_type': self.predicate_name(predicate)
                        }

                        yield from_node_key, to_node_key, self.label, props
                        i += 1

    def predicate_name(self, predicate):
        predicate = str(predicate)
        if predicate == str(OntologyAdapter.HAS_PART):
            return 'has_part'
        elif predicate == str(OntologyAdapter.PART_OF):
            return 'part_of'
        elif predicate == str(OntologyAdapter.SUBCLASS):
            return 'subclass'
        elif predicate == str(OntologyAdapter.DB_XREF):
            return 'dbxref'
        return ''

    @classmethod
    def to_key(cls, node_uri):
        key = str(node_uri).split('/')[-1]
        key = key.replace('#', '.').replace('?', '_')
        key = key.replace('&', '.').replace('=', '_')
        key = key.replace('/', '_').replace('~', '.')
        key = key.replace('_', ':')
        key = key.replace(' ', '')

        if key.replace('.', '').isnumeric():
            key = '{}_{}'.format('number', key)

        return key

    def is_a_restriction_block(self, node):
        node_type = self.get_all_property_values_from_node(node, 'node_types')
        return node_type and node_type[0] == OntologyAdapter.RESTRICTION

    def read_restriction_block(self, node):
        restricted_property = self.get_all_property_values_from_node(node, 'on_property')

        if restricted_property and restricted_property[0] not in OntologyAdapter.RESTRICTION_PREDICATES:
            return None, None

        restriction_predicate = str(restricted_property[0])

        some_values_from = self.get_all_property_values_from_node(node, 'some_values_from')
        if some_values_from:
            return (restriction_predicate, some_values_from[0])

        all_values_from = self.get_all_property_values_from_node(node, 'all_values_from')
        if all_values_from:
            return (restriction_predicate, all_values_from[0])

        return (None, None)

    def filter_nodes_by_subontology(self, nodes, graph):
        node_namespace_lookup = {node: self.get_subontology(node, graph) for node in nodes}
        return {node for node in nodes if node_namespace_lookup[node] in self.subontologies}

    def get_subontology(self, node, graph):
        namespace_predicates = list(graph.objects(subject=node, predicate=OntologyAdapter.NAMESPACE))
        if namespace_predicates:
            return str(namespace_predicates[0])
        return None

    def is_blank(self, node):
        BLANK_NODE = rdflib.term.BNode
        return isinstance(node, BLANK_NODE)

    def clear_cache(self):
        self.cache = {}

    def cache_edge_properties(self):
        self.cache['node_types'] = self.cache_predicate(OntologyAdapter.TYPE)
        self.cache['on_property'] = self.cache_predicate(OntologyAdapter.ON_PROPERTY)
        self.cache['some_values_from'] = self.cache_predicate(OntologyAdapter.SOME_VALUES_FROM)
        self.cache['all_values_from'] = self.cache_predicate(OntologyAdapter.ALL_VALUES_FROM)

    def cache_node_properties(self):
        self.cache['term_names'] = self.cache_predicate(OntologyAdapter.LABEL)
        self.cache['descriptions'] = self.cache_predicate(OntologyAdapter.DESCRIPTION)
        self.cache['node_types'] = self.cache_predicate(OntologyAdapter.TYPE)
        self.cache['related_synonyms'] = self.cache_predicate(OntologyAdapter.RELATED_SYNONYM)
        self.cache['exact_synonyms'] = self.cache_predicate(OntologyAdapter.EXACT_SYNONYM)

    def cache_predicate(self, predicate):
        return {s: str(o) for s, o in self.graph.subject_objects(predicate=predicate)}

    def get_all_property_values_from_node(self, node, prop):
        if node not in self.cache[prop]:
            return []

        values = self.cache[prop][node]
        if not isinstance(values, list):
            values = [values]
        return values