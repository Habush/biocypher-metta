import os
import json
import rdflib
from owlready2 import *
from biocypher_metta.adapters import Adapter


class OntologyAdapter(Adapter):

    ONTOLOGIES = {
        # 'uberon': 'http://purl.obolibrary.org/obo/uberon.owl',
        # 'clo': 'http://purl.obolibrary.org/obo/clo.owl',
        # 'cl': 'http://purl.obolibrary.org/obo/cl.owl',
        # 'hpo': 'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-01-27/hp.owl',
        # 'mondo': 'https://github.com/monarch-initiative/mondo/releases/download/v2023-02-06/mondo.owl',
        'go': 'http://purl.obolibrary.org/obo/go.owl',
        # 'efo': 'https://github.com/EBISPOT/efo/releases/download/current/efo.owl',
        # 'chebi': 'https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl',
        # 'vario': 'http://www.variationontology.org/vario_download/vario.owl.copy',
        # 'orphanet': 'https://www.orphadata.com/data/ontologies/ordo/last_version/ORDO_en_4.3.owl',
        # 'ncit': 'http://purl.obolibrary.org/obo/ncit.owl'
    }

    GO_SUBONTOLGIES = ['molecular_function',
                       'cellular_component', 'biological_process']

    HAS_PART = rdflib.term.URIRef('http://purl.obolibrary.org/obo/BFO_0000051')
    PART_OF = rdflib.term.URIRef('http://purl.obolibrary.org/obo/BFO_0000050')
    SUBCLASS = rdflib.term.URIRef(
        'http://www.w3.org/2000/01/rdf-schema#subClassOf')
    DB_XREF = rdflib.term.URIRef(
        'http://www.geneontology.org/formats/oboInOwl#hasDbXref')

    LABEL = rdflib.term.URIRef('http://www.w3.org/2000/01/rdf-schema#label')
    RESTRICTION = rdflib.term.URIRef(
        'http://www.w3.org/2002/07/owl#Restriction')
    TYPE = rdflib.term.URIRef(
        'http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
    ON_PROPERTY = rdflib.term.URIRef(
        'http://www.w3.org/2002/07/owl#onProperty')
    SOME_VALUES_FROM = rdflib.term.URIRef(
        'http://www.w3.org/2002/07/owl#someValuesFrom')
    ALL_VALUES_FROM = rdflib.term.URIRef(
        'http://www.w3.org/2002/07/owl#allValuesFrom')
    NAMESPACE = rdflib.term.URIRef(
        'http://www.geneontology.org/formats/oboInOwl#hasOBONamespace')
    EXACT_SYNONYM = rdflib.term.URIRef(
        'http://www.geneontology.org/formats/oboInOwl#hasExactSynonym')
    RELATED_SYNONYM = rdflib.term.URIRef(
        'http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym')
    DESCRIPTION = rdflib.term.URIRef(
        'http://purl.obolibrary.org/obo/IAO_0000115')

    PREDICATES = [SUBCLASS, DB_XREF]
    RESTRICTION_PREDICATES = [HAS_PART, PART_OF]

    def __init__(self, type='node', dry_run=False):
        self.type = type
        self.dry_run = dry_run
        if type == 'node':
            self.label = 'ontology_term'
        elif type == 'edge':
            self.label = 'ontology_relationship'
        else:
            raise ValueError('Invalid type. Allowed values: node, edge')

    def get_graph(self, ontology):

        onto = get_ontology(OntologyAdapter.ONTOLOGIES[ontology]).load()
        self.graph = default_world.as_rdflib_graph()
        # self.graph = onto.as_rdflib_graph()
        self.clear_cache()
        return self.graph

    def get_nodes(self):

        for ontology in OntologyAdapter.ONTOLOGIES.keys():
            self.graph = self.get_graph(ontology)
            self.cache_node_properties()

            nodes_in_go_namespaces = self.find_go_nodes(self.graph)
            nodes = nodes_in_go_namespaces.keys()

            i = 0 # dry run is set to true just output the first 1000 nodes
            for node in nodes:
                if i > 100 and self.dry_run:
                    break
                # avoiding blank nodes and other arbitrary node types
                if not isinstance(node, rdflib.term.URIRef):
                    continue

                # term_id = str(node).split('/')[-1]
                term_id = OntologyAdapter.to_key(node)
                props = {
                    # 'uri': str(node),
                    'term_name': ', '.join(self.get_all_property_values_from_node(node, 'term_names')),
                    'description': ' '.join(self.get_all_property_values_from_node(node, 'descriptions')),
                    'synonyms': self.get_all_property_values_from_node(node, 'related_synonyms') +
                    self.get_all_property_values_from_node(node, 'exact_synonyms'),
                    'subontology': nodes_in_go_namespaces.get(node, None),
                    'source': ontology.upper(),
                    'source_url': OntologyAdapter.ONTOLOGIES[ontology],

                }
                i += 1
                yield term_id, self.label, props

    def get_edges(self):
        for ontology in OntologyAdapter.ONTOLOGIES.keys():
            self.graph = self.get_graph(ontology)
            self.cache_edge_properties()
            for predicate in OntologyAdapter.PREDICATES:
                edges = list(self.graph.subject_objects(
                    predicate=predicate, unique=True))
                i = 0  # dry run is set to true just output the first 100 relationships
                for edge in edges:
                    if i > 100 and self.dry_run:
                        break
                    from_node, to_node = edge

                    if self.is_blank(from_node):
                        continue

                    if self.is_blank(to_node) and self.is_a_restriction_block(to_node):
                        restriction_predicate, restriction_node = self.read_restriction_block(
                            to_node)
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

                                # only accepting IDs in the form <ontology>:<ontology_id>
                                if len(str(to_node).split(':')) != 2:
                                    print(
                                        'Unsupported format for xref: ' + str(to_node))
                                    continue

                                to_node_key = str(to_node).replace(':', '_')

                                if from_node_key == to_node_key:
                                    print('Skipping self xref for: ' + from_node_key)
                                    continue
                            else:
                                print('Ignoring non literal xref: {}'.format(str(to_node)))
                                continue

                        predicate_name = self.predicate_name(predicate)
                        if predicate_name == 'dbxref': continue #TODO should we skip dbxref edges?
                        key = '{}_{}_{}'.format(
                            from_node_key,
                            predicate_key,
                            to_node_key
                        )
                        props = {
                            'rel_type': self.predicate_name(predicate),
                            'source': ontology.upper(),
                            'source_url': OntologyAdapter.ONTOLOGIES[ontology],
                        }

                        yield key, from_node_key, to_node_key, self.label, props
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

    # "http://purl.obolibrary.org/obo/CLO_0027762#subclass?id=123" => "CLO_0027762.subclass_id=123"
    # "12345" => "number_12345" - there are cases where URIs are just numbers, e.g. HPO

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

    # Example of a restriction block:
    # <rdfs:subClassOf>
    #     <owl:Restriction>
    #         <owl:onProperty rdf:resource="http://purl.obolibrary.org/obo/RO_0001000"/>
    #         <owl:someValuesFrom rdf:resource="http://purl.obolibrary.org/obo/CL_0000056"/>
    #     </owl:Restriction>
    # </rdfs:subClassOf>
    # This block must be interpreted as the triple (s, p, o):
    # (parent object, http://purl.obolibrary.org/obo/RO_0001000, http://purl.obolibrary.org/obo/CL_0000056)

    def is_a_restriction_block(self, node):
        node_type = self.get_all_property_values_from_node(node, 'node_types')
        return node_type and node_type[0] == OntologyAdapter.RESTRICTION

    def read_restriction_block(self, node):
        restricted_property = self.get_all_property_values_from_node(
            node, 'on_property')

        # assuming a restriction block will always contain only one `owl:onProperty` triple
        if restricted_property and restricted_property[0] not in OntologyAdapter.RESTRICTION_PREDICATES:
            return None, None

        restriction_predicate = str(restricted_property[0])

        # returning the pair (owl:onProperty value, owl:someValuesFrom or owl:allValuesFrom value)
        # assuming a owl:Restriction block in a rdf:subClassOf will contain only one `owl:someValuesFrom` or `owl:allValuesFrom` triple
        some_values_from = self.get_all_property_values_from_node(
            node, 'some_values_from')
        if some_values_from:
            return (restriction_predicate, some_values_from[0])

        all_values_from = self.get_all_property_values_from_node(
            node, 'all_values_from')
        if all_values_from:
            return (restriction_predicate, all_values_from[0])

        return (None, None)

    def find_go_nodes(self, graph):
        # subontologies are defined as `namespaces`
        nodes_in_namespaces = list(graph.subject_objects(
            predicate=OntologyAdapter.NAMESPACE))

        node_namespace_lookup = {}
        for n in nodes_in_namespaces:
            node = n[0]
            namespace = str(n[1])
            if namespace in OntologyAdapter.GO_SUBONTOLGIES:
                node_namespace_lookup[node] = namespace

        return node_namespace_lookup

    def is_blank(self, node):
        # a BNode according to rdflib is a general node (as a 'catch all' node) that doesn't have any type such as Class, Literal, etc.
        BLANK_NODE = rdflib.term.BNode

        return isinstance(node, BLANK_NODE)


    # it's faster to load all subject/objects beforehand
    def clear_cache(self):
        self.cache = {}

    def cache_edge_properties(self):
        self.cache['node_types'] = self.cache_predicate(OntologyAdapter.TYPE)
        self.cache['on_property'] = self.cache_predicate(OntologyAdapter.ON_PROPERTY)
        self.cache['some_values_from'] = self.cache_predicate(
            OntologyAdapter.SOME_VALUES_FROM)
        self.cache['all_values_from'] = self.cache_predicate(
            OntologyAdapter.ALL_VALUES_FROM)

    def cache_node_properties(self):
        self.cache['term_names'] = self.cache_predicate(OntologyAdapter.LABEL)
        self.cache['descriptions'] = self.cache_predicate(OntologyAdapter.DESCRIPTION)
        self.cache['related_synonyms'] = self.cache_predicate(
            OntologyAdapter.EXACT_SYNONYM)
        self.cache['exact_synonyms'] = self.cache_predicate(
            OntologyAdapter.RELATED_SYNONYM)

    def cache_predicate(self, predicate):
        return list(self.graph.subject_objects(predicate=predicate))

    def get_all_property_values_from_node(self, node, collection):
        values = []
        for subject_object in self.cache[collection]:
            subject, object = subject_object
            if subject == node:
                values.append(str(object))
        return values
