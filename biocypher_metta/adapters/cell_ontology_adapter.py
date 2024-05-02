from biocypher_metta.adapters import Adapter
import gzip
import csv
import urllib.parse

# Example Cell Ontology Data from the dataset
# Class ID,Preferred Label,Synonyms,Definitions,Obsolete,CUI,Semantic Types,Parents
# http://purl.obolibrary.org/obo/UBERON_0009623,spinal nerve root,spinal root|spinal neural root|root of spinal nerve,"The paired bundles of nerve fibers entering and leaving the spinal cord at each segment. The dorsal and ventral nerve roots join to form the mixed segmental spinal nerves. The dorsal roots are generally afferent, formed by the central projections of the spinal (dorsal root) ganglia sensory cells, and the ventral roots efferent, comprising the axons of spinal motor and autonomic preganglionic neurons. There are, however, some exceptions to this afferent/efferent rule.",false,,,http://purl.obolibrary.org/obo/UBERON_0002211
# http://purl.obolibrary.org/obo/UBERON_0009621,tail somite,,A somite that is part of a tail.,false,,,http://purl.obolibrary.org/obo/UBERON_0002329
# http://purl.obolibrary.org/obo/UBERON_0009622,pronephric proximal straight tubule,proximal straight tubules,A proximal straight tubule that is part of a pronephros.,false,,,http://purl.obolibrary.org/obo/UBERON_0005310|http://purl.obolibrary.org/obo/UBERON_0001290
# http://purl.obolibrary.org/obo/PR_000001244,long-wave-sensitive opsin 1,RCP|red-sensitive opsin|red cone photoreceptor pigment|ROP|OPN1LW,An animal opsin that is a translation product of the human OPN1LW gene or a 1:1 ortholog thereof.|Category=gene.,false,,,http://purl.obolibrary.org/obo/PR_000001119
# http://purl.obolibrary.org/obo/PR_000001243,melanopsin,Mopn|MOP|opsin-4|OPN4,An animal opsin that is a translation product of the human OPN4 gene or a 1:1 ortholog thereof.|Category=gene.,false,,,http://purl.obolibrary.org/obo/PR_000001119
# http://purl.obolibrary.org/obo/GO_0045970,negative regulation of juvenile hormone catabolic process,inhibition of juvenile hormone catabolic process|negative regulation of juvenile hormone degradation|negative regulation of juvenile hormone catabolism|negative regulation of juvenile hormone breakdown|downregulation of juvenile hormone catabolic process|down-regulation of juvenile hormone catabolic process|down regulation of juvenile hormone catabolic process,"Any process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways resulting in the breakdown of juvenile hormone.",false,,,http://purl.obolibrary.org/obo/GO_0065007|http://purl.obolibrary.org/obo/GO_0045952|http://purl.obolibrary.org/obo/GO_0045928|http://purl.obolibrary.org/obo/GO_0009895


class CellOntologyAdapter(Adapter):
    """
    Adapter for Cell Ontology dataset
    """

    def __init__(self, filepath, write_properties, add_provenance, label='cell'):
        self.filepath = filepath
        self.label = label
        self.source = "Cell Ontology"
        self.source_url = "https://bioportal.bioontology.org/ontologies/CL"
        super(CellOntologyAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.DictReader(f, delimiter=",", quotechar='"')

            for row in reader:
                try:
                    node_id = row['Class ID'].split('/')[-1]
                    node_url = urllib.parse.unquote(row['Class ID'])
                    node_label = row['Preferred Label']
                    is_obsolete = row['Obsolete'].lower() == "true"
                
                    props = {}
                    if self.write_properties:
                        props['preferred_label'] = node_label
                        props['is_obsolate'] = is_obsolete
                        props['url'] = node_url

                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url
                
                    yield node_id, self.label, props
                except Exception as e:
                    print(f"Error processing row: {row}. Error: {e}")

    def get_edges(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.DictReader(f, delimiter=",", quotechar='"')

            for row in reader:
                try:
                    node_id = row['Class ID'].split('/')[-1]
                    parent_urls = [urllib.parse.unquote(p) for p in row['Parents'].split('|') if p]

                    props = {}
                    if self.write_properties:
                        props['parent_urls'] = parent_urls

                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    for parent_url in parent_urls:
                        parent_id = parent_url.split('/')[-1]
                        yield parent_id, node_id, self.label, props
                except Exception as e:
                    print(f"Error processing row: {row}. Error: {e}")