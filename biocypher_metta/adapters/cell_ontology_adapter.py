from biocypher_metta.adapters import Adapter
import gzip
import csv
import urllib.parse

# Example Cell Ontology Data from the dataset
# Class ID,Preferred Label,Synonyms,Definitions,Obsolete,CUI,Semantic Types,Parents
# http://purl.obolibrary.org/obo/UBERON_0009623	spinal nerve root	spinal root|spinal neural root|root of spinal nerve	The paired bundles of nerve fibers entering and leaving the spinal cord at each segment. The dorsal and ventral nerve roots join to form the mixed segmental spinal nerves. The dorsal roots are generally afferent, formed by the central projections of the spinal (dorsal root) ganglia sensory cells, and the ventral roots efferent, comprising the axons of spinal motor and autonomic preganglionic neurons. There are, however, some exceptions to this afferent/efferent rule.	FALSE			http://purl.obolibrary.org/obo/UBERON_0002211									
# http://purl.obolibrary.org/obo/UBERON_0009621	tail somite		A somite that is part of a tail.	FALSE			http://purl.obolibrary.org/obo/UBERON_0002329									
# http://purl.obolibrary.org/obo/UBERON_0009622	pronephric proximal straight tubule	proximal straight tubules	A proximal straight tubule that is part of a pronephros.	FALSE			http://purl.obolibrary.org/obo/UBERON_0005310|http://purl.obolibrary.org/obo/UBERON_0001290									
# http://purl.obolibrary.org/obo/UBERON_0007238	1st arch maxillary component			FALSE			http://purl.obolibrary.org/obo/UBERON_0002050				
# http://purl.obolibrary.org/obo/UBERON_0007239	tunica media of artery	tunica media (arteriae)|arterial media	A tunica media that is part of a artery.	FALSE			http://purl.obolibrary.org/obo/UBERON_0002522				


class CellOntologyAdapter(Adapter):
    """
    Adapter for Biomedical Ontology dataset
    """

    def __init__(self, filepath, write_properties, add_provenance, label="cell", type="cell"):
        self.filepath = filepath
        self.label = label
        self.type = type
        self.source = "Cell Ontology"
        self.source_url = "https://bioportal.bioontology.org/ontologies/CL"
        super(CellOntologyAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.DictReader(f, delimiter=",", quotechar='"')

            for row in reader:
                node_id = row['Class ID'].split('/')[-1]
                node_url = urllib.parse.unquote(row['Class ID'])
                node_label = row['Preferred Label']
                is_obsolete = row['Obsolete'].lower() == "true"
                parent_urls = [urllib.parse.unquote(p) for p in row['Parents'].split('|') if p]
                properties = {
                    "label": node_label,
                    "is_obsolete": is_obsolete,
                    "url": node_url,
                    "parent_urls": parent_urls
                }
                yield node_id, self.label, properties

    def get_edges(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.DictReader(f, delimiter=",", quotechar='"')

            for row in reader:
                node_id = row['Class ID'].split('/')[-1]
                parent_urls = [urllib.parse.unquote(p) for p in row['Parents'].split('|') if p]

                if len(parent_urls) > 1:
                    parent_ids = [parent_url.split('/')[-1] for parent_url in parent_urls]
                    if parent_ids:
                        parent_ids_str = ' '.join(parent_ids)
                    yield node_id, parent_ids_str, self.label, {}
                elif len(parent_urls) == 1:
                    parent_id = parent_urls[0].split('/')[-1]
                    yield node_id, parent_id, self.label, {}