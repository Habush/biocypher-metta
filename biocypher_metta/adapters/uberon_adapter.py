import gzip
import csv
import urllib.parse
from biocypher_metta.adapters import Adapter

# Example UBERON Data from the dataset
# Class ID,Preferred Label,Synonyms,Definitions,Obsolete,CUI,Semantic Types,Parents
# http://purl.obolibrary.org/obo/UBERON_0009623,spinal nerve root,spinal root|spinal neural root|root of spinal nerve,"The paired bundles of nerve fibers entering and leaving the spinal cord at each segment. The dorsal and ventral nerve roots join to form the mixed segmental spinal nerves. The dorsal roots are generally afferent, formed by the central projections of the spinal (dorsal root) ganglia sensory cells, and the ventral roots efferent, comprising the axons of spinal motor and autonomic preganglionic neurons. There are, however, some exceptions to this afferent/efferent rule.",false,,,http://purl.obolibrary.org/obo/UBERON_0002211,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,http://purl.obolibrary.org/obo/UBERON_0002240|http://purl.obolibrary.org/obo/UBERON_0001780,,,,,,ZFA:0005578|UMLS:C0037940|SCTID:362436008|NCIT:C12791|MESH:D013126|GAID:717|FMA:14060|DHBA:146035120|BTO:0000883,,,"The paired bundles of nerve fibers entering and leaving the spinal cord at each segment. The dorsal and ventral nerve roots join to form the mixed segmental spinal nerves. The dorsal roots are generally afferent, formed by the central projections of the spinal (dorsal root) ganglia sensory cells, and the ventral roots efferent, comprising the axons of spinal motor and autonomic preganglionic neurons. There are, however, some exceptions to this afferent/efferent rule.",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,http://purl.obolibrary.org/obo/UBERON_0002240|http://purl.obolibrary.org/obo/UBERON_0001780,,,ZFA has part_of to spinal cord but this causes spatial disjointness violations (CNS/PNS),,,,,,,,,,,,,,,,,,,,,,,,spinal root|spinal neural root|root of spinal nerve,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,uberon,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,spinal nerve root,,,UBERON:0009623,,,,,UBERON:0009623,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,https://github.com/obophenotype/uberon/issues/2644,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
# http://purl.obolibrary.org/obo/UBERON_0009629,coccygeal nerve,nervus coccygeus|coccygeal spinal nerve,The coccygeal nerve is the spinal nerve that corresponds to the coccyx bone.,false,,,http://purl.obolibrary.org/obo/UBERON_0015212|http://purl.obolibrary.org/obo/UBERON_0001780,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,neuronames:2369|Wikipedia:Coccygeal_nerve|UMLS:C0228967|SCTID:287477006|NCIT:C32333|FMA:5863,,,The coccygeal nerve is the spinal nerve that corresponds to the coccyx bone.,,,http://upload.wikimedia.org/wikipedia/commons/d/d1/Gray802.png,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,http://purl.obolibrary.org/obo/UBERON_0005845,,,,,,,,,,,,,,,,,,,,,,,,,,,coccygeal spinal nerve,,,,,,,,,,,,,,,,,,,,,,,,,,,nervus coccygeus,,,,,,,,,,,,,,uberon,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,coccygeal nerve,,,UBERON:0009629,,,,,UBERON:0009629,,,,,,,,,,,,,,,,,http://purl.obolibrary.org/obo/UBERON_0000010,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
# http://purl.obolibrary.org/obo/UBERON_0009621,tail somite,,A somite that is part of a tail.,false,,,http://purl.obolibrary.org/obo/UBERON_0002329,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,XAO:0003143|EMAPA:16860|AAO:0010383,,,A somite that is part of a tail.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,uberon,,,,,,,,,,,,http://purl.obolibrary.org/obo/UBERON_0002415,,http://purl.obolibrary.org/obo/UBERON_0002415,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,tail somite,,,UBERON:0009621,,,,,UBERON:0009621,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
# http://purl.obolibrary.org/obo/UBERON_0009624,lumbar nerve,nervi lumbales|nervus lumbalis|lumbar spinal nerve,The lumbar nerves are the five spinal nerves emerging from the lumbar vertebrae. They are divided into posterior and anterior divisions.,false,,,http://purl.obolibrary.org/obo/UBERON_0001780,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,neuronames:1730|Wikipedia:Lumbar_nerves|UMLS:C0228897|SCTID:361600002|NCIT:C33015|FMA:5861,,,The lumbar nerves are the five spinal nerves emerging from the lumbar vertebrae. They are divided into posterior and anterior divisions.,,,http://upload.wikimedia.org/wikipedia/commons/d/d1/Gray802.png,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,http://purl.obolibrary.org/obo/UBERON_0002792,,,,,,,,,,,,,,,,,,,,,,,,,,,nervus lumbalis|lumbar spinal nerve,,,,,,,,,,,,,,,,,,,,,,,,,,,nervi lumbales,,,,,,,,,,,,,,uberon,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,lumbar nerve,,,UBERON:0009624,,,,,UBERON:0009624,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

class UberonAdapter(Adapter):
    """
    Adapter for UBERON dataset
    """

    def __init__(self, filepath, write_properties, add_provenance, label="uberon"):
        self.filepath = filepath
        self.label = label
        self.source = "Uberon"
        self.source_url = "https://bioportal.bioontology.org/ontologies/UBERON"
        super(UberonAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, "rt") as f:
            reader = csv.DictReader(f, delimiter=",", quotechar='"')

            for row in reader:
                try:
                    node_id = row['Class ID'].split('/')[-1]
                    node_url = urllib.parse.unquote(row['Class ID'])
                    node_label = row['Preferred Label']
                    is_obsolete = row['Obsolete'].lower() == "true"
                    synonyms = [s for s in row['Synonyms'].split('|') if s]

                    props = {}
                    if self.write_properties:
                        props['preferred_label'] = node_label
                        props['synonyms'] = synonyms
                        props['is_obsolete'] = is_obsolete
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

    