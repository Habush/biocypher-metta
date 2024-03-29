from collections import defaultdict
import csv
import gzip
from biocypher_metta.adapters import Adapter

# Example RNAcentral bed input file:
# chr1	10244	10273	URS000035F234_9606	0	-	10244	10273	63,125,151	2	19,5	0,24	.	piRNA	PirBase
# chr1	14396	29331	URS0002033896_9606	0	-	14396	29331	63,125,151	8	433,69,970,198,510,147,99,4594	0,573,1399,2461,2836,3518,3871,10341	.	lncRNA	ENA,GeneCards
# chr1	14404	119856	URS0001A132AA_9606	0	-	14404	119856	63,125,151	12	425,69,515,159,198,510,147,99,154,904,150,211	0,565,1391,2202,2453,2828,3510,3863,10333,76321,77686,105241	.	lncRNA	ENA,GeneCards,MalaCards

# Example RNAcentral rfam annotation input file:
# URS0000000006_1317357	GO:0003735	Rfam:RF02541
# URS0000000006_1317357	GO:0005840	Rfam:RF02541
# URS0000000008_381046	GO:0030533	Rfam:RF00005

class RNACentralAdapter(Adapter):
    INDEX = {'chr': 0, 'coord_start': 1, 'coord_end': 2, 'id': 3, 'rna_type': 13}

    def __init__(self, filepath, rfam_filepath, write_properties, add_provenance):
        self.filepath = filepath
        self.rfam_filepath = rfam_filepath
        self.label = 'non_coding_rna'

        self.source = 'RNAcentral'
        self.version = '24'
        self.source_url = 'https://rnacentral.org/downloads'
        super(RNACentralAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        rfam_dictionary = defaultdict(lambda: defaultdict(set))

        with gzip.open(self.rfam_filepath, 'rt') as input:
            reader = csv.reader(input, delimiter='\t')
            for line in reader:
                rna_id, go_term, rfam = line
                if not rna_id.endswith('_9606'):
                    continue
                rna_id = rna_id.split('_')[0]
                rfam_list = [item.split(":")[1] for item in rfam.split("|")]
                rfam_dictionary[rna_id]['rfam'].update(rfam_list)
                rfam_dictionary[rna_id]['go'].add(go_term)
        
        with gzip.open(self.filepath, 'rt') as input:
            reader = csv.reader(input, delimiter='\t')
            for line in reader:
                rna_id = line[RNACentralAdapter.INDEX['id']].split('_')[0]
                props = {
                    'chr': line[RNACentralAdapter.INDEX['chr']],
                    'start': int(line[RNACentralAdapter.INDEX['coord_start']])+1, # +1 since it is 0 indexed coordinate
                    'end': int(line[RNACentralAdapter.INDEX['coord_end']])+1,
                    'rna_type': line[RNACentralAdapter.INDEX['rna_type']], 
                }
                if rna_id in rfam_dictionary:
                    props['rfam'] = list(rfam_dictionary[rna_id]['rfam'])
                    props['go_term'] = list(rfam_dictionary[rna_id]['go'])

                yield rna_id, self.label, props

