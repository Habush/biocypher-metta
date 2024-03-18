import csv
import gzip
import pickle
from biocypher_metta.adapters import Adapter
# Example EPD bed input file:
# chr1 959245 959305 NOC2L_1 900 - 959245 959256
# chr1 960583 960643 KLHL17_1 900 + 960632 960643
# chr1 966432 966492 PLEKHN1_1 900 + 966481 966492
# chr1 976670 976730 PERM1_1 900 - 976670 976681

class EPDAdapter(Adapter):
    INDEX = {'chr' : 0, 'coord_start' : 1, 'coord_end' : 2, 'gene_id' : 3}

    def __init__(self, filepath, hgnc_to_ensembl_map, label='promoter', delimiter=' '):
        self.filepath = filepath
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.label = label
        self.delimiter = delimiter

        self.source = 'EPD'
        self.version = '006'
        self.source_url = 'https://epd.expasy.org/ftp/epdnew/H_sapiens/'

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                chr = line[EPDAdapter.INDEX['chr']]
                coord_start = line[EPDAdapter.INDEX['coord_start']]
                coord_end = line[EPDAdapter.INDEX['coord_end']]
                gene_id = line[EPDAdapter.INDEX['gene_id']]

                gene_id = gene_id.split('_')[0]
                gene_id = self.hgnc_to_ensembl_map.get(gene_id, gene_id)

                promoter_id = chr+":"+coord_start+"-"+coord_end
                props = {'chr' : chr, 'start' : coord_start, 'end' : coord_end, 'gene' : gene_id}

                yield promoter_id, self.label, props
