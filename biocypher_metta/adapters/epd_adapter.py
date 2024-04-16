import csv
import gzip
import pickle
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_regulatory_region_id
# Example EPD bed input file:
##CHRM Start  End   Id  Score Strand -  -
# chr1 959245 959305 NOC2L_1 900 - 959245 959256
# chr1 960583 960643 KLHL17_1 900 + 960632 960643
# chr1 966432 966492 PLEKHN1_1 900 + 966481 966492
# chr1 976670 976730 PERM1_1 900 - 976670 976681

class EPDAdapter(Adapter):
    INDEX = {'chr' : 0, 'coord_start' : 1, 'coord_end' : 2, 'gene_id' : 3}

    def __init__(self, filepath, hgnc_to_ensembl_map, write_properties, add_provenance, type='promoter', label='promoter', delimiter=' '):
        self.filepath = filepath
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.type = type
        self.label = label
        self.delimiter = delimiter

        self.source = 'EPD'
        self.version = '006'
        self.source_url = 'https://epd.expasy.org/ftp/epdnew/H_sapiens/'

        super(EPDAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                chr = line[EPDAdapter.INDEX['chr']]
                coord_start = int(line[EPDAdapter.INDEX['coord_start']]) + 1 # +1 since it is 0 indexed coordinate
                coord_end = int(line[EPDAdapter.INDEX['coord_end']]) + 1
                promoter_id = build_regulatory_region_id(chr, coord_start, coord_end)

                props = {}
                if self.write_properties:
                    props['chr'] = chr
                    props['start'] = coord_start
                    props['end'] = coord_end

                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield promoter_id, self.label, props

    def get_edges(self):
        with gzip.open(self.filepath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                chr = line[EPDAdapter.INDEX['chr']]
                coord_start = int(line[EPDAdapter.INDEX['coord_start']]) + 1 # +1 since it is 0 indexed coordinate
                coord_end = int(line[EPDAdapter.INDEX['coord_end']]) + 1
                gene_id = line[EPDAdapter.INDEX['gene_id']].split('_')[0]
                ensembl_gene_id = self.hgnc_to_ensembl_map.get(gene_id, None)
                if ensembl_gene_id is None:
                    continue
                promoter_id = build_regulatory_region_id(chr, coord_start, coord_end)

                props = {}
                if self.write_properties:
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield promoter_id, ensembl_gene_id, self.label, props
