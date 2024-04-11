import csv
import gzip
import pickle
from biocypher_metta.adapters import Adapter
from liftover import get_lifter
# Example dbSuper tsv input files:
# chrom	 start	 stop	 se_id	 gene_symbol	 cell_name	 rank
# chr1	120485363	120615071	SE_00001	NOTCH2	Adipose Nuclei	1
# chr13	110838836	111112228	SE_00002	COL4A2	Adipose Nuclei	2
# chr1	145206326	145293008	SE_00003	NOTCH2NL	Adipose Nuclei	3
# chr5	158117077	158371526	SE_00004	EBF1	Adipose Nuclei	4

class DBSuperAdapter(Adapter):
    INDEX = {'chr' : 0, 'coord_start' : 1, 'coord_end' : 2, 'se_id' : 3, 'gene_id' : 4}

    def __init__(self, filepath, hgnc_to_ensembl_map, write_properties, add_provenance, type='super enhancer', label='super_enhancer', delimiter='\t'):
        self.filePath = filepath
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.type = type
        self.label = label
        self.delimiter = delimiter
        self.genome_reference_converter = get_lifter('hg19', 'hg38')

        self.source = 'dbSuper'
        self.version = ''
        self.source_url = 'https://asntech.org/dbsuper/download.php'

        super(DBSuperAdapter, self).__init__(write_properties, add_provenance)

    def convert_to_hg38(self, chr,  pos):
        try:
            chr_no = chr.replace('chr', '').replace('ch', '')
            converted = self.genome_reference_converter.query(chr_no, pos)[0][1]
            return int(converted)
        except:
            return False

    def get_nodes(self):
        with gzip.open(self.filePath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            next(reader)
            for line in reader:
                chr = line[DBSuperAdapter.INDEX['chr']]
                start_position = self.convert_to_hg38(chr, int(line[DBSuperAdapter.INDEX['coord_start']]))
                end_position = self.convert_to_hg38(chr, int(line[DBSuperAdapter.INDEX['coord_end']]))
                if not start_position or not end_position:
                    continue
                se_id = line[DBSuperAdapter.INDEX['se_id']]
                
                props = {}
                if self.write_properties:
                    props['chr'] = chr
                    props['start'] = start_position
                    props['end'] = end_position
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield se_id, self.label, props

    def get_edges(self):
        with gzip.open(self.filePath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            next(reader)
            for line in reader:
                se_id = line[DBSuperAdapter.INDEX['se_id']]
                gene_id = line[DBSuperAdapter.INDEX['gene_id']]
                ensembl_gene_id = self.hgnc_to_ensembl_map.get(gene_id, gene_id)
                
                props = {}
                if self.write_properties:
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield se_id, ensembl_gene_id, self.label, props