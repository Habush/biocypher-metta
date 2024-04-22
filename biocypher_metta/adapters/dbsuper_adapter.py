import csv
import gzip
import pickle
from biocypher_metta.adapters import Adapter

from biocypher_metta.adapters.helpers import build_regulatory_region_id, check_genomic_location, convert_genome_reference
# Example dbSuper tsv input files:
# chrom	 start	 stop	 se_id	 gene_symbol	 cell_name	 rank
# chr1	120485363	120615071	SE_00001	NOTCH2	Adipose Nuclei	1
# chr13	110838836	111112228	SE_00002	COL4A2	Adipose Nuclei	2
# chr1	145206326	145293008	SE_00003	NOTCH2NL	Adipose Nuclei	3
# chr5	158117077	158371526	SE_00004	EBF1	Adipose Nuclei	4

class DBSuperAdapter(Adapter):
    INDEX = {'chr' : 0, 'coord_start' : 1, 'coord_end' : 2, 'se_id' : 3, 'gene_id' : 4}

    def __init__(self, filepath, hgnc_to_ensembl_map, write_properties, add_provenance, 
                 type='super enhancer', label='super_enhancer', delimiter='\t',
                 chr=None, start=None, end=None):
        self.filePath = filepath
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.type = type
        self.label = label
        self.delimiter = delimiter
        self.chr = chr
        self.start = start
        self.end = end

        self.source = 'dbSuper'
        self.version = ''
        self.source_url = 'https://asntech.org/dbsuper/download.php'

        super(DBSuperAdapter, self).__init__(write_properties, add_provenance)


    def get_nodes(self):
        with gzip.open(self.filePath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            next(reader)
            for line in reader:
                se_id = line[DBSuperAdapter.INDEX['se_id']]
                chr = line[DBSuperAdapter.INDEX['chr']]
                start_hg19 = int(line[DBSuperAdapter.INDEX['coord_start']]) + 1 # +1 since it is 0-based genomic coordinate
                end_hg19 = int(line[DBSuperAdapter.INDEX['coord_end']]) + 1
                start = convert_genome_reference(chr, start_hg19)
                end = convert_genome_reference(chr, end_hg19)
                
                if not start or not end:
                    continue
                se_region_id = build_regulatory_region_id(chr, start, end)
                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    props = {}
                    if self.write_properties:
                        props['se_id'] = se_id
                        props['chr'] = chr
                        props['start'] = start
                        props['end'] = end
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield se_region_id, self.label, props

    def get_edges(self):
        with gzip.open(self.filePath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            next(reader)
            for line in reader:
                gene_id = line[DBSuperAdapter.INDEX['gene_id']]
                ensembl_gene_id = self.hgnc_to_ensembl_map.get(gene_id, None)
                chr = line[DBSuperAdapter.INDEX['chr']]
                start_hg19 = int(line[DBSuperAdapter.INDEX['coord_start']]) + 1 # +1 since it is 0-based genomic coordinate
                end_hg19 = int(line[DBSuperAdapter.INDEX['coord_end']]) + 1
                start = convert_genome_reference(chr, start_hg19)
                end = convert_genome_reference(chr, end_hg19)
                
                if not ensembl_gene_id or not start or not end:
                    continue
                se_region_id = build_regulatory_region_id(chr, start, end)
                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    props = {}
                    if self.write_properties:
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield se_region_id, ensembl_gene_id, self.label, props