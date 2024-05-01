import csv
import gzip
import pickle

from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_regulatory_region_id, check_genomic_location
# Example PEREGRINE input files:

# PEREGRINEenhancershg38
# chr1	99534632	99534837	1
# chr1	99638888	99639130	2
# chr1	99673438	99673767	5

# enhancer_gene_link_18.tsv
# enhancer	gene	linkID	assay	tissue	p-value	eQTL_SNP_ID	score	Zscore	CDFscore
# 1	HUMAN|HGNC=15846|UniProtKB=Q9NP74	1	3	64		
# 1	HUMAN|HGNC=15846|UniProtKB=Q9NP74	1	3	65		
# 1	HUMAN|HGNC=15846|UniProtKB=Q9NP74	1	3	66		
# 1	HUMAN|HGNC=15846|UniProtKB=Q9NP74	1	3	67		

# PEREGRINE datasources
# 1	FANTOM
# 2	FANTOM
# 87657	Ensembl
# 87659	Ensembl
# EH37E0436909	ENCODE
# EH37E0436910	ENCODE

class PEREGRINEAdapter(Adapter):
    ALLOWED_TYPES = ['enhancer', 'enhancer to gene association']
    ALLOWED_LABELS = ['enhancer', 'enhancer_gene']
    ALLOWED_KEYS = []
    INDEX = {'enhancer': 0, 'gene': 1, 'tissue': 4, 'score': 7, 'chr': 0, 'start': 1, 'end': 2, 'id': 3}

    def __init__(self, enhancers_file, enhancer_gene_link, 
                 source_file, hgnc_ensembl_map, 
                 tissue_ontology_map, write_properties, add_provenance, 
                 type='enhancer', label='enhancer', delimiter='\t',
                 chr=None, start=None, end=None):
        
        self.enhancers_file = enhancers_file
        self.enhancer_gene_link = enhancer_gene_link
        self.source_file = source_file
        self.hgnc_ensembl_map = pickle.load(open(hgnc_ensembl_map, 'rb'))
        self.tissue_ontology_map = pickle.load(open(tissue_ontology_map, 'rb'))
        self.type = type
        self.label = label
        self.delimiter = delimiter
        self.chr = chr
        self.start = start
        self.end = end

        self.source = 'PEREGRINE'
        self.version = ''
        self.source_url = 'https://www.peregrineproj.org/'

        super(PEREGRINEAdapter, self).__init__(write_properties, add_provenance)

    def handle_gene(self, gene):
        gene = gene.split('|')[1]
        gene = ':'.join(gene.split('='))
        return gene
        
    def get_nodes(self):
        enhancer_info = {}
        with gzip.open(self.enhancers_file, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                chr, start, end, enhancer_id = line
                enhancer_info[enhancer_id] = {
                    'chr': chr,
                    'start': int(start),
                    'end': int(end),
                }

        source_map = {}
        with gzip.open(self.source_file, 'rt') as f:
            for line in f:
                enhancer_id, source = line.strip().split('\t')
                source_map[enhancer_id] = source

        for enhancer_id, info in enhancer_info.items():
            data_source = source_map[enhancer_id]
            chr = info['chr']
            start = info['start']
            end = info['end']
            enhancer_region_id = build_regulatory_region_id(chr, start, end)
            props = {}
            if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                if self.write_properties:
                    props['id'] = enhancer_id
                    props['chr'] = chr
                    props['start'] = start
                    props['end'] = end
                    props['data_source'] = data_source

                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url
                
                yield enhancer_region_id, self.label, props

    def get_edges(self):
        enhancer_id_map = {}
        with gzip.open(self.enhancers_file, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                id = line[self.INDEX['id']]
                chr = line[self.INDEX['chr']]
                start = int(line[self.INDEX['start']])
                end = int(line[self.INDEX['end']])
                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    region_id = build_regulatory_region_id(chr, start, end)
                    enhancer_id_map[id] = region_id
        
        with gzip.open(self.enhancer_gene_link, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            next(reader)    # Skip header
            for line in reader:
                id = line[self.INDEX['enhancer']]
                if id not in enhancer_id_map:
                    continue
                enhancer_region_id = enhancer_id_map[id]
                gene_hgnc_id = self.handle_gene(line[self.INDEX['gene']])
                if gene_hgnc_id not in self.hgnc_ensembl_map:
                    continue
                gene = self.hgnc_ensembl_map[gene_hgnc_id]
                tissue_id = line[self.INDEX['tissue']]
                score = None
                if self.INDEX['score'] < len(line):
                    score = line[self.INDEX['score']]
                if tissue_id not in self.tissue_ontology_map:
                    continue

                props = {}
                if self.write_properties:
                    props['biological_context'] = self.tissue_ontology_map[tissue_id][0]
                    if score:
                        props['score'] = score

                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url
                
                yield enhancer_region_id, gene, self.label, props
            
