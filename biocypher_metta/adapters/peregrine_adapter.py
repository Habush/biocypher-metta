import csv
import gzip
import pickle

from biocypher_metta.adapters import Adapter
# Example PEREGRINE input files:

# PEREGRINEenhancershg38
# chr1	99534632	99534837	1
# chr1	99638888	99639130	2
# chr1	99673438	99673767	5

# PREGRINE scores
# enhancer	gene	HepG2score	HepG2Zscore	CDF_HepG2	HCT116score	HCT116Zscore	CDF_HCT116	K562score	K562Zscore	CDF_K562	MCF7score	MCF7Zscore	CDF_MCF7
# 1	HUMAN|HGNC=15846|UniProtKB=Q9NP74	0.14231604149762	0.364411398538007	0.799181060202044	0.1576516981884	-0.337757438945174	0.873174506724474	0.0344599217044222	-0.337757438945174	0.0254839982033111	0.17554256461814	0.733623642220683	0.906677165599342
# 1	HUMAN|HGNC=321|UniProtKB=P35573	0.0430095105578992	-0.563353003049424	0.241575024295348	0.053969581878424	-0.242133777354291	0.369753681325124	0.0462475325608546	-0.242133777354291	0.504521945621637	0.133478871732277	0.259622696593164	0.71453419341604

# PEREGRINE datasources
# 1	FANTOM
# 2	FANTOM
# 87657	Ensembl
# 87659	Ensembl
# EH37E0436909	ENCODE
# EH37E0436910	ENCODE

class PEREGRINEAdapter(Adapter):
  ALLOWED_TYPES = ['enhancer', 'peregrine enhancer to gene association']
  ALLOWED_LABELS = ['enhancer', 'enhancer_regulates']
  ALLOWED_KEYS = []
  INDEX = {'enhancer': 0,	'gene': 1,	'HepG2score': 2,	'HepG2Zscore': 3,	'CDF_HepG2': 4,	'HCT116score': 5,	'HCT116Zscore': 6,	'CDF_HCT116': 7,	'K562score': 8,	'K562Zscore': 9,	'CDF_K562': 10,	'MCF7score': 11,	'MCF7Zscore': 12,	'CDF_MCF7': 13}

  def __init__(self, enhancers_file, enhancer_gene_link_scores, source_file, hgnc_ensembl_map, hgnc_symbol_map, type='enhancer', label='enhancer', delimiter='\t'):
    self.enhancers_file = enhancers_file
    self.enhancer_gene_link_scores = enhancer_gene_link_scores
    self.source_file = source_file
    self.hgnc_ensembl_map = pickle.load(open(hgnc_ensembl_map, 'rb'))
    self.hgnc_symbol_map = pickle.load(open(hgnc_symbol_map, 'rb'))
    self.type = type
    self.label = label
    self.delimiter = delimiter

    self.source = 'PEREGRINE'
    self.version = ''
    self.source_url = 'https://www.peregrineproj.org/'

  def handel_gene(self, gene):
    gene = gene.split('|')[1]
    gene = ':'.join(gene.split('='))
    if gene in self.hgnc_ensembl_map:
      gene = self.hgnc_ensembl_map[gene]
    else:
      gene = self.hgnc_symbol_map.get(gene, gene)

    return gene

  def get_nodes(self):
    enhancer_info = {}
    with gzip.open(self.enhancers_file, 'rt') as f:
      reader = csv.reader(f, delimiter=self.delimiter)
      for line in reader:
        chr, start, end, enhancer_id = line
        enhancer_info[enhancer_id] = {
          'chr': chr,
          'start': start,
          'end': end
        }

    target_genes = {}
    with gzip.open(self.enhancer_gene_link_scores, 'rt') as f:
      reader = csv.reader(f, delimiter=self.delimiter)
      next(reader)  # Skip header
      for line in reader:
        enhancer_id = line[self.INDEX['enhancer']]
        gene = line[self.INDEX['gene']]
        if enhancer_id not in target_genes:
          target_genes[enhancer_id] = set()
        target_genes[enhancer_id].add(self.handel_gene(gene))

    source_map = {}
    with gzip.open(self.source_file, 'rt') as f:
      for line in f:
        enhancer_id, source = line.strip().split('\t')
        source_map[enhancer_id] = source

    for enhancer_id, info in enhancer_info.items():
      gene_list = target_genes.get(enhancer_id, [])
      data_source = source_map.get(enhancer_id, 'NA')
      props = {
          'chr': info['chr'],
          'start': info['start'],
          'end': info['end'],
      }
      if gene_list:
        props['genes'] = list(gene_list)
      props['data_source'] = data_source
      
      yield enhancer_id, self.label, props

  def get_edges(self):
    with gzip.open(self.enhancer_gene_link_scores, 'rt') as f:
      reader = csv.reader(f, delimiter=self.delimiter)
      next(reader)  # Skip header
      for line in reader:
        enhancer_id = line[self.INDEX['enhancer']]
        gene = self.handel_gene(line[self.INDEX['gene']])
        uniprotkb = line[self.INDEX['gene']].split('|')[2].split('=')[1]
        HepG2score	= line[self.INDEX['HepG2score']]
        HepG2Zscore	= line[self.INDEX['HepG2Zscore']]
        CDF_HepG2	= line[self.INDEX['CDF_HepG2']]
        HCT116score	= line[self.INDEX['HCT116score']]
        HCT116Zscore= line[self.INDEX['HCT116Zscore']]	
        CDF_HCT116	= line[self.INDEX['CDF_HCT116']]
        K562score	= line[self.INDEX['K562score']]
        K562Zscore	= line[self.INDEX['K562Zscore']]
        CDF_K562	= line[self.INDEX['CDF_K562']]
        MCF7score	= line[self.INDEX['MCF7score']]
        MCF7Zscore	= line[self.INDEX['MCF7Zscore']]
        CDF_MCF7= line[self.INDEX['CDF_MCF7']]

        props = {
          'uniprotkb': uniprotkb,
          'HepG2': ['HepG2score:'+HepG2score, 'HepG2Zscore:'+HepG2Zscore, 'CDF_HepG2:'+CDF_HepG2],
          'HCT116': ['HCT116score:'+HCT116score, 'HCT116Zscore:'+HCT116Zscore, 'CDF_HCT116:'+CDF_HCT116],
          'K562': ['K562score:'+K562score, 'K562Zscore:'+K562Zscore, 'CDF_K562:'+CDF_K562],
          'MCF7': ['MCF7score:'+MCF7score, 'MCF7Zscore:'+MCF7Zscore, 'CDF_MCF7:'+CDF_MCF7]
        }
        
      
        yield enhancer_id, gene, self.label, props
      
