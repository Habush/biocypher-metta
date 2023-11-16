import gzip
from biocypher_metta.adapters import Adapter
# Example genocde vcf input file:
# ##description: evidence-based annotation of the human genome (GRCh38), version 42 (Ensembl 108)
# ##provider: GENCODE
# ##contact: gencode-help@ebi.ac.uk
# ##format: gtf
# ##date: 2022-07-20
# chr1    HAVANA  gene    11869   14409   .       +       .       gene_id "ENSG00000290825.1"; gene_type "lncRNA"; gene_name "DDX11L2"; level 2; tag "overlaps_pseudogene";
# chr1    HAVANA  transcript      11869   14409   .       +       .       gene_id "ENSG00000290825.1"; transcript_id "ENST00000456328.2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
# chr1    HAVANA  exon    11869   12227   .       +       .       gene_id "ENSG00000290825.1"; transcript_id "ENST00000456328.2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
# chr1    HAVANA  exon    12613   12721   .       +       .       gene_id "ENSG00000290825.1"; transcript_id "ENST00000456328.2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";


class GencodeGeneAdapter(Adapter):
    ALLOWED_KEYS = ['gene_id', 'gene_type', 'gene_name',
                    'transcript_id', 'transcript_type', 'transcript_name', 'hgnc_id']
    INDEX = {'chr': 0, 'type': 2, 'coord_start': 3, 'coord_end': 4, 'info': 8}

    def __init__(self, filepath=None, gene_alias_file_path=None, chr='all', dry_run=False):

        self.filepath = filepath
        self.chr = chr
        self.label = 'gene'
        self.dataset = 'gencode_gene'
        self.gene_alias_file_path = gene_alias_file_path
        self.dry_run = dry_run

    def parse_info_metadata(self, info):
        parsed_info = {}
        for key, value in zip(info, info[1:]):
            if key in GencodeGeneAdapter.ALLOWED_KEYS:
                parsed_info[key] = value.replace('"', '').replace(';', '')
        return parsed_info

    # the gene alias dict will use both ensembl id and hgnc id as key
    def get_gene_alias(self):
        alias_dict = {}
        with gzip.open(self.gene_alias_file_path, 'rt') as input:
            next(input)
            for line in input:
                (tax_id, gene_id, symbol, locus_tag, synonyms, dbxrefs, chromosome, map_location, description, type_of_gene, symbol_from_nomenclature_authority,
                 full_name_from_nomenclature_authority, Nomenclature_status, Other_designations, Modification_date, Feature_type) = line.split('\t')

                split_dbxrefs = dbxrefs.split('|')
                hgnc = ''
                ensembl = ''
                for ref in split_dbxrefs:
                    if ref.startswith('HGNC:'):
                        hgnc = ref[5:]
                    if ref.startswith('Ensembl:'):
                        ensembl = ref[8:]
                if ensembl or hgnc:
                    complete_synonyms = []
                    complete_synonyms.append(symbol)
                    for i in synonyms.split('|'):
                        complete_synonyms.append(i)
                    if hgnc:
                        complete_synonyms.append(hgnc)
                    for i in Other_designations.split('|'):
                        complete_synonyms.append(i)
                    complete_synonyms.append(
                        symbol_from_nomenclature_authority)
                    complete_synonyms.append(
                        full_name_from_nomenclature_authority)
                    complete_synonyms = list(set(complete_synonyms))
                    if '-' in complete_synonyms:
                        complete_synonyms.remove('-')
                    if ensembl:
                        alias_dict[ensembl] = complete_synonyms
                    if hgnc:
                        alias_dict[hgnc] = complete_synonyms

        return alias_dict

    def get_nodes(self):
        alias_dict = self.get_gene_alias()
        for line in open(self.filepath, 'r'):
            if line.startswith('#'):
                continue
            split_line = line.strip().split()
            if split_line[GencodeGeneAdapter.INDEX['type']] == 'gene':
                info = self.parse_info_metadata(
                    split_line[GencodeGeneAdapter.INDEX['info']:])
                gene_id = info['gene_id']
                id = gene_id.split('.')[0]
                alias = alias_dict.get(id)
                if not alias:
                    hgnc_id = info.get('hgnc_id')
                    if hgnc_id:
                        alias = alias_dict.get(hgnc_id)
                if gene_id.endswith('_PAR_Y'):
                    id = id + '_PAR_Y'

                props = {
                    # 'gene_id': gene_id, # TODO should this be included?
                    'gene_type': info['gene_type'],
                    'chr': split_line[GencodeGeneAdapter.INDEX['chr']],
                    # the gtf file format is [1-based,1-based], needs to convert to BED format [0-based,1-based]
                    'start': int(split_line[GencodeGeneAdapter.INDEX['coord_start']]) - 1,
                    'end': int(split_line[GencodeGeneAdapter.INDEX['coord_end']]),
                    'gene_name': info['gene_name'],
                    'synonyms': alias,
                    'source': 'GENCODE',
                    'version': 'v43',
                    'source_url': 'https://www.gencodegenes.org/human/'
                }

                yield id, self.label, props
