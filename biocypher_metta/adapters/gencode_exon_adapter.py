import gzip
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import check_genomic_location

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


class GencodeExonAdapter(Adapter):
    ALLOWED_KEYS = ['gene_id', 'transcript_id', 'transcript_type', 'transcript_name', 'exon_number', 'exon_id']
    INDEX = {'chr': 0, 'type': 2, 'coord_start': 3, 'coord_end': 4, 'info': 8}

    def __init__(self, write_properties, add_provenance, filepath=None,
                 chr=None, start=None, end=None):
        self.filepath = filepath
        self.chr = chr
        self.start = start
        self.end = end
        self.label = 'exon'
        self.dataset = 'gencode_exon'
        self.source = 'GENCODE'
        self.version = 'v44'
        self.source_url = 'https://www.gencodegenes.org/human/'

        super(GencodeExonAdapter, self).__init__(write_properties, add_provenance)

    def parse_info_metadata(self, info):
        parsed_info = {}
        for key, value in zip(info, info[1:]):
            if key in GencodeExonAdapter.ALLOWED_KEYS:
                parsed_info[key] = value.replace('"', '').replace(';', '')
        return parsed_info

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as input:
            for line in input:
                if line.startswith('#'):
                    continue
                split_line = line.strip().split()
                if split_line[GencodeExonAdapter.INDEX['type']] == 'exon':
                    info = self.parse_info_metadata(
                        split_line[GencodeExonAdapter.INDEX['info']:])
                    gene_id = info['gene_id'].split('.')[0]
                    transcript_id = info['transcript_id'].split('.')[0]
                    exon_id = info['exon_id'].split('.')[0]
                    chr = split_line[GencodeExonAdapter.INDEX['chr']]
                    start = int(split_line[GencodeExonAdapter.INDEX['coord_start']])
                    end = int(split_line[GencodeExonAdapter.INDEX['coord_end']])
                    props = {}
                    try:
                        if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                            if self.write_properties: 
                                props = {
                                    'gene_id': gene_id,
                                    'transcript_id': transcript_id,
                                    'chr': chr,
                                    'start': start,
                                    'end': end,
                                    'exon_number': int(info.get('exon_number', -1)), 
                                    'exon_id': exon_id
                                }
                                if self.add_provenance:
                                    props['source'] = self.source
                                    props['source_url'] = self.source_url
                                    
                            yield exon_id, self.label, props
                    except:
                        print(
                            f'fail to process for label to load: {self.label}, type to load: {self.type}, data: {line}')

