# Author Abdulrahman S. Omar <xabush@singularitynet.io>

from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import check_genomic_location, build_regulatory_region_id

## Example data:
# 1|chr1|800000|1350000,SAMD11|Ensembl:ENSG00000187634|HGNC:SAMD11;NOC2L|Ensembl:ENSG00000188976|HGNC:NOC2L;KLHL17|Ensembl:ENSG00000187961|HGNC:KLHL17;PLEKHN1|Ensembl:ENSG00000187583|HGNC:PLEKHN1;PERM1|Ensembl:ENSG00000187642|HGNC:PERM1;HES4|Ensembl:ENSG00000188290|HGNC:HES4;ISG15|Ensembl:ENSG00000187608|HGNC:ISG15;AGRN|Ensembl:ENSG00000188157|HGNC:AGRN;RNF223|Ensembl:ENSG00000237330|HGNC:RNF223;C1orf159|Ensembl:ENSG00000131591|HGNC:C1orf159;TTLL10|Ensembl:ENSG00000162571|HGNC:TTLL10;TNFRSF18|Ensembl:ENSG00000186891|HGNC:TNFRSF18;TNFRSF4|Ensembl:ENSG00000186827|HGNC:TNFRSF4;SDF4|Ensembl:ENSG00000078808|HGNC:SDF4;B3GALT6|Ensembl:ENSG00000176022|HGNC:B3GALT6;C1QTNF12|Ensembl:ENSG00000184163|HGNC:C1QTNF12;UBE2J2|Ensembl:ENSG00000160087|HGNC:UBE2J2;SCNN1D|Ensembl:ENSG00000162572|HGNC:SCNN1D;ACAP3|Ensembl:ENSG00000131584|HGNC:ACAP3;PUSL1|Ensembl:ENSG00000169972|HGNC:PUSL1;INTS11|Ensembl:ENSG00000127054|HGNC:INTS11;CPTP|Ensembl:ENSG00000224051|HGNC:CPTP;TAS1R3|Ensembl:ENSG00000169962|HGNC:TAS1R3;DVL1|Ensembl:ENSG00000107404|HGNC:DVL1

class TADMapAdapter(Adapter):
    """
    Adapter for Topologically Associated Domain (TAD) data.
    TAD are contiguous segments of the genome where the genomic elements are in frequent contact with each other.
    Source : TADMap https://cb.csail.mit.edu/cb/tadmap/
    """

    def __init__(self, filepath, type, chr=None, start=None, end=None):
        """
        :type filepath: str
        :type type: str
        :type chr: str
        :type start: int
        :type end: int
        :param filepath: path to the TAD file
        :param type: type of TAD data
        :param chr: chromosome name
        :param start: start position
        :param end: end position
        """
        self.filepath = filepath
        self.type = type
        self.dataset = 'tad'
        self.source = 'TADMap'
        self.source_url = 'https://cb.csail.mit.edu/cb/tadmap/'

        self.chr = chr
        self.start = start
        self.end = end

        if self.type == 'node':
            self.label = "tad"
        elif self.type == 'edge':
            self.label = "in_tad_region"

        super(TADMapAdapter, self).__init__()


    def get_nodes(self):
        """
        :return: generator of TAD nodes
        """
        with open(self.filepath, 'r') as tad_file:
            next(tad_file) # skip header
            for row in tad_file:
                row = row.strip().split(',')
                loc_info = row[0].split('|')
                genes_info = row[1].split(';')
                chr = loc_info[1]
                start = loc_info[2]
                end = loc_info[3]
                genes = []
                for gene in genes_info:
                    try:
                        gene = gene.split('|')
                        gene = gene[1].split(':')[1]
                        genes.append(gene)
                    except IndexError:
                        continue

                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    _id = build_regulatory_region_id(chr, start, end)
                    _props = {
                        'chr': chr,
                        'start': int(start),
                        'end': int(end),
                        'genes': genes
                    }

                    yield _id, self.label, _props

    def get_edges(self):
        """
        :return: generator of TAD edges
        """
        with open(self.filepath, 'r') as tad_file:
            next(tad_file)
            for row in tad_file:
                row = row.strip().split(',')
                loc_info = row[0].split('|')
                genes_info = row[1].split(';')
                chr = loc_info[1]
                start = loc_info[2]
                end = loc_info[3]
                genes = []
                for gene in genes_info:
                    try:
                        gene = gene.split('|')
                        gene = gene[1].split(':')[1]
                        genes.append(gene)
                    except IndexError:
                        continue

                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    _source = build_regulatory_region_id(chr, start, end)
                    for gene in genes:
                        _id = f"{self.label}-{_source}-{gene}"
                        yield _id, gene, _source, self.label, {}