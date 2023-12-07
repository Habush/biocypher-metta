import gzip
import csv
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_regulatory_region_id, check_genomic_location
# cCRE,all input file has 10 columns: chromsome, start, end, ID, score (all 0), strand (NA), start, end, color, biochemical_activity
# There are 8 types of biochemical_activity:
# pELS - proximal Enhancer-ike signal
# CA → chromatin accessible
# dELS - distal Enhancer-like signal
# TF → TF binding
# CA-CTCF
# CA-TF
# CA-H3K4me3
# PLS

# Below is example data:
# chr1    10033   10250   EH38E2776516    0       .       10033   10250   255,167,0       pELS
# chr1    10385   10713   EH38E2776517    0       .       10385   10713   255,167,0       pELS
# chr1    16097   16381   EH38E3951272    0       .       16097   16381   0,176,240       CA-CTCF
# chr1    17343   17642   EH38E3951273    0       .       17343   17642   190,40,229      CA-TF
# chr1    29320   29517   EH38E3951274    0       .       29320   29517   6,218,147       CA


class CCREAdapter(Adapter):
    DATASET = 'regulatory_region'

    BIOCHEMICAL_DESCRIPTION = {
        'pELS': 'proximal Enhancer-like signal',
        'CA': 'chromatin accessible',
        'dELS': 'distal Enhancer-like signal',
        'TF': 'TF binding',
        'CA-CTCF': 'chromatin accessible + CTCF binding',
        'CA-TF': 'chromatin accessible + TF binding',
        'CA-H3K4me3': 'chromatin accessible + H3K4me3 high signal',
        'PLS': 'Promoter-like signal'
    }

    def __init__(self, filepath, chr=None, start=None, end=None):
        """
        :type filepath: str
        :type chr: str
        :type start: int
        :type end: int
        :param filepath: path to the cCRE file
        :param chr: chromosome name
        :param start: start position
        :param end: end position
        """
        self.filepath = filepath
        self.dataset = CCREAdapter.DATASET
        self.source = 'ENCODE_SCREEN (ccREs)'
        self.source_url = 'https://www.encodeproject.org/files/ENCFF420VPZ/'

        self.chr = chr
        self.start = start
        self.end = end

        super(CCREAdapter, self).__init__()

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            label = 'regulatory_region'

            for row in reader:
                try:
                    if check_genomic_location(self.chr, self.start, self.end, row[0], row[1], row[2]):
                        #description = CCREAdapter.BIOCHEMICAL_DESCRIPTION.get(row[9])
                        # _id = row[3] #TODO use the position as in id
                        _id = build_regulatory_region_id(row[0], row[1], row[2])
                        _props = {
                            'chr': row[0],
                            'start': row[1],
                            'end': row[2],
                            'biochemical_activity': str(row[9]),
                            #'biochemical_activity_description': description,
                        }
                        yield _id, label, _props

                except:
                    print(f'fail to process: {row}')
                    pass
