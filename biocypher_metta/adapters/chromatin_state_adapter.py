# Author Abdulrahman S. Omar <xabush@singularitynet.io>
import csv
import gzip
import os.path

from biocypher_metta.adapters import Adapter
import pickle
from biocypher_metta.adapters.helpers import check_genomic_location, build_regulatory_region_id


class ChromatinStateAdapter(Adapter):

    def __init__(self, filepath, tissue_id_map, chr=None, start=None, end=None):
        """
        :type filepath: str
        :type tissue_id_map: str
        :type chr: str
        :type start: int
        :type end: int
        :param filepath: path to the directory containing epigenomic data
        :param tissue_id_map: path to the tissue id map (pickle file)
        :param chr: chromosome name
        :param start: start position
        :param end: end position
        """
        self.filepath = filepath

        assert os.path.isdir(self.filepath), "The path to the directory containing epigenomic data is not directory"

        self.source = 'Roadmap Epigenomics Project'
        self.source_url = 'https://personal.broadinstitute.org/jernst/MODEL_IMPUTED12MARKS/'

        self.chr = chr
        self.start = start
        self.end = end
        self.tissue_id_map = pickle.load(open(tissue_id_map, "rb"))

        self.label = "regulatory_region"

        super(ChromatinStateAdapter, self).__init__()


    def get_nodes(self):

        for file_name in os.listdir(self.filepath):
            tissue_id = file_name.split("_")[0]
            tissue_name = self.tissue_id_map[tissue_id]
            with gzip.open(os.path.join(self.filepath, file_name), "rt") as fp:
                reader = csv.reader(fp, delimiter='\t')
                for row in reader:
                    try:
                        if check_genomic_location(self.chr, self.start, self.end, row[0], row[1], row[2]):
                            _id = build_regulatory_region_id(row[0], row[1], row[2])
                            _props = {
                                'chr': row[0],
                                'start': int(row[1]),
                                'end': int(row[2]),
                                'tissue': tissue_name,
                                'biochemical_activity': row[3].split("_")[1],
                                'source': self.source,
                                'source_url': self.source_url
                            }
                            yield _id, self.label, _props

                    except:
                        print(f'fail to process: {row}')