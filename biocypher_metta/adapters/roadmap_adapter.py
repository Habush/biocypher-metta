# Author Abdulrahman S. Omar <xabush@singularitynet.io>
import csv
import gzip
import os.path

from biocypher_metta.adapters import Adapter
import biocypher._logger as logger


class RoadMapAdapter(Adapter):

    def __init__(self, filepath):
        """
        :type filepath: str
        :param filepath: path to the directory containing epigenomic data
        """
        self.filepath = filepath

        assert os.path.isdir(self.filepath), "The path to the directory containing epigenomic data is not directory"

        self.source = 'Roadmap Epigenomics Project'
        self.source_url = ['https://forgedb.cancer.gov/api/forge2.erc2-chromatin15state-all/v1.0/forge2.erc2'
                           '-chromatin15state-all.{0-9}.forgedb.csv.gz',
                           "https://forgedb.cancer.gov/api/forge2.erc2-H3-all/v1.0/forge2.erc2-H3-all.{"
                           "0-9}.forgedb.csv.gz",
                           "https://forgedb.cancer.gov/api/forge2.erc2-DHS/v1.0/forge2.erc2-DHS.forgedb.csv.gz"] # {0-9} indicates this dataset is split into
        # 10 parts
        self.label = "regulatory_region"

        super(RoadMapAdapter, self).__init__()


    def get_nodes(self):

        for file_name in os.listdir(self.filepath):
            with gzip.open(os.path.join(self.filepath, file_name), "rt") as fp:
                next(fp)
                reader = csv.reader(fp, delimiter=',')
                for row in reader:
                    try:
                        _id = row[0]
                        _props = {
                            'biological_context': row[2],
                            'tissue': row[3],
                            'biochemical_activity': row[4]
                        }
                        yield _id, self.label, _props

                    except Exception as e:
                        logger.get_logger(f"error while parsing row: {row}, error: {e} skipping...")