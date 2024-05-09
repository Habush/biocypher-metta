# Author Abdulrahman S. Omar <xabush@singularitynet.io>
import csv
import gzip
import os.path

from biocypher_metta.adapters import Adapter
import biocypher._logger as logger

class Forge2TFAdapter(Adapter):

    def __init__(self, filepath, write_properties, add_provenance):
        """
        :type filepath: str
        :param filepath: path to the directory containing epigenomic data
        """
        self.filepath = filepath

        assert os.path.isdir(self.filepath), "The path to the directory containing epigenomic data is not directory"

        self.source = 'Forge2'
        self.source_url = ['https://forgedb.cancer.gov/api/forge2.erc2-chromatin15state-all/v1.0/forge2.erc2-chromatin15state-all.{0-9}.forgedb.csv.gz',
                        "https://forgedb.cancer.gov/api/forge2.erc2-H3-all/v1.0/forge2.erc2-H3-all.{0-9}.forgedb.csv.gz",
                        "https://forgedb.cancer.gov/api/forge2.erc2-DHS/v1.0/forge2.erc2-DHS.forgedb.csv.gz"] # {0-9} indicates this dataset is split into
        # 10 parts
        self.label = "regulatory_region"

        super(Forge2TFAdapter, self).__init__(write_properties, add_provenance)

