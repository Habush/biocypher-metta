# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher_metta.adapters import Adapter
import csv
import gzip
from biocypher_metta.adapters.helpers import check_genomic_location, build_regulatory_region_id
from biocypher._logger import logger

#Example CADD Data
# rsid,chromosome,position,reference_allele,alternate_allele,raw_cadd_score,phred_score
# rs10,chr7,92383888,A,C,0.223125,6.177
# rs10,chr7,92383888,A,G,0.245354,6.489
# rs10,chr7,92383888,A,T,0.227741,6.243
# rs1000000,chr12,126890980,G,A,0.042237,3.295
# rs1000000,chr12,126890980,G,C,-0.017686,2.365

class CADDAdapter(Adapter):
    """
    Adapter for CADD data
    """
    def __init__(self, filepath, dbsnp_rsid_map,
                 write_properties, add_provenance,  
                 chr=None, start=None, end=None):
        self.file_path = filepath
        self.dbsnp_rsid_map = dbsnp_rsid_map
        self.chr = chr
        self.start = start
        self.end = end

        self.label = "sequence_variant"
        self.source = "CADD"
        self.source_url = "https://forgedb.cancer.gov/api/cadd/v1.0/cadd.forgedb.csv.gz"

        super(CADDAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.file_path, "rt") as fp:
            next(fp)
            reader = csv.reader(fp, delimiter=",")
            for row in reader:
                try:
                    rsid = row[0]
                    pos = self.dbsnp_rsid_map[rsid]["pos"]
                    chr = row[1]
                    ref = row[3]
                    alt = row[4]
                    _props = {}
                    if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                        if self.write_properties:
                            _props = {
                                'chr': chr,
                                'pos': pos,
                                'rsid': rsid,
                                'ref': ref,
                                'alt': alt,
                                'raw_cadd_score': float(row[5]),
                                'phred_score': float(row[6])
                            }
                            if self.add_provenance:
                                _props['source'] = self.source
                                _props['source_url'] = self.source_url
                        yield rsid, self.label, _props
                except KeyError as e:
                    logger.error(f"rsid {rsid} not found in dbsnp_rsid_map, skipping...")
                    continue

    def get_edges(self):
        pass