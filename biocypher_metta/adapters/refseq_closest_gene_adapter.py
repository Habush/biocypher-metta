# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher_metta.adapters import Adapter
import pickle
import csv
import gzip
from biocypher_metta.adapters.helpers import check_genomic_location, build_regulatory_region_id
from biocypher._logger import logger

#Example RefSeq Closest Gene Data
# rsid,chromosome,start_position,end_position,gene_chromosome,gene_start_position,gene_end_position,gene_symbol
# rs10,chr7,92383887,92383888,chr7,92244452,92462637,CDK6
# rs1000000,chr12,126890979,126890980,chr12,126918772,126918773,LINC02350
# rs10000003,chr4,57561646,57561647,chr4,57514884,57523982,HOPX
# rs10000005,chr4,85161557,85161558,chr4,85220321,85220322,LINC02994

class RefSeqClosestGeneAdapter(Adapter):
    """
    Adapter for RefSeq Closest Gene data
    """
    def __init__(self, filepath, hgnc_to_ensembl_map, dbsnp_rsid_map,
                 write_properties, add_provenance,
                 chr=None, start=None, end=None):
        self.file_path = filepath
        self.dbsnp_rsid_map = dbsnp_rsid_map
        self.chr = chr
        self.start = start
        self.end = end
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))

        self.label = "closest_gene"
        self.source = "RefSeq Closest Gene"
        self.source_url = "https://forgedb.cancer.gov/api/closest_gene/v1.0/closest_gene.forgedb.csv.gz"

        super(RefSeqClosestGeneAdapter, self).__init__(write_properties, add_provenance)

    def get_edges(self):
        with gzip.open(self.file_path, "rt") as fp:
            next(fp)
            reader = csv.reader(fp, delimiter=",")
            for row in reader:
                try:
                    rsid = row[0]
                    chr = row[1]
                    pos = self.dbsnp_rsid_map[rsid]["pos"]
                    if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                        try:
                            source_id = rsid
                            target_id = self.hgnc_to_ensembl_map[(row[7]).strip()]
                            distance = int(row[5]) + 1 - int(pos)
                            props = {}
                            if self.write_properties:
                                props['chr'] = chr
                                props['pos'] = pos
                                props['distance'] = distance
                                if self.add_provenance:
                                    props['source'] = self.source
                                    props['source_url'] = self.source_url

                            yield source_id, target_id, self.label, props

                        except Exception as e:
                            logger.error(f"error while parsing row: {row}, error: {e} skipping...")
                            continue
                except KeyError as e:
                    logger.error(f"rsid {rsid} not found in dbsnp_rsid_map, skipping...")
                    continue