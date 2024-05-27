# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher_metta.adapters import Adapter
import pickle
import csv
import gzip
from biocypher_metta.adapters.helpers import check_genomic_location
from biocypher._logger import logger

#Example ABC Data
# rsid,chromosome,start_position,end_position,bait_chromosome,bait_start_position,bait_end_position,name,class,activity_base,target_gene,target_gene_tss,target_gene_expression,target_gene_promoter_activity_quantile,target_gene_is_expressed,distance,is_self_promoter,powerlaw_contact,powerlaw_contact_reference,hic_contact,hic_contact_pl_scaled,hic_contact_pl_scaled_adj,hic_pseudocount,abc_score_numerator,abc_score,powerlaw_score_numerator,powerlaw_score,cell_type,base_overlap
# rs10000009,chr4,71048952,71048953,chr4,71048249,71049305,intergenic|chr4:71048099-71049455,intergenic,18.942371,SULT1E1,70725870,NA,0.459997,True,322907,False,0.003269,0.003286,0.001376,0.001383,0.002599,0.001216,0.049232,0.043128,0.062253,0.024298,HepG2-Roadmap,1
# rs10000009,chr4,71048952,71048953,chr4,71048249,71049305,intergenic|chr4:71048099-71049455,intergenic,18.942371,UGT2A3,69817509,NA,0.766091,True,1231268,False,0.001013,0.001026,0.000301,0.000305,0.001318,0.001013,0.024967,0.015579,0.019429,0.004419,HepG2-Roadmap,1
# rs10000009,chr4,71048952,71048953,chr4,71048249,71049305,intergenic|chr4:71048099-71049455,intergenic,18.942371,ODAM,71062243,NA,0.410746,True,13466,False,0.052694,0.052138,0.00758,0.0075,0.008716,0.001216,0.165096,0.136843,0.987621,0.307587,HepG2-Roadmap,1

COL_DICT = {'rsid': 0, 'chromosome': 1, 'start_position': 2, 'end_position': 3, 'bait_chromosome': 4, 'bait_start_position': 5, 'bait_end_position': 6, 'name': 7, 'class': 8, 'activity_base': 9, 'target_gene': 10, 'target_gene_tss': 11, 'target_gene_expression': 12, 'target_gene_promoter_activity_quantile': 13, 'target_gene_is_expressed': 14, 'distance': 15, 'is_self_promoter': 16, 'powerlaw_contact': 17, 'powerlaw_contact_reference': 18, 'hic_contact': 19, 'hic_contact_pl_scaled': 20, 'hic_contact_pl_scaled_adj': 21, 'hic_pseudocount': 22, 'abc_score_numerator': 23, 'abc_score': 24, 'powerlaw_score_numerator': 25, 'powerlaw_score': 26, 'cell_type': 27, 'base_overlap': 28}

class ABCAdapter(Adapter):
    """
    Adapter for Activity-By-Contact (ABC) data from Fulco CP et.al 2019
    """
    def __init__(self, filepath, type, hgnc_to_ensembl_map, tissue_to_ontology_id_map,
                 dbsnp_rsid_map, write_properties, add_provenance,
                 chr=None, start=None, end=None):
        self.file_path = filepath
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.tissue_to_ontology_id_map = pickle.load(open(tissue_to_ontology_id_map, 'rb'))
        self.dbsnp_rsid_map = dbsnp_rsid_map
        self.chr = chr
        self.start = start
        self.end = end

        assert type in ["node", "edge"], f"type parameter should be either node or edge, got {type}"

        if type == "node":
            self.label = "regulatory_region"
        else:
            self.label = "regulatory_region_gene"
        self.source = "ABC"
        self.source_url = "https://forgedb.cancer.gov/api/abc/v1.0/abc.forgedb.csv.gz"
        super(ABCAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.file_path, "rt") as fp:
            next(fp)
            reader = csv.reader(fp, delimiter=",")
            for row in reader:
                try:
                    rsid = row[COL_DICT['rsid']]
                    chr = row[COL_DICT['chromosome']]
                    pos = self.dbsnp_rsid_map[rsid]["pos"]
                    if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                        _props = {
                            'chr': chr,
                            'start': pos,
                            'end': pos,
                            'biochemical_activity': 'DNase I hypersensitive',
                            'biological_context': self.tissue_to_ontology_id_map[row[COL_DICT['cell_type']]]
                        }
                        yield rsid, self.label, _props
                except KeyError as e:
                    logger.error(f"rsid {rsid} not found in dbsnp_rsid_map, skipping...")
                    continue



    def get_edges(self):
        with gzip.open(self.file_path, "rt") as fp:
            next(fp)
            reader = csv.reader(fp, delimiter=",")
            for row in reader:
                try:
                    rsid = row[COL_DICT['rsid']]
                    chr = row[COL_DICT['chromosome']]
                    pos = self.dbsnp_rsid_map[rsid]
                    if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                        try:
                            _source = rsid
                            _target = self.hgnc_to_ensembl_map[(row[COL_DICT['target_gene']]).strip()]
                            props = {
                                "abc_score": row[COL_DICT['abc_score']],
                                "biological_context": self.tissue_to_ontology_id_map[row[COL_DICT['cell_type']]]
                            }

                            yield _source, _target, self.label, props
                        except Exception as e:
                            print(f"error while parsing row: {row}, error: {e} skipping...")
                            continue
                except KeyError as e:
                    logger.error(f"rsid {rsid} not found in dbsnp_rsid_map, skipping...")
                    continue
