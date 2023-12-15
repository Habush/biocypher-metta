import csv
import gzip
import json
import os
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_variant_id, to_float, check_genomic_location


# Example TOPLD input data file:

# SNP1,SNP2,Uniq_ID_1,Uniq_ID_2,R2,Dprime,+/-corr
# 5031031,5032123,5031031:C:T,5032123:G:A,0.251,0.888,+
# 5031031,5063457,5031031:C:T,5063457:G:C,0.443,0.832,+

# Example TOPLD annotation file:

# Position,rsID,MAF,REF,ALT,Uniq_ID,VEP_ensembl_Gene_Name,VEP_ensembl_Consequence,CADD_phred,fathmm_XF_coding_or_noncoding,FANTOM5_enhancer_expressed_tissue_cell
# 5031031,rs1441313282,0.010486891385767793,C,T,5031031:C:T,FP565260.3|FP565260.3|FP565260.3|FP565260.3|FP565260.3,"intron_variant|intron_variant|intron_variant|intron_variant,NMD_transcript_variant|intron_variant",2.135,.,.


class TopLDAdapter(Adapter):

    def __init__(self, filepath, chr, ancestry='SAS'):
        self.file_path = filepath
        self.chr = chr
        self.ancestry = ancestry
        self.label = "in_ld_with"
        super(TopLDAdapter, self).__init__()

    def get_edges(self):
        with gzip.open(self.file_path, 'rt') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                try:
                    var1_pos = int(row[0])
                    var2_pos = int(row[1])
                    var1_ref, var1_alt = row[2].split(":")[1], row[2].split(":")[2]
                    var2_ref, var2_alt = row[3].split(":")[1], row[3].split(":")[2]
                    _source = build_variant_id(self.chr, var1_pos, var1_ref, var1_alt)
                    _target = build_variant_id(self.chr, var2_pos, var2_ref, var2_alt)

                    props = {
                        'r2': to_float(f"{row[6]}{row[4]}"),
                        'd_prime': float(row[5]),
                        'ancestry': self.ancestry,
                        'source': 'TopLD',
                        'source_url': 'http://topld.genetics.unc.edu/'
                    }

                    yield '', _source, _target, self.label, props

                except Exception as e:
                    print("Couldn't process row:", row)
                    print(e)
