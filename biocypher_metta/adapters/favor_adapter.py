from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_variant_id, to_float, check_genomic_location
import json
import os
import csv
import gzip

# FIELDS = ['chromosome', 'start_position',
#      'ref_vcf', 'alt_vcf', 'aloft_value', 'aloft_description',
#     'apc_conservation', 'apc_conservation_v2', 'apc_epigenetics_active', 'apc_epigenetics',
#     'apc_epigenetics_repressed', 'apc_epigenetics_transcription', 'apc_local_nucleotide_diversity',
#     'apc_local_nucleotide_diversity_v2', 'apc_local_nucleotide_diversity_v3', 'apc_mappability', 'apc_micro_rna',
#     'apc_mutation_density', 'apc_protein_function', 'apc_protein_function_v2', 'apc_protein_function_v3',
#     'apc_proximity_to_coding', 'apc_proximity_to_coding_v2', 'apc_proximity_to_tsstes', 'apc_transcription_factor',
#     'bravo_an', 'bravo_af', 'filter_status', 'clnsig', 'clnsigincl', 'clndn', 'clndnincl', 'clnrevstat', 'origin',
#     'clndisdb', 'clndisdbincl', 'geneinfo', 'polyphen2_hdiv_score', 'polyphen2_hvar_score', 'mutation_taster_score',
#     'mutation_assessor_score', 'metasvm_pred', 'fathmm_xf', 'funseq_value', 'funseq_description',
#     'genecode_comprehensive_categoty', 'af_total', 'af_asj_female', 'af_eas_female', 'af_afr_male', 'af_female',
#     'af_fin_male', 'af_oth_female', 'af_ami', 'af_oth', 'af_male', 'af_ami_female', 'af_afr', 'af_eas_male', 'af_sas',
#     'af_nfe_female', 'af_asj_male', 'af_raw', 'af_oth_male', 'af_nfe_male', 'af_asj', 'af_amr_male', 'af_amr_female',
#     'af_amr_sas_female', 'af_fin', 'af_afr_female', 'af_sas_male', 'af_amr', 'af_nfe', 'af_eas', 'af_ami_male',
#     'af_fin_female', 'sift_cat', 'sift_val', 'polyphen_cat', 'polyphen_val', 'cadd_rawscore', 'cadd_phred',
#     'refseq_category', 'tg_afr', 'tg_all', 'tg_amr', 'tg_eas', 'tg_eur', 'tg_sas'
# ]


FIELDS = {'chromosome': 3, 'start_position': 4, 'ref_vcf': 9, 'alt_vcf': 10, 'aloft_value': 11, 'aloft_description': 12,
          'apc_conservation': 13, 'apc_conservation_v2': 14,
          'apc_epigenetics_active': 15, 'apc_epigenetics': 16, 'apc_epigenetics_repressed': 17,
          'apc_epigenetics_transcription': 18,
          'apc_local_nucleotide_diversity': 19, 'apc_local_nucleotide_diversity_v2': 20,
          'apc_local_nucleotide_diversity_v3': 21,
          'apc_mappability': 22, 'apc_micro_rna': 23, 'apc_mutation_density': 24, 'apc_protein_function': 25,
          'apc_protein_function_v2': 26,
          'apc_protein_function_v3': 27, 'apc_proximity_to_coding': 28, 'apc_proximity_to_coding_v2': 29,
          'apc_proximity_to_tsstes': 30,
          'apc_transcription_factor': 31, 'bravo_an': 32, 'bravo_af': 33, 'filter_status': 34, 'clnsig': 38,
          'clnsigincl': 39, 'clndn': 40,
          'clndnincl': 41, 'clnrevstat': 42, 'origin': 43, 'clndisdb': 44, 'clndisdbincl': 45, 'geneinfo': 46,
          'polyphen2_hdiv_score': 47,
          'polyphen2_hvar_score': 48, 'mutation_taster_score': 49, 'mutation_assessor_score': 50, 'metasvm_pred': 51,
          'fathmm_xf': 52,
          'funseq_value': 53, 'funseq_description': 54, 'af_total': 60, 'af_asj_female': 61, 'af_eas_female': 62,
          'af_afr_male': 63,
          'af_female': 64, 'af_fin_male': 65, 'af_oth_female': 66, 'af_ami': 67, 'af_oth': 68, 'af_male': 69,
          'af_ami_female': 70,
          'af_afr': 71, 'af_eas_male': 72, 'af_sas': 73, 'af_nfe_female': 74, 'af_asj_male': 75, 'af_raw': 76,
          'af_oth_male': 77,
          'af_nfe_male': 78, 'af_asj': 79, 'af_amr_male': 80, 'af_amr_female': 81, 'af_fin': 83, 'af_afr_female': 84,
          'af_sas_male': 85,
          'af_amr': 86, 'af_nfe': 87, 'af_eas': 88, 'af_ami_male': 89, 'af_fin_female': 90, 'sift_cat': 96,
          'sift_val': 97, 'polyphen_cat': 98,
          'polyphen_val': 99, 'cadd_rawscore': 161, 'cadd_phred': 162, 'refseq_category': 174, 'tg_afr': 179,
          'tg_all': 180, 'tg_amr': 181, 'tg_eas': 182, 'tg_eur': 183, 'tg_sas': 184}

class FavorAdapter(Adapter):
    # Originally 1-based coordinate system
    # Converted to 0-based

    WRITE_THRESHOLD = 1000000

    def __init__(self, write_properties, add_provenance, 
                 filepath=None, chr=None, start=None, end=None):
        self.filepath = filepath
        self.chr = chr
        self.start = start
        self.end = end
        self.label = "sequence_variant"
        self.source = "FAVOR"
        self.source_url = "http://favor.genohub.org/"

        super(FavorAdapter, self).__init__(write_properties, add_provenance)

    def convert_freq_value(self, value):
        if value == '.':
            value = 0

        try:
            value = to_float(value)
        except:
            pass

        return value

    # only selecting FREQ value from INFO data
    def parse_annotation(self, row):
        annotations = {}
        exclude = ["chromosome", "start_position", "ref_vcf", "alt_vcf"]

        for k, v in FIELDS.items():
            if k not in exclude:
                annotations[k] = self.convert_freq_value(row[v])

        return annotations

    def get_nodes(self):

        with open(self.filepath, 'r') as f:
            next(f)
            reader = csv.reader(f, delimiter=',')

            for row in reader:

                chr = "chr" + row[FIELDS["chromosome"]]
                pos = int(row[FIELDS["start_position"]])

                if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                    id = build_variant_id(
                        chr, pos,
                        row[FIELDS["ref_vcf"]],
                        row[FIELDS["alt_vcf"]])
                    props = {}
                    if self.write_properties:
                        props = {
                            # '_key': id,
                            'chr': chr,
                            'start': pos,
                            'end': pos,
                            # 'rsid': [row[FIELDS["rsid"], #TODO uncomment when rsid is available
                            'ref': row[FIELDS["ref_vcf"]],
                            'alt': row[FIELDS["alt_vcf"]],
                            'annotation': self.parse_annotation(row),
                        }
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    # TODO add a simple heuristics to resolve conflicting rsids appear close to each other in data
                    #  files when the data becomes available

                    yield id, self.label, props

