import csv
import os
import pickle
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_variant_id, to_float, check_genomic_location
import gzip

# Example QTEx eQTL input file:
# variant_id      gene_id tss_distance    ma_samples      ma_count        maf     pval_nominal    slope   slope_se        pval_nominal_threshold  min_pval_nominal        pval_beta
# chr1_845402_A_G_b38     ENSG00000225972.1       216340  4       4       0.0155039       2.89394e-06     2.04385 0.413032        2.775e-05       2.89394e-06     0.00337661
# chr1_920569_G_A_b38       ENSG00000225972.1       291507  4       4       0.0155039       1.07258e-05     1.92269 0.415516        2.775e-05       2.89394e-06     0.00337661


class GTExEQTLAdapter(Adapter):
    # 1-based coordinate system

    def __init__(self, filepath, tissue_names=None, chr=None, start=None, end=None):
        """
        :type filepath: str
        :type tissue_names: str
        :type chr: str
        :type start: int
        :type end: int
        :param filepath: path to the directory containing eQTL data from GTEx
        :param tissue_names: tissue names to be used as biological context. If None, then all tissues are imported
        :param chr: chromosome name
        :param start: start position
        :param end: end position
        """
        self.filepath = filepath

        assert os.path.isdir(self.filepath), "The path to the directory containing eQTL data is not directory"
        self.filepath = filepath
        self.tissue_names = tissue_names
        self.chr = chr
        self.start = start
        self.end = end
        self.label = 'eqtl'
        self.source = 'GTEx'
        self.source_url = 'https://www.gtexportal.org/home/datasets'
        self.version = 'v8'


        super(GTExEQTLAdapter, self).__init__()

    def get_edges(self):

        for file_name in os.listdir(self.filepath):
            if "egenes" in file_name: #skip other files
                tissue_name = file_name.split(".")[0]
                if self.tissue_names is None or tissue_name in self.tissue_names:
                    with gzip.open(os.path.join(self.filepath, file_name), 'rt') as qtl:
                        next(qtl) # skip header
                        qtl_csv = csv.reader(qtl, delimiter='\t')
                        for row in qtl_csv:
                            try:
                                chr, pos, ref_seq, alt_seq, assembly_code = row[11].split('_')
                                pos = int(pos)
                                if assembly_code != 'b38':
                                    print('Unsuported assembly: ' + assembly_code)
                                    continue

                                variant_id = row[18]
                                if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                                    _source = variant_id
                                    _target = row[0].split('.')[0]
                                    _props = {
                                        #add re factor
                                        'maf': to_float(row[21]),
                                        'slope': to_float(row[24]),
                                        'p_value': to_float(row[27]),
                                        'q_value': to_float(row[28]),
                                        'biological_context': tissue_name
                                    }

                                    yield _source, _target, self.label, _props
                            except Exception as e:
                                print(row)
                                print(e)
