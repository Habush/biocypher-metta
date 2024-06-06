import csv
import os
import pickle
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import to_float, check_genomic_location
from biocypher._logger import logger
import gzip

# Example QTEx eQTL input file:
# variant_id      gene_id tss_distance    ma_samples      ma_count        maf     pval_nominal    slope   slope_se        pval_nominal_threshold  min_pval_nominal        pval_beta
# chr1_845402_A_G_b38     ENSG00000225972.1       216340  4       4       0.0155039       2.89394e-06     2.04385 0.413032        2.775e-05       2.89394e-06     0.00337661
# chr1_920569_G_A_b38       ENSG00000225972.1       291507  4       4       0.0155039       1.07258e-05     1.92269 0.415516        2.775e-05       2.89394e-06     0.00337661

# description for column headers can be found here: 
# https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/README_eQTL_v8.txt

# Example *.egenes.txt.gz input files
# gene_id	gene_name	gene_chr	gene_start	gene_end	strand	num_var	beta_shape1	beta_shape2	true_df	pval_true_df	variant_id	tss_distance	chr	variant_pos	ref	alt	num_alt_per_site	rs_id_dbSNP151_GRCh38p7	minor_allele_samples	minor_allele_count	maf	ref_factor	pval_nominal	slope	slope_se	pval_perm	pval_beta	qval	pval_nominal_threshold	log2_aFC	log2_aFC_lower	log2_aFC_upper
# ENSG00000227232.5	WASH7P	chr1	14410	29553	-	1364	1.02984	294.487	455.958	6.29063e-08	chr1_64764_C_T_b38	35211	chr1	64764	C	T	1	rs769952832	70	71	0.0611015	1	1.01661e-08	0.586346	0.100677	9.999e-05	1.32112e-05	1.01141e-05	0.000505559	0.584194	0.435298	0.744545
# ENSG00000268903.1	RP11-34P13.15	chr1	135141	135895	-	1863	1.04872	330.017	441.174	0.00088883	chr1_103147_C_T_b38	-32748	chr1	103147	C	T	1	rs866355763	18	18	0.0154905	1	0.000347332	-0.612097	0.169958	0.241904	0.2337	0.0810742	0.000472534	-1.823931	-4.015491	-0.676333

class GTExEQTLAdapter(Adapter):
    # 1-based coordinate system

    def __init__(self, filepath, gtex_tissue_ontology_map,
                 write_properties, add_provenance, 
                 tissue_names=None, chr=None, start=None, end=None):
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
        self.gtex_tissue_ontology_map = pickle.load(open(gtex_tissue_ontology_map, 'rb'))
        self.tissue_names = tissue_names
        self.chr = chr
        self.start = start
        self.end = end
        self.label = 'gtex_variant_gene'
        self.source = 'GTEx'
        self.source_url = 'https://www.gtexportal.org/home/datasets'
        self.version = 'v8'


        super(GTExEQTLAdapter, self).__init__(write_properties, add_provenance)

    def get_edges(self):
        for file_name in os.listdir(self.filepath):
            if "egenes" in file_name: #skip other files
                tissue_name = file_name.split(".")[0]
                logger.info(f"Importing tissue: {tissue_name}")
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
                                    _props = {}
                                    if self.write_properties:
                                        _props = {
                                            #add re factor
                                            'maf': to_float(row[21]),
                                            'slope': to_float(row[24]),
                                            'p_value': to_float(row[27]),
                                            'q_value': to_float(row[28]),
                                            'biological_context': self.gtex_tissue_ontology_map[tissue_name]
                                        }
                                        if self.add_provenance:
                                            _props['source'] = self.source
                                            _props['source_url'] = self.source_url

                                    yield _source, _target, self.label, _props
                            except Exception as e:
                                print(row)
                                print(e)
