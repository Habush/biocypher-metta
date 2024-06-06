import os
import pickle
import csv
from biocypher_metta.adapters import Adapter

# Example TF motif file from HOCOMOCO (e.g. ATF1_HUMAN.H11MO.0.B.pwm), which adastra used.
# Each pwm (position weight matrix) is a N x 4 matrix, where N is the length of the TF motif.
# >ATF1_HUMAN.H11MO.0.B
# 0.13987841615191168	0.3608848685361665	-0.4517428977281205	-0.25006525443490907
# 0.3553616298110878	-0.7707228823699511	0.3663777685862727	-0.4032800768485152
# 0.338606483241983	-0.8414816866844361	0.7066714674964639	-1.975404324094267
# -2.4841943814369576	-3.324736791140142	-1.9754043240942667	1.3195988707291637
# -2.3936647653320064	-2.9606438896517475	1.328010168941795	-2.4841943814369576
# 1.3217083373824634	-2.9606438896517475	-2.5837428525823363	-2.0963708769895044
# -1.9754043240942667	1.191348924768384	-1.4888815002078464	-1.0666727869454955
# -2.2340103118057417	-2.5837428525823363	1.2450934398879747	-1.0666727869454955
# -1.1380405061628662	0.23796182230991783	-2.5837428525823363	0.8481840389725303
# 0.13987841615191168	0.6170180100710398	-0.5426512454816383	-0.8788317538331962
# 0.7561011054759478	-0.7707228823699511	-0.2914989252431338	-0.4151773801942997


class HoCoMoCoMotifAdapter(Adapter):
    def __init__(self, filepath, annotation_file, hgnc_to_ensembl_map,
                 write_properties, add_provenance):

        self.filepath = filepath
        assert os.path.isdir(self.filepath), f"{self.filepath} is not a directory"
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.model_tf_path = annotation_file

        self.label = 'motif'
        self.source = 'HOCOMOCOv11'
        self.source_url = 'hocomoco11.autosome.org/motif/'

        self.load_model_tf_mapping()

        super(HoCoMoCoMotifAdapter, self).__init__(write_properties, add_provenance)

    def load_model_tf_mapping(self):
        self.model_tf_map = {}  # e.g. key: 'ANDR_HUMAN'; value: 'P10275'
        with open(self.model_tf_path, 'r') as fp:
            next(fp)
            reader = csv.reader(fp, delimiter='\t')
            for row in reader:
                model = row[0].strip()
                tf = row[1].strip()
                self.model_tf_map[model] = tf
    def get_nodes(self):
        for filename in os.listdir(self.filepath):
            if filename.endswith('.pwm'):
                model_name = filename.replace('.pwm', '')
                pwm = {"pmw_A": [], "pmw_C": [], "pmw_G": [], "pmw_T": []}
                with open(self.filepath + '/' + filename, 'r') as pwm_file:
                    next(pwm_file)
                    reader = csv.reader(pwm_file, delimiter='\t')
                    for row in reader:
                        pwm["pmw_A"].append(float(row[0]))
                        pwm["pmw_C"].append(float(row[1]))
                        pwm["pmw_G"].append(float(row[2]))
                        pwm["pmw_T"].append(float(row[3]))

                length = len(pwm["pmw_A"])

                tf_name = self.model_tf_map.get(model_name)
                _id = self.hgnc_to_ensembl_map.get(tf_name)
                if _id is None:
                    continue

                props = {}
                if self.write_properties:
                    props = {
                        'tf_name': tf_name,
                        'pwm_A': pwm["pmw_A"],
                        'pwm_C': pwm["pmw_C"],
                        'pwm_G': pwm["pmw_G"],
                        'pwm_T': pwm["pmw_T"],
                        'length': length
                    }
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield _id, self.label, props
