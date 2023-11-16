import requests
from requests.adapters import HTTPAdapter, Retry
from requests.exceptions import JSONDecodeError
from tqdm import tqdm
from biocypher_metta.adapters import Adapter


# Example pathway input file:
# R-GGA-199992	trans-Golgi Network Vesicle Budding	Gallus gallus
# R-HSA-164843	2-LTR circle formation	Homo sapiens
# R-HSA-73843	5-Phosphoribose 1-diphosphate biosynthesis	Homo sapiens
# R-HSA-1971475	A tetrasaccharide linker sequence is required for GAG synthesis	Homo sapiens
# R-HSA-5619084	ABC transporter disorders	Homo sapiens


class ReactomePathwayAdapter(Adapter):

    def __init__(self, filepath=None, dry_run=False):

        self.filepath = filepath
        self.label = 'pathway'
        self.dataset = 'pathway'

    def get_nodes(self):
        session = requests.Session()
        retries = Retry(total=5, backoff_factor=1,
                        status_forcelist=[500, 502, 503, 504])
        session.mount('https://', HTTPAdapter(max_retries=retries))

        with open(self.filepath) as input:
            for line in input:
                id, name, species = line.strip().split('\t')
                if species == 'Homo sapiens':
                    # try:
                    #     query = 'https://reactome.org/ContentService/data/query/' + id #TODO this takes a long time
                    #     response = session.get(query)
                    #     data = response.json()
                    # except JSONDecodeError as e:
                    #     print(
                    #         f'Can not query for {query}. The status code is {response.status_code}. The text is {response.text}')
                    #     raise JSONDecodeError()
                    # id_version = data['stIdVersion']
                    # is_in_disease = data['isInDisease']
                    # name_aliases = data['name']
                    # is_top_level_pathway = False
                    # if data['className'] == 'TopLevelPathway':
                    #     is_top_level_pathway = True
                    props = {
                            'name': name,
                            # 'id_version': id_version,
                            # 'is_in_disease': is_in_disease,
                            # 'name_aliases': name_aliases,
                            # 'is_top_level_pathway': is_top_level_pathway
                            }
                    # if is_in_disease:
                    #     disease = data.get('disease')
                    #     disease_ontology_terms = []
                    #     for d in disease:
                    #         disease_ontology_term = 'ontology_terms/' + \
                    #             d['databaseName'] + '_' + d['identifier']
                    #         disease_ontology_terms.append(
                    #             disease_ontology_term)
                    #     props.update(
                    #         {'disease_ontology_terms': disease_ontology_terms}
                    #     )
                    # go_biological_process = data.get('goBiologicalProcess')
                    # if go_biological_process:
                    #     props.update(
                    #         {
                    #             'go_biological_process': go_biological_process['databaseName'] + '_' + go_biological_process['accession']
                    #         }
                    #     )
                    yield id, self.label, props