import os
import gzip
import json
import hashlib

from Bio.UniProt.GOA import gafiterator

from biocypher_metta.adapters import Adapter

# GAF files are defined here: https://geneontology.github.io/docs/go-annotation-file-gaf-format-2.2/
#
# Example:
# !gaf-version: 2.2
# !
# !generated-by: GOC
# !
# !date-generated: 2023-04-02T12:17
# !
# ...
# !
# !=================================
# !
# !Documentation about this header can be found here: https://github.com/geneontology/go-site/blob/master/docs/gaf_validation.md
# !
# UniProtKB	A0A024RBG1	NUDT4B	enables	GO:0003723	GO_REF:0000043	IEA	UniProtKB-KW:KW-0694	F	Diphosphoinositol polyphosphate phosphohydrolase NUDT4B	NUDT4B	protein	taxon:9606	20230306	UniProt
# UniProtKB	A0A024RBG1	NUDT4B	enables	GO:0046872	GO_REF:0000043	IEA	UniProtKB-KW:KW-0479	F	Diphosphoinositol polyphosphate phosphohydrolase NUDT4B	NUDT4B	protein	taxon:9606	20230306	UniProt
# UniProtKB	A0A024RBG1	NUDT4B	located_in	GO:0005829	GO_REF:0000052	IDA		C	Diphosphoinositol polyphosphate phosphohydrolase NUDT4B	NUDT4B	protein	taxon:9606	20161204	HPA
# UniProtKB	A0A075B6H7	IGKV3-7	involved_in	GO:0002250	GO_REF:0000043	IEA	UniProtKB-KW:KW-1064	P	Probable non-functional immunoglobulin kappa variable 3-7	IGKV3-7	protein	taxon:9606	20230306	UniProt
# UniProtKB	A0A075B6H7	IGKV3-7	located_in	GO:0005886	GO_REF:0000044	IEA	UniProtKB-SubCell:SL-0039	C	Probable non-functional immunoglobulin kappa variable 3-7	IGKV3-7	protein	taxon:9606	20230306	UniProt


# RNA Central file example:
#
# URS0000000055	ENSEMBL_GENCODE	ENST00000585414	9606	lncRNA	ENSG00000226803.9
# URS00000000C9	ENSEMBL_GENCODE	ENST00000514011	9606	lncRNA	ENSG00000248309.9
# URS00000000FD	ENSEMBL_GENCODE	ENST00000448543	9606	lncRNA	ENSG00000234279.2
# URS0000000351	ENSEMBL_GENCODE	ENST00000452009	9606	lncRNA	ENSG00000235427.1
# URS00000005D1	ENSEMBL_GENCODE	ENST00000563639	9606	lncRNA	ENSG00000260457.2
# URS0000000787	ENSEMBL_GENCODE	ENST00000452952	9606	lncRNA	ENSG00000206142.9
# URS0000000AA1	ENSEMBL_GENCODE	ENST00000615750	9606	lncRNA	ENSG00000277089.4
# URS0000000C0D	ENSEMBL_GENCODE	ENST00000582841	9606	lncRNA	ENSG00000265443.1
# URS0000000CF3	ENSEMBL_GENCODE	ENST00000414886	9606	lncRNA	ENSG00000226856.9

class GAFAdapter(Adapter):
    DATASET = 'gaf'
    RNACENTRAL_ID_MAPPING_PATH = './samples/rnacentral_ensembl_gencode.tsv.gz'
    SOURCES = {
        'human': 'http://geneontology.org/gene-associations/goa_human.gaf.gz',
        'human_isoform': 'http://geneontology.org/gene-associations/goa_human_isoform.gaf.gz',
        'rna': 'http://geneontology.org/gene-associations/goa_human_rna.gaf.gz',
        'rnacentral': 'https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/ensembl_gencode.tsv'
    }

    def __init__(self, filepath, write_properties, add_provenance, gaf_type='human'):
        if gaf_type not in GAFAdapter.SOURCES.keys():
            raise ValueError('Invalid type. Allowed values: ' +
                             ', '.join(GAFAdapter.SOURCES.keys()))

        self.filepath = filepath
        self.dataset = GAFAdapter.DATASET
        self.type = gaf_type
        self.label = "go_gene_product"
        self.source = "GO"
        self.source_url = GAFAdapter.SOURCES[gaf_type]

        super(GAFAdapter, self).__init__(write_properties, add_provenance)

    def load_rnacentral_mapping(self):
        self.rnacentral_mapping = {}
        with gzip.open(GAFAdapter.RNACENTRAL_ID_MAPPING_PATH, 'rt') as mapping_file:
            for annotation in mapping_file:
                mapping = annotation.split('\t')
                self.rnacentral_mapping[mapping[0] +
                                        '_' + mapping[3]] = mapping[2]

    def get_edges(self):

        if self.type == 'rna':
            self.load_rnacentral_mapping()

        with gzip.open(self.filepath, 'rt') as input_file:
            for annotation in gafiterator(input_file):
                source = annotation['GO_ID']
                target = annotation['DB_Object_ID']

                if self.type == 'rna':
                    transcript_id = self.rnacentral_mapping.get(
                        annotation['DB_Object_ID'])
                    if transcript_id is None:
                        continue
                    target = transcript_id
                props = {}
                if self.write_properties:
                    props = {
                        'qualifier': annotation['Qualifier'],
                        'db_reference': annotation['DB:Reference'],
                        'evidence': annotation['Evidence']
                    }
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield source, target, self.label, props

