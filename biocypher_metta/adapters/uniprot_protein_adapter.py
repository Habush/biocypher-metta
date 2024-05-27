import gzip
import json
import os
from biocypher_metta.adapters import Adapter
from Bio import SwissProt


# Data file is uniprot_sprot_human.dat.gz and uniprot_trembl_human.dat.gz at https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/.
# We can use SeqIO from Bio to read the file.
# Each record in file will have those attributes: https://biopython.org/docs/1.75/api/Bio.SeqRecord.html
# id, name will be loaded for protein. Ensembl IDs(example: Ensembl:ENST00000372839.7) in dbxrefs will be used to create protein and transcript relationship.


class UniprotProteinAdapter(Adapter):
   # ALLOWED_SOURCES = ['UniProtKB/Swiss-Prot', 'UniProtKB/TrEMBL']

    def __init__(self, filepath, write_properties, add_provenance):
        self.filepath = filepath
        self.dataset = 'UniProtKB_protein'
        self.label = 'protein'
        self.source = "Uniprot"
        self.source_url = "https://www.uniprot.org/"

        super(UniprotProteinAdapter, self).__init__(write_properties, add_provenance)
    
    def get_dbxrefs(self, cross_references):
        dbxrefs = []
        for cross_reference in cross_references:
            database_name = cross_reference[0]
            if database_name == 'EMBL':
                for item in cross_reference[1:3]:
                    if item != '-':
                        id = database_name + ':' + item
                        dbxrefs.append(id)
            elif database_name in ['RefSeq', 'Ensembl', 'MANE-Select']:
                for item in cross_reference[1:]:
                    if item != '-':
                        id = database_name + ':' + item.split('. ')[0]
                        dbxrefs.append(id)
            else:
                id = cross_reference[0] + ':' + cross_reference[1]
                dbxrefs.append(id)

        return sorted(list(set(dbxrefs)), key=str.casefold)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as input_file:
            records = SwissProt.parse(input_file)
            for record in records:
                dbxrefs = self.get_dbxrefs(record.cross_references)
                id = record.accessions[0]
                props = {}
                if self.write_properties:
                    props = {
                        'accessions': record.accessions[1:] if len(record.accessions) > 1 else record.accessions[0],
                        'protein_name': record.entry_name.split('_')[0],
                        'synonyms': dbxrefs
                    }
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url
                yield id, self.label, props
