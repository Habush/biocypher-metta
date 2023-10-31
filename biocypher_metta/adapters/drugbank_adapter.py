from __future__ import annotations

from pypath.share import curl, settings

from pypath.inputs import drugbank, unichem

from contextlib import ExitStack
from typing import Union
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time

from biocypher._logger import logger

from enum import Enum, auto

logger.debug(f"Loading module {__name__}.")


class DrugBankNodeField(Enum):
    # primary fields
    DRUGBANK_ID = "drugbank_id"
    NAME = "name"
    CAS_NUMBER = "cas_number"
    GROUPS = "groups"
    GENERAL_REFERENCES = "general_references"
    ATC_CODES = "atc_codes"
    
    # property fields
    INCHI = "InChI"
    INCHIKEY = "InChIKey"
    
    # drugbank external fields
    KEGG_DRUG = "KEGG Drug"
    RXCUI = "RxCUI"
    PHARMGKB = "PharmGKB"
    PDB = "PDB"
    
    # unichem external mappings
    ZINC = "zinc"
    CHEMBL = "chembl"
    BINDINGDB = "bindingdb"
    CLINICALTRIALS = "clinicaltrials"
    CHEBI = "chebi"
    PUBCHEM = "pubchem"
    DRUGCENTRAL = "drugcentral"
    
    @classmethod
    def get_primary_fields(cls):
        """
        Get drugbank primary fields
        """
        return [cls.NAME, cls.CAS_NUMBER, cls.GROUPS, cls.GENERAL_REFERENCES, cls.ATC_CODES, cls.DRUGBANK_ID]
    
    @classmethod
    def get_property_fields(cls):
        """
        Get drugbank property fields
        """
        return [cls.INCHI, cls.INCHIKEY]
    
    @classmethod
    def get_drugbank_external_fields(cls):
        """
        Get drugbank external database cross-references
        """
        return [cls.KEGG_DRUG, cls.RXCUI, cls.PHARMGKB, cls.PDB]
    
    @classmethod
    def get_unichem_external_mappings(cls):
        """
        Get unichem database cross-references
        """
        return [cls.ZINC, cls.CHEMBL, cls.BINDINGDB, cls.CLINICALTRIALS, cls.CHEBI, cls.PUBCHEM, cls.DRUGCENTRAL]
    
class PrimaryNodeIdentifier(Enum):
    DRUGBANK = "drugbank"
    KEGG_DRUG = "KEGG Drug"
    KEGG_LIGAND = "kegg_ligand" # kegg ligand id
    DRUGCENTRAL = "drugcentral"
    CHEMBL = "chembl"
    RXCUI = "RxCUI"
    
class DrugbankDTIEdgeField(Enum):
    ACTIONS = "actions"
    REFERENCES = "references"
    KNOWN_ACTION = "known_action"
    

class DrugbankEdgeType(Enum):
    DRUG_TARGET_INTERACTION = auto()
    DRUG_DRUG_INTERACTION = auto()
    

class DrugBank:
    """
    Class that downloads DRUGBANK drug data using pypath and reformats it to be ready
    for import into BioCypher.
    """
    
    def __init__(self, drugbank_user:str, drugbank_passwd:str, node_fields:Union[list[DrugBankNodeField], None] = None,
                dti_edge_fields: Union[list[DrugbankDTIEdgeField], None] = None, primary_node_id:Union[PrimaryNodeIdentifier, None] = None,
                edge_types: Union[list[DrugbankEdgeType], None] = None,  add_prefix = True, test_mode: bool = False,):
        """
        Args
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
            node_fields: node properties that will be included in graph, if None select all fields belonging to DrugBankNodeField class
            dti_edge_fields: Drug_target interaction edge properties that will included in graph, if None select all fields belonging to DrugbankDTIEdgeField class
            primary_node_id: primary drug node identifier, if None selects drugbank as primary identifier
            edge_types: edge types that will be included in graph, if None select all fields belonging to DrugbankEdgeType class
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of data for testing
        """
        
        self.user = drugbank_user
        self.passwd = drugbank_passwd
        self.add_prefix = add_prefix
        
        # set node fields
        self.set_node_and_edge_fields(node_fields=node_fields, dti_edge_fields=dti_edge_fields)
        
        # set primary id of drug nodes
        self.set_primary_id(primary_node_id=primary_node_id)
        
        # set edge types
        self.set_edge_types(edge_types = edge_types)
        
        # set early_stopping, if test_mode true
        self.early_stopping = None
        if test_mode:
            self.early_stopping = 100
        
    def download_drugbank_data(self, cache=False, debug=False, retries=6,):
        """
        Wrapper function to download drug data from various databases using pypath.

        Args
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())
                
            
            t0 = time()
            
            # download node data
            self.download_drugbank_node_data()
            
            # create drugbank to primary id mapping
            self.set_drugbank_to_primary_id_mapping()
            
            if DrugbankEdgeType.DRUG_TARGET_INTERACTION in self.edge_types:
                # download DTI data
                self.download_drugbank_dti_data()
            
            if DrugbankEdgeType.DRUG_DRUG_INTERACTION in self.edge_types:
                # download DDI data
                self.download_drugbank_ddi_data()
            
            t1 = time()
            logger.info(f'All data is downloaded in {round((t1-t0) / 60, 2)} mins'.upper())
            
            
    def download_drugbank_node_data(self) -> None:
        """
        Wrapper function to download DrugBank drug entries using pypath
        """
        
        logger.debug('Downloading Drugbank drug node data')
        t0 = time()
        
        self.drugbank_data = drugbank.DrugbankFull(user = self.user, passwd = self.passwd)
        
        self.drugbank_drugs_external_ids = self.drugbank_data.drugbank_external_ids_full()
        
        if self.property_node_fields:            
            self.drugbank_properties = self.drugbank_data.drugbank_properties_full()
        
        if self.primary_node_fields:
            fields = self.primary_node_fields.copy()
            self.drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(fields = fields)
            
        if self.unichem_mapping_node_fields:
            self.get_unichem_mappings()
        
        t1 = time()
        logger.info(f'Drugbank drug node data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def download_drugbank_dti_data(self) -> None:
        """
        Wrapper function to download DrugBank DTI data using pypath
        """ 
        
        logger.debug('Downloading Drugbank DTI data')
        t0 = time()
        
        fields = self.dti_edge_fields.copy() + ["drugbank_id", "polypeptide"]
        self.drugbank_dti = self.drugbank_data.drugbank_targets_full(fields = fields)
        
        t1 = time()
        logger.info(f'Drugbank DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def download_drugbank_ddi_data(self) -> None:
        """
        Wrapper function to download DrugBank DDI data using pypath
        """ 
            
        logger.debug('Downloading Drugbank DDI data')
        t0 = time()
        
        self.drugbank_ddi = self.drugbank_data.drugbank_drugs_full(fields = "drug_interactions")
        
        t1 = time()
        logger.info(f'Drugbank DDI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def set_node_and_edge_fields(self, node_fields, dti_edge_fields):
        """
        Sets node and edge properties
        """
        if node_fields:
            self.all_node_fields = [field.value for field in node_fields]
            self.primary_node_fields = [field.value for field in DrugBankNodeField.get_primary_fields() if field in node_fields]
            self.property_node_fields = [field.value for field in DrugBankNodeField.get_property_fields() if field in node_fields]
            self.drugbank_external_node_fields = [field.value for field in DrugBankNodeField.get_drugbank_external_fields() if field in node_fields]
            self.unichem_mapping_node_fields = [field.value for field in DrugBankNodeField.get_unichem_external_mappings() if field in node_fields]
        else:
            self.all_node_fields = [field.value for field in DrugBankNodeField]
            self.primary_node_fields = [field.value for field in DrugBankNodeField.get_primary_fields()]
            self.property_node_fields = [field.value for field in DrugBankNodeField.get_property_fields()]
            self.drugbank_external_node_fields = [field.value for field in DrugBankNodeField.get_drugbank_external_fields()]
            self.unichem_mapping_node_fields = [field.value for field in DrugBankNodeField.get_unichem_external_mappings()]
            
        if dti_edge_fields:
            self.dti_edge_fields = [field.value for field in dti_edge_fields]
        else:
            self.dti_edge_fields = [field.value for field in DrugbankDTIEdgeField]
            
    def set_edge_types(self, edge_types):
        """
        Set edge types
        """
        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [field for field in DrugbankEdgeType]
        
    def get_unichem_mappings(self) -> None:
        """
        Download database mappings from pypath's unichem script and save them in a dictionary
        """
        
        self.unichem_mappings = {}
        
        for field in self.unichem_mapping_node_fields: 
            try:
                _dict = unichem.unichem_mapping("drugbank", field)
                
                # remove set data type from value
                _dict = {k:list(v)[0] for k, v in _dict.items()}
                
                # add unichem_mappings
                self.unichem_mappings[field] = _dict
            
            except TypeError: # in case id1-id2 doesnt exist, try id2-id1
                _dict = unichem.unichem_mapping(field, "drugbank")
                
                # switch values and keys
                _dict = {list(v)[0]:k for k, v in _dict.items()}
                
                # add unichem_mappings
                self.unichem_mappings[field] = _dict
                
    def set_primary_id(self, primary_node_id) -> None:
        """
        Set primary identifier of drug nodes
        """
        if primary_node_id:
            self.primary_id = primary_node_id.value
        else:
            self.primary_id = PrimaryNodeIdentifier.DRUGBANK.value
                            
                
    def set_drugbank_to_primary_id_mapping(self) -> None:
        """
        If identifiers other than drugbank are selected, create a mapping dictionary from drugbank to selected identifier
        """
        unichem_fields = [field.value for field in DrugBankNodeField.get_unichem_external_mappings()]
        drugbank_external_fields = [field.value for field in DrugBankNodeField.get_drugbank_external_fields()]
        
        if self.primary_id in unichem_fields:
            try:
                self.drugbank_to_primary_id = unichem.unichem_mapping("drugbank", self.primary_id)
                self.drugbank_to_primary_id = {k:list(v)[0] for k, v in drugbank_to_primary_id.items()}
                
            except TypeError:
                self.drugbank_to_primary_id = unichem.unichem_mapping(self.primary_id, "drugbank")
                self.drugbank_to_primary_id = {list(v)[0]:k for k, v in self.drugbank_to_primary_id.items()}
                
        elif self.primary_id in drugbank_external_fields:
            self.drugbank_to_primary_id = {}
            
            for drugbank_id, mapping_dict in self.drugbank_drugs_external_ids.items():
                if mapping_dict.get(self.primary_id, None):
                    self.drugbank_to_primary_id[drugbank_id] = mapping_dict[self.primary_id]
                    
        else:
            self.drugbank_to_primary_id = {}
            
    def add_prefix_to_id(self, prefix, identifier, sep=":") -> str:
        """
        Adds prefix to ids
        """
        if self.add_prefix:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
                    
    def get_dti_edges(self, label="drug_targets_protein") -> None:
        """
        Get Drug-target interaction edges
        Args:
            label = label of edges
        """
        
        self.dti_edge_list = []
        
        def get_uniprot_ids(element):
            """
            Extracts Swiss-Prot Uniprot ids
            """
            if element:
                if isinstance(element, tuple) and element[1] == 'Swiss-Prot':
                    return [element[0]]
                elif isinstance(element, list):
                    return [x[0] for x in element if x[1] == "Swiss-Prot"]
            else:
                return None
            
            
        logger.debug('Getting Drugbank DTI edges')
        
        counter = 0
        
        for dti in tqdm(self.drugbank_dti):
            uniprot_ids = get_uniprot_ids(dti.polypeptide)
            
            # if uniprot or drugbank id doesnt exist, skip it
            if not uniprot_ids or not dti.drugbank_id:
                continue
            
            # if primary id is not drugbank and it doesnt have any mapping skip it
            if self.primary_id != "drugbank" and self.drugbank_to_primary_id and not self.drugbank_to_primary_id.get(dti.drugbank_id, None):
                continue
                
            if self.primary_id != "drugbank":
                # if primary id is not drugbank, get mapping for it
                source = self.add_prefix_to_id("something", self.drugbank_to_primary_id[dti.drugbank_id])            
            else:
                source = self.add_prefix_to_id("drugbank", dti.drugbank_id)
                    
            for _id in uniprot_ids:
                # create props dict
                props = {str(k).replace(" ", "_").lower():v for k, v in dti._asdict().items() if k in self.dti_edge_fields}
                
                target = self.add_prefix_to_id("uniprot", _id)
                
                self.dti_edge_list.append((None, source, target, label, props))
                
                counter += 1
            
            if self.early_stopping and counter >= self.early_stopping:
                break
            
    def get_ddi_edges(self, label="drug_interacts_with_drug") -> None:
        """
        Get Drug-Drug interaction edges
        Args:
            label = label of edges
        """
        self.ddi_edge_list = []
        
        counter = 0
        
        logger.debug('Getting Drugbank DDI edges')
        
        for ddi in tqdm(self.drugbank_ddi):
            
            # if drug1 doesnt have any interaction, skip it
            if not ddi.drug_interactions:
                continue
                
            # if drug1's primary id is not drugbank and it doesnt have any mapping skip it
            if self.primary_id != "drugbank" and self.drugbank_to_primary_id and not self.drugbank_to_primary_id.get(ddi.drugbank_id, None):
                continue
                
            if self.primary_id != "drugbank":
                # if drug1's primary id is not drugbank, get mapping for it
                source = self.add_prefix_to_id("something", self.drugbank_to_primary_id[ddi.drugbank_id])             
            else:
                source = self.add_prefix_to_id("drugbank", ddi.drugbank_id)
            
            for pair2 in ddi.drug_interactions:
                
                # if drug2's primary id is not drugbank and it doesnt have any mapping skip it
                if self.primary_id != "drugbank" and self.drugbank_to_primary_id and not self.drugbank_to_primary_id.get(pair2, None):
                    continue
                    
                if self.primary_id != "drugbank":
                    # if drug2's primary id is not drugbank, get mapping for it
                    target = self.add_prefix_to_id("something", self.drugbank_to_primary_id[pair2])  
                else:
                    target = self.add_prefix_to_id("drugbank", pair2)
                
                self.ddi_edge_list.append((None, source, target, label, {}))
                
                counter += 1
                
            if self.early_stopping and counter >= self.early_stopping:
                break
                
                
    def get_drug_nodes(self, label="drug"):
        """
        Get drug nodes
        Args:
            label = label of nodes
        """
        
        self.node_list = []
        
        counter = 0
        
        logger.debug('Getting Drugbank nodes')
        
        for drug in tqdm(self.drugbank_drugs_detailed):
            
            # if node's primary id is not drugbank and it doesnt have any mapping skip it
            if self.primary_id != "drugbank" and self.drugbank_to_primary_id and not self.drugbank_to_primary_id.get(drug.drugbank_id, None):
                continue

            if self.primary_id != "drugbank":
                # if primary id is not drugbank, get mapping for it
                node_id = self.add_prefix_to_id("something", self.drugbank_to_primary_id[drug.drugbank_id])
            else:
                node_id = self.add_prefix_to_id("drugbank", drug.drugbank_id)
            
            # create props dict
            props = {}
            
            # get node properties
            if self.primary_node_fields:
                props = props | {str(k).replace(" ","_").lower():v for k, v in drug._asdict().items() if k in self.primary_node_fields}

            if self.property_node_fields and self.drugbank_properties.get(drug.drugbank_id, None):
                props = props | {str(k).replace(" ","_").lower():v for k, v in self.drugbank_properties[drug.drugbank_id].items() if k in self.property_node_fields}

            if self.drugbank_external_node_fields and self.drugbank_drugs_external_ids.get(drug.drugbank_id, None):
                props = props | {str(k).replace(" ","_").lower():v for k, v in self.drugbank_drugs_external_ids[drug.drugbank_id].items() if k in self.drugbank_external_node_fields}

            if self.unichem_mapping_node_fields:
                for db, mapping_dict in self.unichem_mappings.items():
                    if mapping_dict.get(drug.drugbank_id, None):
                        props[str(db).replace(" ","_").lower()] = mapping_dict[drug.drugbank_id]


            self.node_list.append((node_id, label, props))
            
            counter += 1
            
            if self.early_stopping and counter >= self.early_stopping:
                break
