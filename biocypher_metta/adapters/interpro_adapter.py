from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import interpro

from contextlib import ExitStack
from bioregistry import normalize_curie
from tqdm import tqdm

from time import time
from biocypher._logger import logger

from typing import Union

from enum import Enum

logger.debug(f"Loading module {__name__}.")

class InterProNodeField(Enum):
    """
    Domain node fields in InterPro 
    """
    
    # primary attributes
    PROTEIN_COUNT = "protein_count"
    NAME = "name"
    TYPE = "type"
    PARENT_LIST = "parent_list"
    CHILD_LIST = "child_list"
    
    # member list attributes
    PFAM = "PFAM"
    
    # external attributes
    EC = "EC"
    
    # structural attributes
    PDB = "PDB"
    
    @classmethod
    def get_primary_attributes(cls):
        """
        Returns primary InterPro attributes
        """
        return [cls.PROTEIN_COUNT.value, cls.NAME.value, cls.TYPE.value,
               cls.PARENT_LIST.value, cls.CHILD_LIST.value]
    
    @classmethod
    def get_member_list_attributes(cls):
        """
        Returns external InterPro attributes
        """
        return [cls.PFAM.value]
    
    @classmethod
    def get_external_attributes(cls):
        """
        Returns external InterPro attributes
        """
        return [cls.EC.value]
    
    @classmethod
    def get_structural_attributes(cls):
        """
        Returns structural InterPro attributes
        """
        return [cls.PDB.value]
    
    
class InterProEdgeField(Enum):
    """
    Domain edge fields in InterPro 
    """
    START = "start"
    END = "end"
    

class InterPro:
    """
    Class that downloads InterPro data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, cache=False, debug=False, page_size=150, retries=6, organism=None, add_prefix = True,
                node_fields:Union[InterProNodeField, None] = None, edge_fields:Union[InterProEdgeField, None] = None,
                test_mode: bool = False):
        """        
        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
            page_size: page size of downloaded annotation data
            organism: taxonomy id of selected organism, if None take all available organisms
            add_prefix: if True, add prefix to database identifiers
            node_fields: fields that will be included in domain nodes
            edge_fields: fields that will be included in protein-domain edges
            test_mode: limits amount of data for testing
        """
        self.cache = cache
        self.debug = debug
        self.retries = retries
        self.page_size = page_size
        self.organism = organism
        self.add_prefix = add_prefix
        
        
        # set node and edge fields
        self.set_node_and_edge_fields(node_fields=node_fields, edge_fields=edge_fields)
        
        self.early_stopping = None
        if test_mode:
            self.early_stopping = 100        


    def download_domain_node_data(self) -> None:
        """
        Downloads domain node data from Interpro
        """
        
        # stack pypath context managers
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=self.retries))

            if self.debug:
                stack.enter_context(curl.debug_on())

            if not self.cache:
                stack.enter_context(curl.cache_off())
            
            logger.debug("Started downloading InterPro domain data")
            t0 = time()
            
            # To do: Filter according to tax id??
            self.interpro_entries = interpro.interpro_entries() # returns a list of namedtuples
            self.interpro_structural_xrefs = interpro.interpro_xrefs(db_type = 'structural')
            self.interpro_external_xrefs = interpro.interpro_xrefs(db_type = 'external')
            
            t1 = time()
            logger.info(f'InterPro domain data is downloaded in {round((t1-t0) / 60, 2)} mins')
            
    def download_domain_edge_data(self) -> None:
        """
        Downloads Uniprot annotation data from Interpro
        """
        
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=self.retries))

            if self.debug:
                stack.enter_context(curl.debug_on())

            if not self.cache:
                stack.enter_context(curl.cache_off())
            
            logger.debug("Started downloading InterPro annotation data")
            t0 = time()
            
            if self.organism:                
                # WARNING: decrease page_size parameter if there is a curl error about timeout in the pypath_log
                self.interpro_annotations = interpro.interpro_annotations(page_size = self.page_size, reviewed = True, tax_id = self.organism)
            else:
                self.interpro_annotations = interpro.interpro_annotations(page_size = self.page_size, reviewed = True, tax_id = '')
            
            t1 = time()
            logger.info(f'InterPro annotation data is downloaded in {round((t1-t0) / 60, 2)} mins')

    
    def get_interpro_nodes(self, node_label="domain") -> list:
        """
        Prepares InterPro domain nodes for BioCypher
        Args:
            node_label : label of interpro nodes
        """

        # create list of nodes
        node_list = []

        # define primary and external attributes
        primary_attributes = InterProNodeField.get_primary_attributes()
        member_list_attributes = InterProNodeField.get_member_list_attributes()
        external_attributes = InterProNodeField.get_external_attributes()
        structural_attributes = InterProNodeField.get_structural_attributes()
        
        logger.debug("Creating domain nodes")
        t0 = time()
        
        # set counter for early stopping
        counter = 0
        
        for entry in tqdm(self.interpro_entries):
            props = dict()
            interpro_props = entry._asdict()
            
            domain_id = self.add_prefix_to_id("interpro", entry.interpro_id)
            
            # get primary InterPro attributes
            for element in primary_attributes:
                
                if element in self.node_fields and interpro_props.get(element):
                    if element == "protein_count":
                        props[element.replace(" ", "_").lower()] = int(interpro_props.get(element))
                    else:
                        props[element.replace(" ", "_").lower()] = self.check_length(interpro_props.get(element))
            
            # get member list InterPro attributes
            for element in member_list_attributes:
                if element in self.node_fields and interpro_props.get('member_list').get(element):
                    props[element.replace(" ", "_").lower()] = self.check_length(interpro_props.get('member_list').get(element))
                        
            # get external InterPro attributes            
            for element in external_attributes:
                if element in self.node_fields and self.interpro_external_xrefs.get(entry.interpro_id).get(element):
                    props[element.replace(" ", "_").lower()] = self.check_length(self.interpro_external_xrefs.get(entry.interpro_id).get(element))
            
            # get structural InterPro attributes
            for element in structural_attributes:
                    if element in self.node_fields and self.interpro_structural_xrefs.get(entry.interpro_id).get(element):
                        props[element.replace(" ", "_").lower()] = self.check_length(self.interpro_structural_xrefs.get(entry.interpro_id).get(element))
              
            node_list.append((domain_id, node_label, props))
            
            counter += 1
            
            if self.early_stopping and counter == self.early_stopping:
                break
        
        t1 = time()
        logger.info(f'InterPro nodes created in {round((t1-t0) / 60, 2)} mins')

        return node_list
    
    def get_interpro_edges(self, edge_label="protein_has_domain") -> list:
        """
        Prepares Protein-Domain edges for BioCypher
        Args:
            edge_label: label of protein-domain edge
        """
        
        # create list of edges
        edge_list = []
        
        logger.debug("Creating protein-domain edges")
        t0 = time()
        
        # set counter for early stopping
        counter = 0

        # DOMAIN-PROTEIN EDGES
        for k, v in tqdm(self.interpro_annotations.items()):            
            # k -> uniprot id
            for annotation in v:
                
                interpro_props = annotation._asdict()
                props = dict()
                
                for field in self.edge_fields:
                    if interpro_props.get(field, None):
                        props[field.replace(" ","_").lower()] = self.check_length(interpro_props[field])
                    
                interpro_id = self.add_prefix_to_id("interpro", annotation.interpro_id)
                uniprot_id = self.add_prefix_to_id("uniprot", k)
                
                edge_list.append((None, uniprot_id, interpro_id, edge_label, props))
                
                counter += 1
                
            if self.early_stopping and counter >= self.early_stopping:
                break
                
        t1 = time()
        logger.info(f'InterPro edges created in {round((t1-t0) / 60, 2)} mins')
        
        return edge_list
    
    def check_length(self, element:str) -> str | list:
        """
        If the type of given entry is a list and has just one element returns this one element
        """
        if isinstance(element, list) and len(element) == 1:
            return element[0]
        else:
            return element
        
    def add_prefix_to_id(self, prefix, identifier, sep=":") -> str:
        """
        Adds prefix to ids
        """
        if self.add_prefix:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier      
        
    def set_node_and_edge_fields(self, node_fields, edge_fields) -> None:
        """
        Sets Interpro node and edge fields 
        """
        
        if node_fields:
            self.node_fields = [field.value for field in node_fields]            
        else:
            self.node_fields = [field.value for field in InterProNodeField]
            
        if edge_fields:
            self.edge_fields = [field.value for field in edge_fields]
        else:
            self.edge_fields = [field.value for field in InterProEdgeField]
        
