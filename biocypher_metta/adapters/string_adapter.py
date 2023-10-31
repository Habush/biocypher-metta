import os
import sys
import pandas as pd
import numpy as np
import time
import collections
from typing import Union

from pathlib import Path
from time import time

from pypath.inputs import string
from pypath.inputs import biogrid, uniprot
from pypath.share import curl, settings

from tqdm import tqdm # progress bar

from biocypher._logger import logger
from pypath.resources import urls
from contextlib import ExitStack

from bioregistry import normalize_curie

from enum import Enum

global adapter_name
adapter_name = "string" # might be useful for future
    
class StringEdgeFields(Enum):
    SOURCE = "source"
    COMBINED_SCORE = "combined_score"
    PHYSICAL_COMBINED_SCORE = "physical_combined_score"
    

class STRING:
    def __init__(self, output_dir = None, export_csvs = False, split_output = False, cache=False, debug=False, retries=6,
                organism=9606, string_fields: Union[None, list[StringEdgeFields]] = None, add_prefix = True, test_mode = False):
        """
        Downloads and processes STRING data

            Args:
                export_csvs: Flag for whether or not create csvs of outputs of databases
                split_csvs: whether or not to split output csv files to multiple parts
                cache: if True, it uses the cached version of the data, otherwise
                forces download.
                debug: if True, turns on debug mode in pypath.
                retries: number of retries in case of download error.
                organism: taxonomy id number of selected organism, if it is None, downloads all organism data.
                string_fields: string fields to be used in the graph.
                add_prefix: if True, add prefix to uniprot ids
                test_mode: if True, take small sample from data for testing
                
        """
        
        self.export_csvs = export_csvs
        self.split_output = split_output
        self.swissprots = list(uniprot._all_uniprots("*", True))
        self.cache = cache
        self.debug = debug
        self.retries = retries        
        self.organism = organism
        self.string_fields = string_fields
        self.add_prefix = add_prefix
        self.test_mode = test_mode

        
        if export_csvs:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            self.output_dir = output_dir

    
    def export_dataframe(self, dataframe, data_label):
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # TODO: activate this block after adding n_rows_in_file setting to config file
        # if self.split_output:
        #     n_chunks = round(len(dataframe) / n_rows_in_file)
            # output_path =  os.path.join(self.output_dir, data_label)
            # for id, chunk in  enumerate(np.array_split(dataframe, n_chunks)):
            #     chunk.to_csv(os.path.join(output_path, f"crossbar_ppi_data_{data_label}_{id+1}.csv"), index=False)
        # else:
        output_path = os.path.join(self.output_dir, f"crossbar_ppi_data_{data_label}.csv")
        dataframe.to_csv(output_path, index=False)

        return output_path
    

    def download_string_data(self):
        """
        Wrapper function to download STRING data using pypath; used to access
        settings.
            
        To do: Make arguments of string.string_links_interactions selectable for user.
        """

        t0 = time()

        with ExitStack() as stack:
            
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())

            if self.organism is None:
                string_species = string.string_species()
                self.tax_ids = list(string_species.keys())
            else:
                self.tax_ids = [self.organism]
        
            # map string ids to swissprot ids
            uniprot_to_string = uniprot.uniprot_data("xref_string", "*", True)
            
            self.string_to_uniprot = collections.defaultdict(list)
            for k,v in uniprot_to_string.items():
                for string_id in list(filter(None, v.split(";"))):
                    self.string_to_uniprot[string_id.split(".")[1]].append(k)

            self.string_ints = []
            
            logger.debug("Started downloading STRING data")
            logger.info(f"This is the link of STRING data we downloaded:{urls.urls['string']['links']}. Please check if it is up to date")
            
            # this tax id give an error
            tax_ids_to_be_skipped = ['4565', ]
            
            # it may take around 100 hours to download whole data
            for tax in tqdm(self.tax_ids):
                if tax not in tax_ids_to_be_skipped:
                    # remove proteins that does not have swissprot ids
                    organism_string_ints = [
                        i for i in string.string_links_interactions(ncbi_tax_id=int(tax), score_threshold="high_confidence")
                        if i.protein_a in self.string_to_uniprot and i.protein_b in self.string_to_uniprot]
                    
                    logger.debug(f"Downloaded STRING data with taxonomy id {str(tax)}")
                    
                    if organism_string_ints:
                        self.string_ints.extend(organism_string_ints)
        
        if self.test_mode:
            self.string_ints = self.string_ints[:100]
            
        t1 = time()
        logger.info(f'STRING data is downloaded in {round((t1-t0) / 60, 2)} mins')
                         

    def string_process(self, rename_selected_fields: Union[None, list[str]] = None) -> None:
        """
        Processor function for STRING data. It drops duplicate and reciprocal duplicate protein pairs. In addition, it maps entries to uniprot ids 
        using crossreferences to STRING in the Uniprot data. Also, it filters protein pairs found in swissprot.
        
         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        selected_fields = self.set_edge_fields()
        
        default_field_names = {"source":"source", "combined_score":"string_combined_score", 
                               "physical_combined_score":"string_physical_combined_score"}
        
        self.string_field_new_names = {}
        
        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception("Length of selected_fields variable should be equal to length of rename_selected_fields variable")
            
            for field_old_name, field_new_name in list(zip(selected_fields, rename_selected_fields)):
                self.string_field_new_names[field_old_name] = field_new_name
            
            self.string_field_new_names["uniprot_a"] = "uniprot_a"
            self.string_field_new_names["uniprot_b"] = "uniprot_b"
        else:
            for field_old_name in selected_fields:
                self.string_field_new_names[field_old_name] = default_field_names[field_old_name]
                
            self.string_field_new_names["uniprot_a"] = "uniprot_a"
            self.string_field_new_names["uniprot_b"] = "uniprot_b"
            
        logger.debug("Started processing STRING data")
        t1 = time()
                                                  
        # create dataframe
        string_df = pd.DataFrame.from_records(self.string_ints, columns=self.string_ints[0]._fields)
        
                         
        prot_a_uniprots = []
        for protein in string_df['protein_a']:            
            id_a= (";".join(self.string_to_uniprot[protein]))
                   # if protein in self.string_to_uniprot else None) 
                   # now that we filtered interactions in line 307, we should not get KeyError here
            prot_a_uniprots.append(id_a)
                         
        prot_b_uniprots = []
        for protein in string_df['protein_b']:            
            id_b= (";".join(self.string_to_uniprot[protein]))
                   # if protein in self.string_to_uniprot else None)
                   # now that we filtered interactions in line 307, we should not get KeyError here
            prot_b_uniprots.append(id_b)
                         
        string_df["uniprot_a"] = prot_a_uniprots
        string_df["uniprot_b"] = prot_b_uniprots
        
        string_df.fillna(value=np.nan, inplace=True)
        
        # add source database info
        string_df["source"] = "STRING"
        # filter selected fields
        string_df = string_df[list(self.string_field_new_names.keys())]
        # rename columns
        string_df.rename(columns=self.string_field_new_names, inplace=True)
        
        # filter with swissprot ids
        # we already filtered interactions in line 307, we can remove this part or keep it for a double check
        # string_df = string_df[(string_df["uniprot_a"].isin(self.swissprots)) & (string_df["uniprot_b"].isin(self.swissprots))]
        # string_df.reset_index(drop=True, inplace=True)
                         
        # drop duplicates if same a x b pair exists in b x a format
        # keep the one with the highest combined score
        if "combined_score" in self.string_field_new_names.keys():            
            string_df.sort_values(by=self.string_field_new_names["combined_score"], ascending=False, inplace=True)        
            string_df_unique = string_df.dropna(subset=["uniprot_a", "uniprot_b"]).drop_duplicates(subset=["uniprot_a", "uniprot_b"], keep="first").reset_index(drop=True)
            string_df_unique = string_df_unique[~string_df_unique[["uniprot_a", "uniprot_b", self.string_field_new_names["combined_score"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            string_df_unique = string_df_unique[~string_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

            
        t2 = time()
        logger.info(f'STRING data is processed in {round((t2-t1) / 60, 2)} mins')
        
        if self.export_csvs:
            string_output_path = self.export_dataframe(string_df_unique, "string")
            logger.info(f'Final STRING data is written: {string_output_path}')

        self.final_string_ints = string_df_unique
        
    
    def set_edge_fields(self) -> list:
        """
        Sets string edge fields
        Returns:
            selected field list
        """
        if self.string_fields is None:
            return [field.value for field in StringEdgeFields]
        else:
            return [field.value for field in self.string_fields]
        
    def add_prefix_to_id(self, prefix="uniprot", identifier=None, sep=":") -> str:
        """
        Adds prefix to uniprot id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier 
    
    def get_string_edges(self) -> list:
        """
        Get PPI edges from string data
        """
        
        # create edge list
        edge_list = []
        
        for index, row in tqdm(self.final_string_ints.iterrows(), total=self.final_string_ints.shape[0]):
            _dict = row.to_dict()
            
            _source = self.add_prefix_to_id(identifier = _dict["uniprot_a"])
            _target = self.add_prefix_to_id(identifier = _dict["uniprot_b"])
            
            del _dict["uniprot_a"], _dict["uniprot_b"]
            
            _props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        # if column has multiple entries create list
                        _props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        _props[str(k).replace(" ","_").lower()] = v
           

            edge_list.append((None, _source, _target, "Interacts_With", _props))
            
        return edge_list
