import os
import sys
import pandas as pd
import numpy as np
import collections
from typing import Union

from pathlib import Path
from time import time

from pypath.inputs import intact, uniprot
from pypath.share import curl, settings

from tqdm import tqdm # progress bar

from biocypher._logger import logger
from pypath.resources import urls
from contextlib import ExitStack

from bioregistry import normalize_curie

from enum import Enum

global adapter_name
adapter_name = "intact" # might be useful

class IntactEdgeFields(Enum):
    SOURCE = "source"
    PUBMED_IDS = "pubmeds"
    INTACT_SCORE = "mi_score"
    METHODS = "methods"
    INTERACTION_TYPES = "interaction_types"
    

class IntAct:
    def __init__(self, output_dir = None, export_csvs = False, split_output = False, cache=False, debug=False, retries=6,
                organism=9606, intact_fields: Union[None, list[IntactEdgeFields]] = None, add_prefix = True, test_mode = False, aggregate_pubmed_ids: bool = True,
                aggregate_methods: bool = True):
        """
        Downloads and processes IntAct data

            Args:
                export_csvs: Flag for whether or not create csvs of outputs of databases
                split_csvs: whether or not to split output csv files to multiple parts
                cache: if True, it uses the cached version of the data, otherwise
                forces download.
                debug: if True, turns on debug mode in pypath.
                retries: number of retries in case of download error.
                organism: taxonomy id number of selected organism, if it is None, downloads all organism data.
                intact_fields: intact fields to be used in the graph.
                add_prefix: if True, add prefix to uniprot ids
                test_mode: if True, take small sample from data for testing
                aggregate_pubmed_ids: if True, collects all pubmed ids that belongs to same protein pair
                aggregate_methods = if True, collects all experiemental methods that belongs to same protein pair
        """
        
        self.export_csvs = export_csvs
        self.split_output = split_output
        self.swissprots = list(uniprot._all_uniprots("*", True))
        self.cache = cache
        self.debug = debug
        self.retries = retries        
        self.organism = organism
        self.intact_fields = intact_fields
        self.add_prefix = add_prefix
        self.test_mode = test_mode
        
        self.aggregate_dict = {IntactEdgeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              IntactEdgeFields.METHODS.value:aggregate_methods}

        
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

    def download_intact_data(self):
        """
        Wrapper function to download IntAct data using pypath; used to access
        settings.
            
        To do: Make arguments of intact.intact_interactions selectable for user.
        """
                     
        logger.debug("Started downloading IntAct data")
        logger.info(f"This is the link of IntAct data we downloaded:{urls.urls['intact']['mitab']}. Please check if it is up to date")
        t0 = time()
        
        with ExitStack() as stack:                         
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())

            self.intact_ints = intact.intact_interactions(miscore=0, organism=self.organism, complex_expansion=True, only_proteins=True)
        
        if self.test_mode:
            self.intact_ints = self.intact_ints[:100]
        
        t1 = time()
        
        logger.info(f'IntAct data is downloaded in {round((t1-t0) / 60, 2)} mins')


    def intact_process(self, rename_selected_fields: Union[None, list[str]] = None):
        """
        Processor function for IntAct data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. Also, it filters
        protein pairs found in swissprot.
        
        Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        # 
        selected_fields = self.set_edge_fields()
        
        # if rename_selected_fields is not defined create column names from this dictionary
        default_field_names = {"source":"source", "pubmeds":"pubmed_ids", "mi_score":"intact_score", "methods":"method", "interaction_types":"interaction_type"}
        
        self.intact_field_new_names = {}
        
        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception("Length of selected_fields variable should be equal to length of rename_selected_fields variable")
            
                
            for field_old_name, field_new_name in list(zip(selected_fields, rename_selected_fields)):                    
                self.intact_field_new_names[field_old_name] = field_new_name
            
            self.intact_field_new_names["id_a"] = "uniprot_a"
            self.intact_field_new_names["id_b"] = "uniprot_b"
        
        else:            
            for field_old_name in selected_fields:
                self.intact_field_new_names[field_old_name] = default_field_names[field_old_name]
                
            self.intact_field_new_names["id_a"] = "uniprot_a"
            self.intact_field_new_names["id_b"] = "uniprot_b"        
       
            
        logger.debug("Started processing IntAct data")
        t1 = time()
                         
        # create dataframe            
        intact_df = pd.DataFrame.from_records(self.intact_ints, columns=self.intact_ints[0]._fields)
        
        # turn list columns to string
        for list_column in ["pubmeds", "methods", "interaction_types"]:
            intact_df[list_column] = [';'.join(map(str, l)) for l in intact_df[list_column]]
        
        intact_df.fillna(value=np.nan, inplace=True)
        
        # add source database info
        intact_df["source"] = "IntAct"
        
        # filter selected fields
        intact_df = intact_df[list(self.intact_field_new_names.keys())]
        
        # rename columns
        intact_df.rename(columns=self.intact_field_new_names, inplace=True)
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        intact_df = intact_df[(intact_df["uniprot_a"].isin(self.swissprots)) & (intact_df["uniprot_b"].isin(self.swissprots))]
        intact_df.reset_index(drop=True, inplace=True)
        
        if "pubmeds" in self.intact_field_new_names.keys():            
            # assing pubmed ids that contain unassigned to NaN value 
            intact_df[self.intact_field_new_names["pubmeds"]].loc[intact_df[self.intact_field_new_names["pubmeds"]].astype(str).str.contains("unassigned", na=False)] = np.nan
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the pair with the highest score and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same interaction type with b x a pair, drop b x a pair
        if "mi_score" in self.intact_field_new_names.keys():            
            intact_df.sort_values(by=self.intact_field_new_names["mi_score"], ascending=False, inplace=True)
        
        intact_df_unique = intact_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)
        
        
        def aggregate_fields(element):
            element = "|".join([str(e) for e in set(element.dropna())])
            if not element:
                return np.nan
            else:
                return element
        
        if any(list(self.aggregate_dict.values())):
            agg_field_list = [k for k, v in self.aggregate_dict.items() if v]
            
            agg_dict = {}            
            for k, v in self.intact_field_new_names.items():
                if k in agg_field_list:
                    agg_dict[v] = aggregate_fields
                else:                
                    agg_dict[v] = "first"
            
            intact_df_unique = intact_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate(agg_dict)

            
        #intact_df_unique["pubmed_id"].replace("", np.nan, inplace=True) # replace empty string with NaN
        
        if "interaction_types" in self.intact_field_new_names.keys():
            intact_df_unique = intact_df_unique[~intact_df_unique[["uniprot_a", "uniprot_b", self.intact_field_new_names["interaction_types"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            intact_df_unique = intact_df_unique[~intact_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        
        if self.export_csvs:
            intact_output_path = self.export_dataframe(intact_df_unique, "intact")
            logger.info(f'Final IntAct data is written: {intact_output_path}')

        self.final_intact_ints = intact_df_unique
        
        t2 = time()
        logger.info(f'IntAct data is processed in {round((t2-t1) / 60, 2)} mins')                        
            
    def set_edge_fields(self) -> list:
        """
        Sets intact edge fields
        Returns:
            selected field list
        """
        if self.intact_fields is None:
            return [field.value for field in IntactEdgeFields]
        else:
            return [field.value for field in self.intact_fields]
        
        
    def add_prefix_to_id(self, prefix="uniprot", identifier=None, sep=":") -> str:
        """
        Adds prefix to uniprot id
        """
        if self.add_prefix and identifier:
            return normalize_curie( prefix + sep + str(identifier))
        
        return identifier         
            
        
    def get_intact_edges(self) -> list:
        """
        Get PPI edges from intact data 
        """
        
        # create edge list
        edge_list = []
        
        for index, row in tqdm(self.final_intact_ints.iterrows(), total=self.final_intact_ints.shape[0]):
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
