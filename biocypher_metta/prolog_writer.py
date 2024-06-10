# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher import BioCypher
import pathlib
import os
from biocypher._logger import logger
import networkx as nx

class PrologWriter:

    def __init__(self, schema_config, biocypher_config,
                 output_dir):
        self.schema_config = schema_config
        self.biocypher_config = biocypher_config
        self.output_path = pathlib.Path(output_dir)

        if not os.path.exists(output_dir):
            self.output_path.mkdir()

        self.bcy = BioCypher(schema_config_path=schema_config,
                             biocypher_config_path=biocypher_config)

        self.onotology = self.bcy._get_ontology()
        self.create_edge_types()
        #self.excluded_properties = ["licence", "version", "source"]
        self.excluded_properties = []


    def create_edge_types(self):
        schema = self.bcy._get_ontology_mapping()._extend_schema()
        self.edge_node_types = {}

        for k, v in schema.items():
            if v["represented_as"] == "edge": #(: (label $x $y) (-> source_type target_type
                edge_type = self.convert_input_labels(k)

                # ## TODO fix this in the scheme config
                if isinstance(v["input_label"], list):
                    label = self.convert_input_labels(v["input_label"][0])
                    source_type = self.convert_input_labels(v["source"][0])
                    target_type = self.convert_input_labels(v["target"][0])
                else:
                    label = self.convert_input_labels(v["input_label"])
                    source_type = self.convert_input_labels(v["source"])
                    target_type = self.convert_input_labels(v["target"])
                self.edge_node_types[label.lower()] = {"source": source_type.lower(), "target": target_type.lower()}

    def write_nodes(self, nodes, path_prefix=None, create_dir=True):
        if path_prefix is not None:
            file_path = f"{self.output_path}/{path_prefix}/nodes.pl"
            if create_dir:
                if not os.path.exists(f"{self.output_path}/{path_prefix}"):
                    pathlib.Path(f"{self.output_path}/{path_prefix}").mkdir(parents=True, exist_ok=True)
        else:
            file_path = f"{self.output_path}/nodes.pl"
        with open(file_path, "a") as f:
            for node in nodes:
                out_str = self.write_node(node)
                for s in out_str:
                    f.write(s + "\n")

            f.write("\n")

        logger.info("Finished writing out nodes")

    def write_edges(self, edges, path_prefix=None, create_dir=True):
        if path_prefix is not None:
            file_path = f"{self.output_path}/{path_prefix}/edges.pl"
            if create_dir:
                if not os.path.exists(f"{self.output_path}/{path_prefix}"):
                    pathlib.Path(f"{self.output_path}/{path_prefix}").mkdir(parents=True, exist_ok=True)
        else:
            file_path = f"{self.output_path}/edges.pl"

        with open(file_path, "a") as f:
            for edge in edges:
                out_str = self.write_edge(edge)
                for s in out_str:
                    f.write(s + "\n")

            f.write("\n")

    def write_node(self, node):
        id, label, properties = node
        if "." in label:
            label = label.split(".")[1]
        label = label.lower()
        id = id.lower()
        def_out = f"{self.convert_input_labels(label)}({id})"
        return self.write_property(def_out, properties)

    def write_edge(self, edge):
        source_id, target_id, label, properties = edge
        label = label.lower()
        source_id = source_id.lower()
        target_id = target_id.lower()
        source_type = self.edge_node_types[label]["source"]
        target_type = self.edge_node_types[label]["target"]
        def_out = f"{label}({source_type}({source_id}), {target_type}({target_id}))"
        return self.write_property(def_out, properties)


    def write_property(self, def_out, property):
        out_str = [f"{def_out}."]
        for k, v in property.items():
            if k in self.excluded_properties or v is None or v == "": continue
            if isinstance(v, list):
                prop = "["
                for i, e in enumerate(v):
                    prop += f'{self.check_property(e)}'
                    if i != len(v) - 1: prop += ","
                prop += "]."
            elif isinstance(v, dict):
                prop = f"{k}({def_out})."
                out_str.extend(self.write_property(prop, v))
            else:
                out_str.append(f'{k}({def_out}, {self.check_property(v)}).')
        return out_str

    def check_property(self, prop):
        if isinstance(prop, str):
            if " " in prop:
                prop = prop.replace(" ", "_")

            special_chars = ["(", ")"]
            escape_char = "\\"
            return "".join(escape_char + c if c in special_chars or c == escape_char else c for c in prop)

        return prop

    def convert_input_labels(self, label, replace_char="_"):
        """
        A method that removes spaces in input labels and replaces them with replace_char
        :param label: Input label of a node or edge
        :param replace_char: the character to replace spaces with
        :return:
        """
        return label.replace(" ", replace_char)

    def get_parent(self, G, node):
        """
        Get the immediate parent of a node in the ontology.
        """
        return nx.dfs_preorder_nodes(G, node, depth_limit=2)

    def show_ontology_structure(self):
        self.bcy.show_ontology_structure()

    def summary(self):
        self.bcy.summary()