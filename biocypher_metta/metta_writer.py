# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher import BioCypher
import pathlib
import os
import logging as log
import networkx as nx

class MeTTaWriter:

    def __init__(self, schema_config, biocypher_config,
                 output_dir):
        self.schema_config = schema_config
        self.biocypher_config = biocypher_config
        self.output_path = pathlib.Path(output_dir)

        if not os.path.exists(output_dir):
            log.info(f"Directory {output_dir} doesn't exist. Creating it...")
            self.output_path.mkdir()

        self.bcy = BioCypher(schema_config_path=schema_config,
                             biocypher_config_path=biocypher_config)

        self.onotology = self.bcy._get_ontology()
        self.create_type_hierarchy()

    def create_type_hierarchy(self):
        G = self.onotology._nx_graph
        file_path = f"{self.output_path}/type_defs.metta"
        with open(file_path, "w") as f:
            for node in G.nodes:
                if "mixin" in node: continue
                ancestor = list(self.get_parent(G, node))[-1]
                node = self.convert_input_labels(node)
                ancestor = self.convert_input_labels(ancestor)
                if ancestor == node:
                    f.write(f"(: {node.upper()} Type)\n")

                else:
                    f.write(f"(: {node.upper()} {ancestor.upper()})\n")

            self.create_data_constructors(f)

        log.info("Type hierarchy created successfully.")

    def create_data_constructors(self, file):
        schema = self.onotology.mapping._extend_schema()

        def edge_data_constructor(edge_type, source_type, target_type, label):
            return f"(: ({label.lower()} $x $y) (-> {source_type.upper()} {target_type.upper()} {edge_type.upper()})"

        def node_data_constructor(node_type, node_label):
            return f"(: ({node_label.lower()} $x) (-> $x {node_type.upper()}))"

        for k, v in schema.items():
            if v["represented_as"] == "edge": #(: (label $x $y) (-> source_type target_type
                edge_type = self.convert_input_labels(k)
                label = self.convert_input_labels(v["input_label"])
                source_type = self.convert_input_labels(v["source"])
                target_type = self.convert_input_labels(v["target"])
                out_str = edge_data_constructor(edge_type, source_type, target_type, label)
                file.write(out_str + "\n")

            elif v["represented_as"] == "node":
                label = v["input_label"]
                if not isinstance(label, list):
                    label = [label]

                label = [self.convert_input_labels(l) for l in label]
                node_type = self.convert_input_labels(k)
                for l in label:
                    out_str = node_data_constructor(node_type, l)
                    file.write(out_str + "\n")


    def write_nodes(self, nodes):
        file_path = f"{self.output_path}/nodes.metta"
        with open(file_path, "w") as f:
            for node in nodes:
                out_str = self.write_node(node)
                for s in out_str:
                    f.write(s + "\n")

        log.info("Finished writing out nodes")



    def write_edges(self, edges):
        file_path = f"{self.output_path}/edges.metta"

        with open(file_path, "w") as f:
            for edge in edges:
                out_str = self.write_edge(edge)
                for s in out_str:
                    f.write(s + "\n")


    def write_node(self, node):
        id, label, properties = node
        if "." in label:
            label = label.split(".")[1]
        out_str = []
        def_out = f"({self.convert_input_labels(label)} {id})"
        out_str.append(def_out)
        for k, v in properties.items():
            out_str.append(f"(has-property {k} {def_out} {v})")

        return out_str

    def write_edge(self, edge):
        # label = self.convert_input_labels(edge["label"])
        # source_id = edge["source"]
        # target_id = edge["target"]
        # properties = edge.get_properties()
        id, source_id, target_id, label, properties = edge
        rel = properties["relationship"]
        source_type, target_type = rel["from"], rel["to"]
        del properties["relationship"]
        def_out = f"({label} ({source_type} {source_id}) ({target_type} {target_id}))"
        out_str = [def_out]
        for k, v in properties.items():
            out_str.append(f"(has-property {k} {def_out} {v})")

        return out_str


    def write_property(self, node):
        pass

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