"""
Knowledge graph generation through BioCypher script
"""
from biocypher_metta.metta_writer import *
from biocypher._logger import logger
import typer
import yaml
import importlib #for reflection
from typing_extensions import Annotated

app = typer.Typer()

# Run build
@app.command()
def main(output_dir: Annotated[pathlib.Path, typer.Option(exists=False, file_okay=False, dir_okay=True)],
         adapter_config: Annotated[pathlib.Path, typer.Option(exists=True, file_okay=True, dir_okay=False)]):
    """
    Main function. Call individual adapters to download and process data. Build
    via BioCypher from node and edge data.
    """

    # Start biocypher

    bc = MeTTaWriter(schema_config="config/schema_config.yaml",
                     biocypher_config="config/biocypher_config.yaml",
                     output_dir=output_dir)

    # bc.show_ontology_structure()

    # Run adapters

    with open(adapter_config, "r") as fp:
        try:
            adapters_config = yaml.safe_load(fp)
        except yaml.YAMLError as e:
            logger.error(f"Error while trying to load adapter config")
            logger.error(e)

    for c in adapters_config:
        logger.info(f"Running adapter: {c}")
        adapter_config = adapters_config[c]["adapter"]
        adapter_module = importlib.import_module(adapter_config["module"])
        adapter_cls = getattr(adapter_module, adapter_config["cls"])
        ctr_args = adapter_config["args"]
        adapter = adapter_cls(**ctr_args)
        write_nodes = adapters_config[c]["nodes"]
        write_edges = adapters_config[c]["edges"]
        outdir = adapters_config[c]["outdir"]

        if write_nodes:
            nodes = adapter.get_nodes()
            bc.write_nodes(nodes, path_prefix=outdir)

        if write_edges:
            edges = adapter.get_edges()
            bc.write_edges(edges, path_prefix=outdir)


    logger.info("Done")

if __name__ == "__main__":
    app()
