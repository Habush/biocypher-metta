# BioCypher MeTTa
A project for creating [BioCypher-driven](https://github.com/biocypher/biocypher) knowledge graphs and outputs 
[MeTTa](https://wiki.opencog.org/w/File:MeTTa_Specification.pdf) files.


## ⚙️ Installation (local)
1. Clone this repository.
```{bash}
git clone https://github.com/Habush/biocypher-metta.git
```

2. Install the dependencies using [Poetry](https://python-poetry.org/). (Or feel
free to use your own dependency management system. We provide a `pyproject.toml`
to define dependencies.)
```{bash}
poetry install
```
3. You are ready to go!
```{bash}
poetry shell
python create_knowledge_graph.py
```


## 🛠 Usage

### Structure
The project template is structured as follows:
```
.
│  # Project setup
│
├── LICENSE
├── README.md
├── pyproject.toml
│
│  # Docker setup
│
├── Dockerfile
├── docker
│   ├── biocypher_entrypoint_patch.sh
│   ├── create_table.sh
│   └── import.sh
├── docker-compose.yml
├── docker-variables.env
│
│  # Project pipeline
|── biocypher_metta
│   ├── adapters
│   ├── metta_writer.py
│
├── create_knowledge_graph.py
├── config
│   ├── biocypher_config.yaml
│   ├── biocypher_docker_config.yaml
│   └── schema_config.yaml
```

The main components of the BioCypher pipeline are the
`create_knowledge_graph.py`, the configuration in the `config` directory, and
the adapter module in the `biocypher_metta` directory. The input adapters are used for preprocessing biomedical 
databases and converting them into a BioCypher nodes and edges. The `metta_writer.py` script contains code to convert 
these BioCypher nodes and edges into MeTTa represntation.
