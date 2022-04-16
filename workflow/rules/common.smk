from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: config['outputDir'] + "/config/snakemake_config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["sampleManifest"], sep=",").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")
