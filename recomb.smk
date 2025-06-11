#!python3

import numpy as np
import pandas as pd

import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO


configfile: "config.yaml"


# ------- Rules Section ------- #
localrules:
    all,


rule all:
    input:
        [
            [
                f"results/preprocessed_data/{f}/{e}.bafs.pkl.gz"
                for e in config["families"][f]["embryos"].keys()
            ]
            for f in config["families"].keys()
        ],


rule preprocess_data:
    """Preprocess some of the data into a compressed numpy array format."""
    input:
        maternal_data=lambda wildcards: config["families"][wildcards.family][
            "maternal"
        ]["path"],
        paternal_data=lambda wildcards: config["families"][wildcards.family][
            "paternal"
        ]["path"],
        embryo_data=lambda wildcards: config["families"][wildcards.family]["embryos"][
            wildcards.embryo
        ],
    output:
        baf_pkl="results/preprocessed_data/{family}/{embryo}.bafs.pkl.gz",
    script:
        "scripts/preprocess.py"


# rule mendelian_filtering:
    # """Filtering to non-mendelian violating variants."""
    # input:
        # lambda wildcards: expand(
            # "results/preprocessed_data/{{family}}/{embryo}.bafs.pkl.gz",
            # embryo=config["families"][wildcards.family]["embryos"].keys(),
        # ),
    # output:
        # mendelian_errors="results/mendelian_errors/{family}.mendelian_errors.tsv",
    # script:
        # "scripts/mendelian_errors.py"


# rule est_ploidy_and_params:
    # """Estimate local ploidy and per-chromosome parameters"""
    # input:
        # baf_pkl=rules.preprocess_data.output.baf_pkl,
    # output:
        # "",
    # script:
        # "scripts/baf_hmm.py"


# rule est_crossover_euploid_chrom_trio_heuristic:
    # """Estimate crossover events for euploid chromosomes in trio datasets using the heuristic approach of Coop et al 2007."""
    # input:
        # triplets="results/natera_inference/valid_trios.triplets.txt",
        # baf_pkl=lambda wildcards: define_triplets_baf(
            # mother_id=wildcards.mother, father_id=wildcards.father
        # ),
        # hmm_pkl=lambda wildcards: define_triplets_hmm(
            # mother_id=wildcards.mother, father_id=wildcards.father
        # ),
    # output:
        # est_recomb="results/natera_inference_heuristic/{mother}+{father}.est_recomb.tsv",
    # params:
        # chroms=chroms,
        # use_prev_params=True,
        # ppThresh=0.90,
        # phaseCorrect=True,
    # resources:
        # time="1:00:00",
        # mem_mb="10G",
    # script:
        # "scripts/sibling_rec_est.py"
