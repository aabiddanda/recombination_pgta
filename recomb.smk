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
                f"results/aneuploidy_calls/{f}/{e}.aneuploidy.tsv"
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


rule est_ploidy_and_params:
    """Estimate local ploidy and per-chromosome parameters"""
    input:
        baf_pkl=rules.preprocess_data.output.baf_pkl,
    output:
        hmm_pkl="results/aneuploidy_calls/{family}/{embryo}.hmm.pkl.gz",
    params:
        mother_id=lambda wildcards: config["families"][wildcards.family]["maternal"][
            "id"
        ],
        father_id=lambda wildcards: config["families"][wildcards.family]["paternal"][
            "id"
        ],
        child_id=lambda wildcards: wildcards.embryo,
        unphased=False,
    script:
        "scripts/baf_hmm.py"


def bayes_factor(posteriors, priors=None):
    """Compute Bayes Factors for evidence of specific aneuploidy states."""
    if priors is None:
        priors = np.ones(posteriors.size) / posteriors.size

    assert posteriors.size == priors.size
    assert np.isclose(np.sum(priors), 1.0)
    bfs = np.zeros(posteriors.size)
    for i in range(posteriors.size):
        denom = np.sum(
            [posteriors[j] * priors[i] for j in range(posteriors.size) if j != i]
        )
        bfs[i] = posteriors[i] * (1 - priors[i]) / denom
    return bfs


rule reformat_aneuploidy_calls:
    input:
        hmm_pkl=rules.est_ploidy_and_params.output.hmm_pkl,
    output:
        aneuploidy_tsv="results/aneuploidy_calls/{family}/{embryo}.aneuploidy.tsv",
    run:
        with open(output.aneuploidy_tsv, "w") as out:
            full_hmm_output = pickle.load(gzip.open(input.hmm_pkl, "r"))
            chroms = full_hmm_output.keys()
            out.write(
                "mother\tfather\tchild\tchrom\tsigma_baf\tpi0_baf\t0\t1m\t1p\t2\t3m\t3p\tbf_max\tbf_max_cat\tprop01_exp_het\n"
            )
            cats = np.array(["0", "1m", "1p", "2", "3m", "3p"])
            for c in chroms:
                data = full_hmm_output[c]
                post_vals = np.array([data[x] for x in cats])
                bayes_factor_chrom = bayes_factor(post_vals)
                max_bf = np.max(bayes_factor_chrom)
                max_cat = cats[np.argmax(bayes_factor_chrom)]
                out.write(
                    f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{c}\t{data['sigma_baf']}\t{data['pi0_baf']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\t{max_bf}\t{max_cat}\t{data['prop01_exp_het']}\n"
                )


rule mendelian_filtering:
    """Filtering to mendelian violating variants in disomic embryos."""
    input:
        baf_pkl=lambda wildcards: expand(
            "results/preprocessed_data/{{family}}/{embryo}.bafs.pkl.gz",
            embryo=config["families"][wildcards.family]["embryos"].keys(),
        ),
        aneuploidy_calls=lambda wildcards: expand(
            "results/aneuploidy_calls/{{family}}/{embryo}.aneuploidy.tsv",
            embryo=config["families"][wildcards.family]["embryos"].keys(),
        ),
    output:
        mendelian_errors="results/mendelian_errors/{family}.mendelian_errors.tsv",
    params:
        ppThresh=0.90,
        phaseCorrect=False,
        eps=1e-2,
        err_rate=5e-2,
    resources:
        time="00:05:00",
    script:
        "scripts/mendelian_errors.py"


rule est_crossover_euploid_chrom_trio_heuristic:
    """Estimate crossover events on disomic chromosomes in trio datasets using the heuristic approach of Coop et al 2008."""
    input:
        baf_pkl=lambda wildcards: expand(
            "results/preprocessed_data/{{family}}/{embryo}.bafs.pkl.gz",
            embryo=config["families"][wildcards.family]["embryos"].keys(),
        ),
        hmm_pkl=lambda wildcards: expand(
            "results/aneuploidy_calls/{{family}}/{embryo}.hmm.pkl.gz",
            embryo=config["families"][wildcards.family]["embryos"].keys(),
        ),
        mendelian_errors=rules.mendelian_filtering.output.mendelian_errors,
        aneuploidy_calls=lambda wildcards: expand(
            "results/aneuploidy_calls/{{family}}/{embryo}.aneuploidy.tsv",
            embryo=config["families"][wildcards.family]["embryos"].keys(),
        ),
    output:
        est_recomb="results/natera_inference_heuristic/{family}.est_recomb.tsv",
    params:
        use_prev_params=True,
        ppThresh=0.90,
        phaseCorrect=True,
    resources:
        time="1:00:00",
        mem_mb="10G",
    script:
        "scripts/sibling_rec_est.py"
