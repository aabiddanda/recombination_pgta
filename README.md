# Recombination Calling from PGT-A Data

This pipeline implements the procedures used for calling recombination events in PGT-A datasets. This is a subset of routines used in: https://github.com/mccoy-lab/natera_recomb


## Installation

In order to install the pipeline, you can use a viable `conda` distribution and call the following steps:

```
git clone 
cd recombination_pgta/
conda env create -f env.yaml
conda activate recombination_pgta 
```

## Pipeline 

Briefly - the pipeline proceeds according to the following steps per-chromosome:

1. Quality control of mendelian violating variants across all embryos within a family
2. Phasing chromosomes using minimal haplotype switch criteria between embryos 
3. Parameter estimation for noise parameter within a specific embryo
4. Identification of crossovers at informative markers for each embryo
5. Refinement of crossover positioning using `QuadHMM` for differences in haplotype sharing


### Configuration

The configuration of the pipeline is a critical step as it allows the user to specify family structure for analysis, which have example blocks like: 

```
families:
    family1:
      maternal:
        id: "B"
        path: "data/family1/maternal.autosomes.tsv.gz"
      paternal:
        id: "A"
        path: "data/family1/paternal.autosomes.tsv.gz"
      embryos:
        embryo1: "data/family1/embryos/embryo0.autosomes.tsv.gz"
        embryo2: "data/family1/embryos/embryo1.autosomes.tsv.gz"
        embryo3: "data/family1/embryos/embryo2.autosomes.tsv.gz"
        embryo4: "data/family1/embryos/embryo3.autosomes.tsv.gz"
```

NOTE: the names "embryo1" can be interchanged with other names (e.g. array ID) in case that might help with downstream intersection with metadata. 

You can see examples of a data directory here: https://www.dropbox.com/scl/fi/ncf7b2i9kx56emoeyijaf/data.tar.gz?rlkey=iv70iw56kmdpktk7v19dn1ehr&dl=0

The parental files contain parental haplotypes (though it is not entirely necessary to phase the parents), and the embryo data contain chromosome, position, and B-allele frequency as tab-separated values (TSV) files.

### Running 

Following environment setup and configuration, you can run the pipeline using the following command:

```
snakemake -s recomb.smk -j4 -p -n
```

The `-j 1` flag specifies that there will be 4 parallel jobs running at a time. To really run the pipeline, remove the `-n` flag at the end. 

In the example data and configuration setup - the full pipeline for two families takes approximately 15 minutes on a reasonable laptop computer. 

If the pipeline finishes, recombination calls per-family will be set as ``


## Citation 

If you use any of the code here, please cite the following manuscript: 

```
Common variation in meiosis genes shapes human recombination phenotypes and aneuploidy risk
Sara A. Carioscia, Arjun Biddanda, Margaret R. Starostik, Xiaona Tang, Eva R. Hoffmann, Zachary P. Demko, Rajiv C. McCoy
medRxiv 2025.04.02.25325097; doi: https://doi.org/10.1101/2025.04.02.25325097
```
## Contact

For any questions - please either submit an issue or contact @aabiddanda
