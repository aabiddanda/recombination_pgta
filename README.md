# Recombination Calling from PGT-A Data

This pipeline implements the procedures used for calling recombination events in PGT-A datasets 

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
family1:
    maternal: "mother1"
    paternal: "father1"
    embryos:
      embryo1: "data/embryo1.tsv.gz"
      embryo2: "data/embryo2.tsv.gz"
      embryo3: "data/embryo3.tsv.gz"
```

## Contact

For any questions - please either submit an issue or contact @aabiddanda
