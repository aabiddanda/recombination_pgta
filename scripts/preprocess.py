import gzip as gz
import pickle
import sys

from functools import reduce
import numpy as np
import polars as pl
from tqdm import tqdm

if __name__ == "__main__":
    meta_dict = {}
    mat_df = pl.read_csv(
        snakemake.input["maternal_data"], separator="\t", null_values=["NA"]
    )
    pat_df = pl.read_csv(
        snakemake.input["paternal_data"], separator="\t", null_values=["NA"]
    )
    embryo_df = pl.read_csv(
        snakemake.input["embryo_data"], separator="\t", null_values=["NA"]
    )
    mat_chroms = mat_df["chrom"].unique().to_numpy()
    pat_chroms = pat_df["chrom"].unique().to_numpy()
    embryo_chroms = embryo_df["chrom"].unique().to_numpy()
    chroms = reduce(np.intersect1d, (mat_chroms, pat_chroms, embryo_chroms))
    for c in chroms:
        print(f"Processing {c}...", file=sys.stderr)
        mat_chrom_df = mat_df.filter(pl.col("chrom") == c)[
            ["chrom", "pos", "mat0", "mat1"]
        ]
        pat_chrom_df = pat_df.filter(pl.col("chrom") == c)[
            ["chrom", "pos", "pat0", "pat1"]
        ]
        embryo_chrom_df = embryo_df.filter(pl.col("chrom") == c)["chrom", "pos", "baf"]
        merged_df = mat_chrom_df.join(pat_chrom_df, on=["chrom", "pos"]).join(
            embryo_chrom_df, on=["chrom", "pos"]
        )
        merged_df = merged_df.sort(pl.col("pos"))
        mat_haps = np.vstack(
            [merged_df["mat0"].to_numpy(), merged_df["mat1"].to_numpy()]
        )
        pat_haps = np.vstack(
            [merged_df["pat0"].to_numpy(), merged_df["pat1"].to_numpy()]
        )
        chrom_dict = {
            "baf_embryo": merged_df["baf"].to_numpy(),
            "pos": merged_df["pos"].to_numpy(),
            "mat_haps": mat_haps,
            "pat_haps": pat_haps,
            "aploid": "real_data",
        }
        meta_dict[c] = chrom_dict
    pickle.dump(meta_dict, gz.open(snakemake.output["baf_pkl"], "wb"))
