import numpy as np
import polars as pl
from utils import *


def err_opposite_homozygotes(mat_haps, pat_haps, bafs, eps=5e-2):
    mat_geno = mat_haps.sum(axis=0)
    pat_geno = pat_haps.sum(axis=0)
    opp_homozygotes = ((mat_geno == 0) & (pat_geno == 2)) | (
        (mat_geno == 2) & (pat_geno == 0)
    )
    err_rates = (
        np.sum(np.vstack(bafs)[:, opp_homozygotes] < eps, axis=0)
        + np.sum(np.vstack(bafs)[:, opp_homozygotes] > 1 - eps, axis=0)
    ) / np.vstack(bafs).shape[0]
    err_rate_tot = np.zeros(mat_geno.size)
    err_rate_tot[opp_homozygotes] = err_rates
    return err_rate_tot


def err_het_hom(mat_haps, pat_haps, bafs, eps=5e-2):
    mat_geno = mat_haps.sum(axis=0)
    pat_geno = pat_haps.sum(axis=0)
    mat_het1 = (mat_geno == 1) & (pat_geno == 2)
    mat_het2 = (mat_geno == 1) & (pat_geno == 0)
    err_rates1 = (np.sum(np.vstack(bafs)[:, mat_het1] < eps, axis=0)) / np.vstack(
        bafs
    ).shape[0]
    err_rates2 = (np.sum(np.vstack(bafs)[:, mat_het2] > 1 - eps, axis=0)) / np.vstack(
        bafs
    ).shape[0]
    pat_het1 = (mat_geno == 2) & (pat_geno == 1)
    pat_het2 = (mat_geno == 0) & (pat_geno == 1)
    err_rates_pat1 = (np.sum(np.vstack(bafs)[:, pat_het1] < eps, axis=0)) / np.vstack(
        bafs
    ).shape[0]
    err_rates_pat2 = (
        np.sum(np.vstack(bafs)[:, pat_het2] > 1 - eps, axis=0)
    ) / np.vstack(bafs).shape[0]
    err_rate_tot = np.zeros(mat_geno.size)
    err_rate_tot[mat_het1] = err_rates1
    err_rate_tot[mat_het2] = err_rates2
    err_rate_tot[pat_het1] = err_rates_pat1
    err_rate_tot[pat_het2] = err_rates_pat2
    return err_rate_tot


if __name__ == "__main__":
    # Load in the input BAF datasets
    family_data = load_baf_data(snakemake.input["baf_pkl"])
    fps = [f for f in snakemake.input["aneuploidy_calls"]]
    aneuploidy_df = pl.concat(
        [pl.read_csv(fp, separator="\t", null_values=["NA"]) for fp in fps]
    )
    chroms = aneuploidy_df["chrom"].unique().to_numpy()
    chrom_agg = []
    pos_agg = []
    prop_err_agg = []
    for c in chroms:
        # Reading in the current datasets
        euploid_indivs = euploid_per_chrom(
            aneuploidy_df,
            mother=snakemake.params["mother_id"],
            father=snakemake.params["father_id"],
            names=family_data.keys(),
            chrom=f"{c}",
            pp_thresh=snakemake.params["ppThresh"],
        )
        mat_haps, pat_haps, bafs, real_names, pos = prep_data(
            family_dict=family_data, names=euploid_indivs, chrom=f"{c}"
        )
        if len(bafs) > 0:
            # Detecting which variants have some proportion of errors
            err_opp_hom_r = err_opposite_homozygotes(
                mat_haps, pat_haps, np.vstack(bafs), eps=snakemake.params["eps"]
            )
            err_het_hom_r = err_het_hom(
                mat_haps, pat_haps, np.vstack(bafs), eps=snakemake.params["eps"]
            )
            err_rates_x = err_opp_hom_r + err_het_hom_r
            for p, x in zip(
                pos[err_rates_x >= snakemake.params["err_rate"]],
                err_rates_x[err_rates_x >= snakemake.params["err_rate"]],
            ):
                chrom_agg.append(c)
                pos_agg.append(p)
                prop_err_agg.append(x)
    mendelian_df = pl.from_dict(
        {"chrom": chrom_agg, "pos": pos_agg, "err_rate": prop_err_agg}
    )
    mendelian_df.write_csv(snakemake.output["mendelian_tsv"], separator="\t")
