import gzip as gz
import pickle
import sys
from pathlib import Path
from karyohmm import MetaHMM
import numpy as np
import polars as pl
from tqdm import tqdm


def load_baf_data(baf_pkls):
    """Load in the multiple BAF datasets."""
    family_data = {}
    for fp in baf_pkls:
        embryo_name = Path(fp).stem.split(".")[0]
        with gz.open(fp, "rb") as f:
            data = pickle.load(f)
            family_data[embryo_name] = data
    return family_data


def euploid_per_chrom(
    aneuploidy_df, mother, father, names, chrom="chr1", pp_thresh=0.95
):
    """Return only the euploid embryo names for this chromosome."""
    assert "bf_max_cat" in aneuploidy_df.columns
    assert "mother" in aneuploidy_df.columns
    assert "father" in aneuploidy_df.columns
    assert "child" in aneuploidy_df.columns
    assert "2" in aneuploidy_df.columns
    assert len(names) >= 1
    filt_names = (
        aneuploidy_df.filter(pl.col("mother") == mother)
        .filter(pl.col("father") == father)
        .filter(pl.col("child").is_in(names))
        .filter(pl.col("chrom") == chrom)
        .filter(pl.col("2") >= pp_thresh)["child"]
        .to_numpy()
    )
    if filt_names.size < 1:
        return []
    else:
        return filt_names.tolist()


def prep_data(family_dict, names, chrom="chr21", bad_snp_df=None):
    """Prepare the data for the chromosome to have the same length in BAF across all samples."""
    shared_pos = []
    for k in family_dict.keys():
        if k in names:
            shared_pos.append(family_dict[k][chrom]["pos"])
    bafs = []
    real_names = []
    mat_haps = None
    pat_haps = None
    pos = None
    try:
        collective_pos = list(set(shared_pos[0]).intersection(*shared_pos))
        for k in family_dict.keys():
            if k in names:
                cur_pos = family_dict[k][chrom]["pos"]
                baf = family_dict[k][chrom]["baf_embryo"]
                idx = np.isin(cur_pos, collective_pos)
                bafs.append(baf[idx])
                mat_haps = family_dict[k][chrom]["mat_haps"][:, idx]
                pat_haps = family_dict[k][chrom]["pat_haps"][:, idx]
                real_names.append(k)
        pos = np.sort(collective_pos)
    except IndexError:
        # NOTE: this only happens when we have all embryos as aneuploid for the chromosome.
        pass
    if (bad_snp_df is not None) and (mat_haps is not None):
        bad_chrom_pos = bad_snp_df.filter(pl.col("chrom") == chrom)["pos"].to_numpy()
        idx = np.isin(pos, bad_chrom_pos, invert=True)
        mat_haps = mat_haps[:, idx]
        pat_haps = pat_haps[:, idx]
        pos = pos[idx]
        bafs = [b[idx] for b in bafs]
    return mat_haps, pat_haps, bafs, real_names, pos


def est_genotype_quality(hmm_data, names, chrom="chr21"):
    """Estimate the genotype quality by the local estimate of disomy at each SNP across the sibling embryos."""
    assert len(hmm_data) > 0
    hmm = MetaHMM()
    shared_pos = []
    for k in hmm_data:
        if k in names:
            shared_pos.append(hmm_data[k][chrom]["pos"])
    collective_pos = list(set(shared_pos[0]).intersection(*shared_pos))
    pos = np.sort(collective_pos)
    assert pos.size > 0
    posterior_disomy_agg = []
    for k in hmm_data:
        if k in names:
            cur_pos = hmm_data[k][chrom]["pos"]
            g, karyo = hmm.marginal_posterior_karyotypes(
                hmm_data[k][chrom]["gammas"], hmm_data[k][chrom]["karyotypes"]
            )
            assert g.shape[1] == cur_pos.size
            cur_posterior_disomy = g[np.where(karyo == "2")[0], np.isin(cur_pos, pos)]
            posterior_disomy_agg.append(cur_posterior_disomy)
    posterior_disomy_ = np.vstack(posterior_disomy_agg)
    posterior_disomy = np.min(posterior_disomy_, axis=0)
    assert posterior_disomy.size == pos.size
    return posterior_disomy, pos


def extract_parameters(aneuploidy_df, mother, father, names, chrom):
    """Extract the core parameters for inference from the aneuploidy data frame."""
    filt_df = (
        aneuploidy_df.filter(pl.col("chrom") == chrom)
        .filter(pl.col("mother") == mother)
        .filter(pl.col("father") == father)
    )
    pi0_bafs = np.zeros(len(names))
    sigma_bafs = np.zeros(len(names))
    for i, n in enumerate(names):
        pi0_baf_test = filt_df.filter(pl.col("child") == n)["pi0_baf"].to_numpy()[0]
        sigma_baf_test = filt_df.filter(pl.col("child") == n)["sigma_baf"].to_numpy()[0]
        pi0_bafs[i] = pi0_baf_test
        sigma_bafs[i] = sigma_baf_test
    assert np.all(pi0_bafs > 0)
    assert np.all(sigma_bafs > 0)
    return pi0_bafs, sigma_bafs
