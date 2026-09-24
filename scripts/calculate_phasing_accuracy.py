#!/usr/bin/env python3

import sys
import time
import pandas as pd
from pathlib import Path

def compare_to_first(group):
    # Take the first GT in the block
    first_rep = group["gt_replicate"].iloc[0]
    first_truth = group["gt_truth"].iloc[0]

    # Compare every variant to the first variant
    group["replicate_match_first"] = (
        group["gt_replicate"] == first_rep
    ).map({True: "same", False: "different"})

    group["truth_match_first"] = (
        group["gt_truth"] == first_truth
    ).map({True: "same", False: "different"})

    return group
    
def check_phasing_accuracy(truth, replicate, bed_file):
    cols = [0, 1, 4, 5]
    rep = pd.read_csv(replicate, sep="\t", header=None, usecols=cols)
    truth = pd.read_csv(truth, sep="\t", header=None, usecols=cols)
    rep = rep.iloc[:, :6]
    # Name columns (adjust if needed)
    rep.columns = ["chrom", "pos", "gt_replicate", "ps_replicate"]
    truth.columns = ["chrom", "pos", "gt_truth", "ps_truth"]

    merged = pd.merge(
        rep[["chrom", "pos", "gt_replicate", "ps_replicate"]],
        truth[["chrom", "pos", "gt_truth", "ps_truth"]],
        on=["chrom", "pos"],
        how="inner"
    )

    total_common_variants = len(merged)


    rep_phased = merged["gt_replicate"].str[1] == "|"
    truth_phased = merged["gt_truth"].str[1] == "|"

    # Filter variants which are phased in both truth and replicate
    filtered = merged[rep_phased & truth_phased]
    filtered = filtered.copy()
    filtered["ps_replicate_truth"] = (
        filtered["ps_replicate"].astype(str) + "_" +
        filtered["ps_truth"].astype(str)
    )
    
    # Apply to each ps_replicate_truth group
    #filtered = filtered.groupby("ps_replicate_truth", group_keys=False).apply(compare_to_first)
    g = filtered.groupby("ps_replicate_truth")

    first_rep = g["gt_replicate"].transform("first")
    first_truth = g["gt_truth"].transform("first")

    filtered = filtered.copy()
    filtered["replicate_match_first"] = filtered["gt_replicate"] == first_rep
    filtered["truth_match_first"] = filtered["gt_truth"] == first_truth

    # Flag if truth and replicate has same phasing pattern
    filtered["match_replicate_truth"] = (
        filtered["replicate_match_first"] == filtered["truth_match_first"]
    )

    grouped_df = (
        filtered
        .groupby(["chrom", "ps_replicate_truth"])
        .filter(lambda x: len(x) >= 2)
    )

    counts_df = (
        grouped_df
        .groupby(["chrom", "ps_replicate_truth"])["match_replicate_truth"]
        .value_counts()
        .unstack(fill_value=0)
        .reset_index()
    )
    
    counts_df = counts_df.rename(columns={True: "n_match", False: "n_mismatch"})
    counts_df["correct_phasing"] = counts_df[["n_match", "n_mismatch"]].max(axis=1)
    counts_df["incorrect_phasing"] = counts_df[["n_match", "n_mismatch"]].min(axis=1)

    total_correct = counts_df["correct_phasing"].sum()
    total_incorrect = counts_df["incorrect_phasing"].sum()
    total_phased_variants = total_correct + total_incorrect

    # Fraction of correct phasing globally
    global_fraction_correct = (total_correct * 100)  / total_phased_variants

    return total_correct, total_phased_variants, global_fraction_correct

def main(manifest_replicate_truth, input_dir, output_file):
    with open(output_file, "w") as out:
        out.write("Sample\tCorrect\tTotal\tAccuracy\n")

        with open(manifest_replicate_truth) as f:
            for line in f:
                sample_id, replicate_vcf, truth_vcf = line.strip().split("\t")
                truth_sample = Path(truth_vcf).name.removesuffix(".vcf.gz")

                replicate_ps = f"{input_dir}/ps_common_variants_{sample_id}_and_truth.tsv"
                truth_ps = f"{input_dir}/ps_{truth_sample}_truth.tsv"

                print(f"[START] {sample_id}", flush=True)
                print(f"[FILES] {replicate_ps} {truth_ps}", flush=True)

                start = time.time()

                correct, total, acc = check_phasing_accuracy(
                    truth_ps, replicate_ps, None
                )

                elapsed = time.time() - start

                out.write(f"{sample_id}\t{correct}\t{total}\t{acc:.6f}\n")
                out.flush()  # ensures it's written immediately

                print(
                    f"[DONE] {sample_id} | correct={correct} total={total} acc={acc:.4f} | {elapsed:.2f}s",
                    flush=True
                )

if __name__ == "__main__":
  if len(sys.argv) != 4:
	print("Usage: script.py manifest.txt input_dir output.tsv")
	print("Run extract_phase_sets_to_calculate_phasing_accuracy.sh first to create input files")
	print("manifest.txt: sample_id replicate_vcf truth_vcf")
	sys.exit(1)

  main(sys.argv[1], sys.argv[2], sys.argv[3])
