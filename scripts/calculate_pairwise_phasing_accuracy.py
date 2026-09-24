#!/usr/bin/env python3

import sys
import time
import pandas as pd
from itertools import combinations
from pathlib import Path

def same_haplotype_fast(a, b):
    return a[0] == b[0]  # assumes "A|B"

def calculate_pairwise_phasing_accuracy(truth_data, replicate_data):

    cols = [0, 1, 4, 5]

    rep = pd.read_csv(replicate_data, sep="\t", header=None, usecols=cols)
    truth = pd.read_csv(truth_data, sep="\t", header=None, usecols=cols)

    rep.columns = ["chrom", "pos", "gt_rep", "ps_rep"]
    truth.columns = ["chrom", "pos", "gt_truth", "ps_truth"]

    merged = rep.merge(truth, on=["chrom", "pos"], suffixes=("_rep", "_truth"))

    correct_pairs = 0
    total_pairs = 0

    for _, group in merged.groupby("ps_rep", sort=False):

        if len(group) < 2:
            continue

        gts_rep = group["gt_rep"].to_numpy()
        gts_truth = group["gt_truth"].to_numpy()
        ps_truth = group["ps_truth"].to_numpy()

        for i, j in combinations(range(len(group)), 2):

            if ps_truth[i] != ps_truth[j]:
                continue

            truth_rel = gts_truth[i].split("|")[0] == gts_truth[j].split("|")[0]
            rep_rel   = gts_rep[i].split("|")[0] == gts_rep[j].split("|")[0]

            correct_pairs += (truth_rel == rep_rel)
            total_pairs += 1

    accuracy = correct_pairs / total_pairs if total_pairs else 0
    return correct_pairs, total_pairs, accuracy

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

                correct, total, acc = calculate_pairwise_phasing_accuracy(
                    truth_ps, replicate_ps
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
