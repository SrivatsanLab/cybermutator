#!/usr/bin/env python3
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from Bio.Seq import Seq
from dataclasses import dataclass
from typing import List, Tuple
import argparse
import os
import time

## TODO:
#       - add additional SBS signatures
#       - test with multiple SBS signatures

@dataclass
class GenomicRegion:
    chrom: str
    start: int
    end: int

def parse_region(region_str):
    chrom, coords = region_str.split(":")
    start, end = map(int, coords.split("-"))
    return GenomicRegion(chrom, start, end)

def load_reference_sequence(fasta_path: str, regions: List[GenomicRegion], min_context=1) -> Tuple[str, List[Tuple[str, int]]]:
    ref = pysam.FastaFile(fasta_path)
    sequence_parts = []
    coord_map = []
    for region in regions:
        seq = ref.fetch(region.chrom, region.start, region.end).upper()
        sequence_parts.append(seq)
        coord_map.extend((region.chrom, i) for i in range(region.start, region.end))
    full_sequence = "".join(sequence_parts)
    ref.close()
    if len(full_sequence) != len(coord_map):
        raise ValueError("Length mismatch between sequence and coordinate map.")
    return full_sequence, coord_map

def load_signatures(paths: List[str], weights: List[float], genome: str) -> pd.DataFrame:
    if len(paths) != len(weights):
        raise ValueError("Number of signatures and weights must match")
    if not np.isclose(sum(weights), 1.0):
        raise ValueError("Signature weights must sum to 1.0")
    combined = None
    for path, weight in zip(paths, weights):
        df = pd.read_csv(path, sep="\t", index_col=0)
        df[['anc', 'der']] = df.apply(lambda row: pd.Series([f'{row.name[0]}{row.name[2]}{row.name[-1]}',
                                                             f'{row.name[0]}{row.name[4]}{row.name[-1]}']), axis=1)
        sig_col = [col for col in df.columns if col.startswith("SBS")][-1]
        sig_matrix = df.pivot_table(index='anc', columns='der', values=sig_col, fill_value=0)
        sig_matrix = sig_matrix.reindex(index=sig_matrix.index.union(sig_matrix.columns),
                                        columns=sig_matrix.index.union(sig_matrix.columns),
                                        fill_value=0)
        sig_matrix *= weight
        if combined is None:
            combined = sig_matrix
        else:
            combined = combined.add(sig_matrix, fill_value=0)
    return combined

def add_context(ts, sequence):
    tables = ts.dump_tables()
    tables.sites.clear()
    trinucs = {}
    for i in range(1, int(ts.sequence_length) - 1):
        context = sequence[i - 1:i + 2]
        if context[1] in 'GA':
            context = str(Seq(context).complement())
        trinucs[i] = context
        tables.sites.add_row(position=i, ancestral_state=context[1])
    return tables.tree_sequence(), trinucs

def simulate_mutations(ts, sequence, trinucs, transition_matrix, mu):
    site_probs = np.array([transition_matrix.loc[trinucs[i]].sum() if trinucs[i] in transition_matrix.index else 0 for i in trinucs])
    seq_prob = site_probs.sum()
    mutations = []
    tables = ts.dump_tables()
    for tree in ts.trees():
        for node in tree.nodes():
            branch_len = tree.branch_length(node)
            if branch_len == 0.0:
                continue
            expected_mut = mu * seq_prob * branch_len
            n = np.random.poisson(expected_mut)
            for _ in range(n):
                site_index = np.random.choice(list(trinucs.keys()), p=site_probs / seq_prob)
                mut_probs = transition_matrix.loc[trinucs[site_index]]
                der = np.random.choice(mut_probs.index, p=mut_probs / mut_probs.sum())
                tables.mutations.add_row(site=site_index, node=node, derived_state=der[1])
    return tables.tree_sequence()

def compute_VAF(ts):
    records = []
    for var in ts.variants():
        site = var.site.position
        anc = var.site.ancestral_state
        allele_counts = np.bincount(var.genotypes, minlength=len(var.alleles))
        total_alleles = allele_counts.sum()
        der = [a for a in var.alleles if a != anc]
        if len(der) == 0:
            continue
        else: 
            for a in der:
                der_count = np.sum(var.genotypes == var.alleles.index(a))
                anc_count = total_alleles - der_count
                vaf = der_count / total_alleles
                records.append([
                    site,
                    anc,
                    anc_count,
                    a,
                    der_count,
                    vaf
                ])
    return pd.DataFrame(
        records,
        columns=["site", "anc", "anc_count", "der", "der_count", "VAF"]
    )

def plot_VAF(vafs, outpath):
    plt.figure()
    plt.hist(vafs, bins=25, range=(0, 1), edgecolor='black')
    plt.xlabel("Variant Allele Frequency (VAF)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(outpath)

def main():
    parser = argparse.ArgumentParser(description='Create your own cybermutator to simulate a hypermutator.')
    parser.add_argument("--cells", type=int, default=1000,
                        help='The number of cells in your simulated sample')
    parser.add_argument("--sequence", type=str,
                        help='Optional: the sequence you would like to simulate mutations on. Leave empty to randomly generate a sequence, or query a ref genome for a sequence.')
    parser.add_argument("--seq_len", type=int,
                        help='Length of randomly generated sequence to simulate mutations on.')
    parser.add_argument("--genome", type=str, default="mm10",
                        help='Specify genome to generate random sequence with appropriate GC content. Currently only mm10 is supported')
    parser.add_argument("--genome_fasta", type=str,
                        help='Path to indexed fasta to query ')
    parser.add_argument("--regions", nargs="+", type=parse_region,
                        help='region or list of regions (separated by a space) of provided reference genome to simulate mutations within, in the following format: chr1:200-300')
    parser.add_argument("--growth_model", type=str, default="exponential",
                        help='currently only exponential is supported')
    parser.add_argument("--coalescent_model", type=str, default="hudson",
                        help='Model provided to msprime.sim_ancestry')
    parser.add_argument("--time", type=int, default=30,
                        help='Number of simulated generations / cell divisions')
    parser.add_argument("--N0", type=int, default=1,
                        help='Founder population size (use 1 to simulate development from a single cell zygote)')
    parser.add_argument("--Nt", type=int, default=10000,
                        help='Population size at time of simulated sampling.')
    parser.add_argument("--seed", type=float,
                        help='seed for random number generation. Defaults to current time if none provided.')
    parser.add_argument("--Mu", type=float, default=2e-6,
                        help='Overall mutation rate.')
    parser.add_argument("--sbs_signatures", nargs="+", default=["cybermutator/SBS/v3.3_SBS10a_PROFILE.txt", "cybermutator/SBS/v3.3_SBS10b_PROFILE.txt"],
                        help='Path or list of paths to COSMIC SBS signatures in default tsv format. The two PolE-P286R signatures (SBS10a and SBS10b) are provided.')
    parser.add_argument("--sbs_weights", nargs="+", type=float, default=[0.5, 0.5],
                        help='SBS signature weights, used if multiple signatures provided.')
    parser.add_argument("--outdir", type=str, default='cybermutator_results',
                        help='Directory for outputs. Will be created if it does not exist')
    parser.add_argument("--name", type=str, default="cybermutator",
                        help='name for your simulation, will be the prefix of output files.')
    args = parser.parse_args()

    if args.seed is None:
        args.seed = time.time()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    nucleotides = list("ACGT")
    probs = {"mm10": [0.21, 0.21, 0.29, 0.29]}  # default GC content

    if args.sequence:
        sequence = args.sequence
    elif args.regions and args.genome_fasta:
        sequence, coordmap = load_reference_sequence(args.genome_fasta, args.regions)
    elif args.seq_len:
        p = probs.get(args.genome, [0.25]*4)
        sequence = "".join(np.random.choice(nucleotides, p=p, size=args.seq_len + 2))
    else:
        raise ValueError("Must provide either --sequence or --regions + --genome_fasta or --seq_len")

    ts = msprime.sim_ancestry(samples=args.cells, sequence_length=len(sequence), recombination_rate=0, ploidy=2, random_seed=args.seed)
    ts, trinucs = add_context(ts, sequence)
    transition_matrix = load_signatures(args.sbs_signatures, args.sbs_weights, args.genome)
    ts_mut = simulate_mutations(ts, sequence, trinucs, transition_matrix, mu=args.Mu)

    vaf_df = compute_VAF(ts_mut)
    vaf_csv_path = os.path.join(args.outdir, f"{args.name}.csv")
    vaf_png_path = os.path.join(args.outdir, f"{args.name}.png")
    vaf_df.to_csv(vaf_csv_path, index=False)
    plot_VAF(vaf_df["VAF"], vaf_png_path)

if __name__ == "__main__":
    main()