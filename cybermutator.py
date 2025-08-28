#!/usr/bin/env python3
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from Bio.Seq import Seq
from dataclasses import dataclass
from typing import List, Tuple
import argparse
import os
from tqdm import tqdm
sns.set_theme(style="ticks")

## TODO:
#       - add additional SBS signatures
#       - test with multiple SBS signatures
#       - genomic coordinates in VAF dataframe

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
    regions = [parse_region(reg) for reg in regions]
    for region in regions:
        seq = ref.fetch(region.chrom, region.start, region.end).upper()
        sequence_parts.append(seq)
        coord_map.extend((region.chrom, i) for i in range(region.start, region.end))
    full_sequence = "".join(sequence_parts)
    ref.close()
    return full_sequence, coord_map

def load_signatures(paths: List[str], weights: List[float], genome: str) -> pd.DataFrame:
    if len(paths) != len(weights):
        raise ValueError("Number of signatures and weights must match")
    if not np.isclose(sum(weights), 1.0):
        raise ValueError("Signature weights must sum to 1.0")
    combined = None

    for path, weight in zip(paths, weights):
        print(path, weight)
        # Infer signature name from filename (e.g., "v3.3_SBS10a_PROFILE.txt" → "SBS10a")
        sig_name = f'{os.path.basename(path).split("_")[1]}_{genome}'  # robust if formatted as SBS10a_PROFILE.txt
        df = pd.read_csv(path, sep="\t", index_col=0)

        # Extract ancestral and derived context (e.g., T[C>T]T → anc: TCT, der: TAT)
        df[['anc', 'der']] = df.apply(lambda row: pd.Series([
            f'{row.name[0]}{row.name[2]}{row.name[-1]}',  # ancestral trinuc
            f'{row.name[0]}{row.name[4]}{row.name[-1]}'   # derived trinuc
        ]), axis=1)

        # Pivot to context × mutation matrix
        sig_matrix = df.pivot_table(index='anc', columns='der', values=sig_name, fill_value=0)

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
        trinucs[i-1] = context
        tables.sites.add_row(position=i-1, ancestral_state=context)
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
                mut_probs = mut_probs[mut_probs > 0]  # filter zero-prob alleles
                if mut_probs.empty:
                    continue  # no valid mutations from this context
                der = np.random.choice(mut_probs.index, p=mut_probs / mut_probs.sum())
                mutations.append((site_index, node, der))
                
    mutations.sort(key=lambda x: (x[0], -ts.node(x[1]).time))
    tables.mutations.clear()
    for site, node, derived in mutations:
        tables.mutations.add_row(site=site, node=node, derived_state=derived)
    return tables.tree_sequence()

def compute_VAF(ts, trinucs, rep=1, coordmap = None):
    '''
    Computes variant allele frequency of all simulated mutations.
    Returns a dataframe with columns:
    "rep", "site", "anc", "anc_count", "der", "der_count", "VAF"
    where 
    rep = simulation replicate
    site = numerical position in provided or randomly generated sequence, or genomic coordinates in reference genome
    anc = "ancestral" allele - the allele at that site in the sequence or reference genome
    anc_count = count of ancestral alleles in the simulation
    der = derived allele - the simulated mutation
    der_count = count of derived allele in the simulation 
    VAF = der_count / (count of all alleles at that site)
    '''

    records = []
    for n,var in enumerate(ts.variants()):
        if coordmap:
            pos = var.site.position
            anc = trinucs[pos]
            site = f'{coordmap[int(pos)][0]}:{coordmap[int(pos)][1]}' #if coordmap provided, include genomic coordinates in output rather than arbitrary index
        else:
            site = var.site.position
            anc = trinucs[site]
        allele_counts = np.bincount(var.genotypes, minlength=len(var.alleles))
        total_alleles = allele_counts.sum()
        der = [a for a in var.alleles if a != anc]
        # print(anc, der)
        if len(der) == 0:
            continue
        else:
            for a in der:
                der_count = np.sum(var.genotypes == var.alleles.index(a))
                anc_count = total_alleles - der_count
                vaf = der_count / total_alleles
                records.append([
                    rep,
                    site,
                    anc,
                    anc_count,
                    a,
                    der_count,
                    vaf
                ])
    return pd.DataFrame(
        records,
        columns=["rep", "site", "anc", "anc_count", "der", "der_count", "VAF"]
    )

# def plot_VAF(vafs, bins=25, outpath=None):
#     '''
#     Takes the VAF column of the dataframe from compute_VAF and plots a histogram.
#     '''
#     plt.figure()
#     plt.hist(vafs, bins=bins, range=(0, 1), edgecolor='black')
#     plt.xlabel("Variant Allele Frequency (VAF)")
#     plt.ylabel("Count")
#     plt.tight_layout()
#     if outpath:
#         plt.savefig(outpath)

def plot_VAF(vafs, name=None, outpath=None):
    f, ax = plt.subplots(figsize=(7, 5))
    sns.despine(f)
    sns.histplot(vafs['VAF'], 
                # binrange=(0,0.1), 
                log_scale=(True,False),
                bins=50, ax=ax)
    plt.xlabel("VAF", fontweight='bold', fontsize=12)
    plt.ylabel("Count of Sites", fontweight='bold', fontsize=12)
    if name:
        plt.title(name,fontweight='bold', fontsize=12)
    if outpath:
        f.savefig(outpath, dpi=600, bbox_inches='tight')

def plot_counts(mut_counts, name=None, outpath=None):

    f, ax = plt.subplots(figsize=(7, 5))
    sns.despine(f)
    sns.histplot(mut_counts, bins=max(mut_counts), ax=ax)
    # ax.set_yscale('log')
    plt.xlabel("Count of Mutations", fontweight='bold', fontsize=12)
    plt.ylabel("Count of Replicates", fontweight='bold', fontsize=12)
    if name:
        plt.title(name,fontweight='bold', fontsize=12)
    if outpath:
        f.savefig(outpath, dpi=600, bbox_inches='tight')

def get_spectra(vafs):
    types = ['ACA>AAA', 'ACC>AAC', 'ACG>AAG', 'ACT>AAT', 'ACA>AGA', 'ACC>AGC',
             'ACG>AGG', 'ACT>AGT', 'ACA>ATA', 'ACC>ATC', 'ACG>ATG', 'ACT>ATT',
             'ATA>AAA', 'ATC>AAC', 'ATG>AAG', 'ATT>AAT', 'ATA>ACA', 'ATC>ACC',
             'ATG>ACG', 'ATT>ACT', 'ATA>AGA', 'ATC>AGC', 'ATG>AGG', 'ATT>AGT',
             'CCA>CAA', 'CCC>CAC', 'CCG>CAG', 'CCT>CAT', 'CCA>CGA', 'CCC>CGC',
             'CCG>CGG', 'CCT>CGT', 'CCA>CTA', 'CCC>CTC', 'CCG>CTG', 'CCT>CTT',
             'CTA>CAA', 'CTC>CAC', 'CTG>CAG', 'CTT>CAT', 'CTA>CCA', 'CTC>CCC',
             'CTG>CCG', 'CTT>CCT', 'CTA>CGA', 'CTC>CGC', 'CTG>CGG', 'CTT>CGT',
             'GCA>GAA', 'GCC>GAC', 'GCG>GAG', 'GCT>GAT', 'GCA>GGA', 'GCC>GGC',
             'GCG>GGG', 'GCT>GGT', 'GCA>GTA', 'GCC>GTC', 'GCG>GTG', 'GCT>GTT',
             'GTA>GAA', 'GTC>GAC', 'GTG>GAG', 'GTT>GAT', 'GTA>GCA', 'GTC>GCC',
             'GTG>GCG', 'GTT>GCT', 'GTA>GGA', 'GTC>GGC', 'GTG>GGG', 'GTT>GGT',
             'TCA>TAA', 'TCC>TAC', 'TCG>TAG', 'TCT>TAT', 'TCA>TGA', 'TCC>TGC',
             'TCG>TGG', 'TCT>TGT', 'TCA>TTA', 'TCC>TTC', 'TCG>TTG', 'TCT>TTT',
             'TTA>TAA', 'TTC>TAC', 'TTG>TAG', 'TTT>TAT', 'TTA>TCA', 'TTC>TCC',
             'TTG>TCG', 'TTT>TCT', 'TTA>TGA', 'TTC>TGC', 'TTG>TGG', 'TTT>TGT']
    
    vafs['type']=vafs['anc']+'>'+vafs['der']
    spectra = vafs['type'].value_counts()
    spectra = spectra.reindex(types, fill_value=0).to_frame(name="count")
    spectra = spectra.T
    return spectra

def standardize_context(context):
    ref = context[1]
    if ref in {"A", "G"}:
        # Reverse complement both sides
        tri_from, tri_to = context.split(">")
        tri_from_c = str(Seq(seq).complement())
        tri_to_c = str(Seq(seq).complement())
        return f"{tri_from_c}>{tri_to_c}"
    else:
        return context

def plot_spectra(spectra, name=None, outpath=None):

    spectra.columns = [standardize_context(mut) for mut in spectra.columns.to_list()]
    
    titles = [
        'C>A',
        'C>G',
        'C>T',
        'T>A',
        'T>C',
        'T>G'
    ]
    mut_types = {t:[] for t in titles}

    for mut in spectra.columns.to_list():
        # print(mut)
        # C/G mutations
        if mut[1]=='C' and mut[-2]=='A':
            mut_types['C>A'].append(mut)
        elif mut[1]=='C' and mut[-2]=='G':
            mut_types['C>G'].append(mut)
        elif mut[1]=='C' and mut[-2]=='T':
            mut_types['C>T'].append(mut)

        # A/T mutations
        elif mut[1]=='T' and mut[-2]=='A':
            mut_types['T>A'].append(mut)
        elif mut[1]=='T' and mut[-2]=='C':
            mut_types['T>C'].append(mut)
        elif mut[1]=='T' and mut[-2]=='G':
            mut_types['T>G'].append(mut)

        hex_cols = [
            '#438CFD',
            '#01182E',
            '#ff1A5E',
            '#E6E6E6',
            '#80F15D',
            '#FFDB58'
        ]

    data = spectra.iloc[0]
    fig, axes = plt.subplots(1,6,figsize=(24,6), sharey=True, gridspec_kw={'wspace': 0})
    for i, ax in enumerate(axes):
        x = mut_types[titles[i]]
        ax.bar(x, data[x], color=hex_cols[i])
        ax.set_title(f"{titles[i]}",
                     fontweight='bold', 
                     fontsize=16
                    )
        ax.set_xticklabels(x,rotation=90)

        rect = patches.Rectangle(
            (0.01, 0.95),  # Bottom left corner of the rectangle (relative to axis coordinates)
            0.98, 0.05,     # Width and height (relative to axis coordinates)
            transform=ax.transAxes,  # Use axis coordinates
            color=hex_cols[i],         # Color of the patch
            clip_on=False            # Ensure the patch extends outside the axis
        )
        ax.add_patch(rect)
        ax.spines["top"].set_visible(False)
        ax.set_ylim(max(data))
        if i == 0:
            # Add y-axis label only to the first subplot
            ax.set_ylabel("Mutations", fontweight='bold', fontsize=16)
        else:
            # Remove the y-axis line and ticks for other subplots
            ax.yaxis.set_visible(False)
            ax.spines["left"].set_visible(False)

    axes[0].set_ylim(0, max(data)+0.1*max(data))
    # axes[0].set_ylim(0, 350)
    if name:
        fig.suptitle(name,fontweight='bold', fontsize=20)
    if outpath:
        fig.savefig(outpath, dpi=600, bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser(description='Create your own cybermutator to simulate a hypermutator.')
    parser.add_argument("--cells", type=float, default=1000,
                        help='Float. The number of cells in your simulated sample')
    parser.add_argument("--reps", type=int, default=100,
                        help='Integer. Number of simulation replicates')
    parser.add_argument("--sequence", type=str,
                        help='Optional: the sequence you would like to simulate mutations on. Leave empty to randomly generate a sequence, or query a ref genome for a sequence.')
    parser.add_argument("--seq_len", type=int,
                        help='Integer. Length of randomly generated sequence to simulate mutations on.')
    parser.add_argument("--new_seq_per_rep", action='store_true',
                        help='Include flag if you want to generate a new random sequence each replicate. Overrides other seq arguments.')
    parser.add_argument("--genome", type=str, default="mm10",
                        help='Specify genome to generate random sequence with appropriate GC content. Currently only mm10 is supported')
    parser.add_argument("--genome_fasta", type=str, default='/shared/biodata/reference/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa',
                        help='Path to indexed fasta to query ')
    parser.add_argument("--regions", nargs="+", type=parse_region,
                        help='region or list of regions (separated by spaces) of provided reference genome to simulate mutations within, in the following format: chr1:200-300')
    parser.add_argument("--growth_model", type=str, default="exponential",
                        help='currently only exponential is supported')
    parser.add_argument("--coalescent_model", type=str, default="hudson",
                        help='Model provided to msprime.sim_ancestry')
    parser.add_argument("--Ne", type=float, default=1e6,
                        help='Integery. Population size at time of simulated sampling.')
    parser.add_argument("--Mu", type=float, default=2e-6,
                        help='Float. Overall mutation rate. default=2e-6')
    parser.add_argument("--sbs_signatures", nargs="+", default=["cybermutator/SBS/v3.3_SBS5_PROFILE.txt","cybermutator/SBS/v3.3_SBS10a_PROFILE.txt", "cybermutator/SBS/v3.3_SBS10b_PROFILE.txt"],
                        help='Path or list of paths (separated by spaces) to COSMIC SBS signatures in default tsv format. The two PolE-P286R signatures (SBS10a and SBS10b) and SBS5 are provided and used by default.')
    parser.add_argument("--sbs_weights", nargs="+", type=float, default=[0.5, 0.25, 0.25],
                        help='floats separated by spaces. SBS signature weights, used if multiple signatures provided.')
    parser.add_argument("--outdir", type=str, default='cybermutator_results',
                        help='Directory for outputs. Will be created if it does not exist')
    parser.add_argument("--name", type=str, default="cybermutator",
                        help='name for your simulation, will be the prefix of output files.')
    parser.add_argument("--save", type=bool, default=True,
                        help='True or False. Whether to save outputs. Default = True')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    nucleotides = list("ACGT")
    probs = {"mm10": [0.21, 0.21, 0.29, 0.29]}  # default GC content

    if args.sequence:
        sequence = args.sequence
    elif args.regions and args.genome_fasta:
        sequence, coordmap = load_reference_sequence(args.genome_fasta, args.regions)
    elif args.seq_len:
        p = probs[args.genome]
        sequence = "".join(np.random.choice(nucleotides, p=p, size=args.seq_len + 2))
    else:
        raise ValueError("Must provide either --sequence or --regions + --genome_fasta or --seq_len")

    if args.growth_model == 'exponential':
        growth_rate = np.log(2)
        demography = msprime.Demography()
        demography.add_population(
            initial_size=args.Ne,
            growth_rate=growth_rate
        )

    vafs = []
    for rep in tqdm(range(0,args.reps)):
        if args.new_seq_per_rep:
            if not args.seq_len:
                raise ValueError("Must provide --seq_len with --new_seq_per_rep")
            p = probs[args.genome]
            sequence = "".join(np.random.choice(nucleotides, p=p, size=args.seq_len + 2))
        
        ts = msprime.sim_ancestry(
            samples=args.cells,
            model=args.coalescent_model,
            demography=demography,
            sequence_length=len(sequence),
            recombination_rate=0,
            ploidy=2,
        )
        ts, trinucs = add_context(ts, sequence)
        transition_matrix = load_signatures(args.sbs_signatures, args.sbs_weights, args.genome)
        ts_mut = simulate_mutations(ts, sequence, trinucs, transition_matrix, mu=2e-6)

        if args.regions != None:
            vaf_df = compute_VAF(ts_mut, trinucs, coordmap=coordmap, rep=rep)
        else:
            vaf_df = compute_VAF(ts_mut, trinucs, rep=r)
        vafs.append(vaf_df)

    vafs = pd.concat(vafs)
    spectra = get_spectra(vafs)

    if args.save == True:
        vaf_csv_path = os.path.join(args.outdir, f"{args.name}.csv")
        vaf_png_path = os.path.join(args.outdir, f"{args.name}_vaf.png")
        mut_png_path = os.path.join(args.outdir, f"{args.name}_mut.png")
        spectra_png_path = os.path.join(args.outdir, f"{args.name}_spectra.png")
        vafs.to_csv(vaf_csv_path, index=False)
        plot_VAF(vafs["VAF"], name=args.name, outpath=vaf_png_path)
        plot_counts(vafs['rep'].value_counts(), name=args.name, outpath=mut_png_path)
        plot_spectra(spectra, name=args.name, outpath=spectra_png_path)

if __name__ == "__main__":
    main()