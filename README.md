![A cybermutator mouse](assets/cybermutator_mouse.png)

### **`cybermutator.py`** helps you create your own cybermutator to simulate a hypermutator. 

It uses `msprime` to simulate a lineage tree of single cells in the developing animal, then uses a custom mutation simulator to add somatic mutations based on COSMIC SBS signatures.

    -h, --help            show this help message and exit
    --cells CELLS         The number of cells in your simulated sample
    --sequence SEQUENCE   Optional: the sequence you would like to simulate mutations on. Leave empty to randomly generate a sequence, or
                            query a ref genome for a sequence.
    --seq_len SEQ_LEN     Length of randomly generated sequence to simulate mutations on.
    --genome GENOME       Specify genome to generate random sequence with appropriate GC content. Currently only mm10 is supported
    --genome_fasta GENOME_FASTA
                            Path to indexed fasta to query
    --regions REGIONS [REGIONS ...]
                            region or list of regions (separated by a space) of provided reference genome to simulate mutations within, in the
                            following format: chr1:200-300
    --growth_model GROWTH_MODEL
                            currently only exponential is supported
    --coalescent_model COALESCENT_MODEL
                            Model provided to msprime.sim_ancestry
    --time TIME           Number of simulated generations / cell divisions
    --N0 N0               Founder population size (use 1 to simulate development from a single cell zygote)
    --Nt NT               Population size at time of simulated sampling.
    --seed SEED           seed for random number generation. Defaults to current time if none provided.
    --Mu MU               Overall mutation rate.
    --sbs_signatures SBS_SIGNATURES [SBS_SIGNATURES ...]
                            Path or list of paths to COSMIC SBS signatures in default tsv format. The two PolE-P286R signatures (SBS10a and
                            SBS10b) are provided.
    --sbs_weights SBS_WEIGHTS [SBS_WEIGHTS ...]
                            SBS signature weights, used if multiple signatures provided.
    --outdir OUTDIR       Directory for outputs. Will be created if it does not exist
    --name NAME           name for your simulation, will be the prefix of output files.
