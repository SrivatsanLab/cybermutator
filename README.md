![A cybermutator mouse](assets/cybermutator_mouse.png)

### **`cybermutator.py`** helps you create your own cybermutator to simulate a hypermutator. 

It uses `msprime` to simulate a lineage tree of single cells in the developing animal, then uses a custom mutation simulator to add somatic mutations based on COSMIC SBS signatures.

    options:
    -h, --help            show this help message and exit
    --genome GENOME
    --genome_fasta GENOME_FASTA
    --regions REGIONS [REGIONS ...]
    --sbs_signatures SBS_SIGNATURES [SBS_SIGNATURES ...]
    --sbs_weights SBS_WEIGHTS [SBS_WEIGHTS ...]
    --outdir OUTDIR
    --name NAME
    --Mu MU
    --N0 N0
    --Nt NT
    --time TIME
    --seed SEED