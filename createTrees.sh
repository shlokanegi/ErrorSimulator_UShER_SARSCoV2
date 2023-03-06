'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Execute line by line, or blocks of code as required
Purpose: Use Taxonium (CoV2Tree.org) to visualize trees with errors.
Input file requirements: MAT.pb (real and inferred tree), metadata.TSV
Conda environment requirements: usher-env, bte
'''

pip install taxoniumtools

# Generating error.VCF for parameters [-r 0.5 -rev 0 -ad 0]
python3 errorSimulator.py \
    -t ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    -ref ~/data/newrefB.1.1.1.fa \
    -r 0.5 -rev 0 -ad 0

# Recreate inferred collapsed tree using using error.VCF
conda activate usher-env
usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usher_output/ -o usher_output/subtree_errors.pb
cd ~/data/usher_output/
matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O
matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk

# Calculate RFD between real and inferred tree
conda activate bte
python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.nwk \
    -t2 subtree_errors_collapsed.nwk

# Convert inferred tree into JSONL format, View on taxonium
usher_to_taxonium \
    --input ~/data/usher_output/subtree_errors_collapsed.pb \
    --output ~/data/taxonium/inferredtree_ran_B.1.1.1-taxonium.jsonl \
    --metadata ~/data/taxonium/metadata_errors.tsv \
    --genbank ~/data/taxonium/hu1.gb \
    --columns strain,random_errors,reversions,amplicon_dropouts \
    --title "Inferred Tree with random error rate = 0.5, rfd = 0.08"

# Convert real tree into JSONL format, View on taxonium
usher_to_taxonium \
    --input ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    --output ~/data/taxonium/realtree_ran_B.1.1.1-taxonium.jsonl \
    --metadata ~/data/taxonium/metadata_errors.tsv \
    --genbank ~/data/taxonium/hu1.gb \
    --columns strain,random_errors,reversions,amplicon_dropouts \
    --title "Real Tree"

# Examples of colour configs
#https://taxonium.org/?config={%22colorMapping%22:{%22Yes%22:[255,0,0],%22No%22:[150,255,150]}}
#https://taxonium.org/?config={%22colorMapping%22:{%22Yes%22:[255,0,255],%22No%22:[200,200,200]}}
