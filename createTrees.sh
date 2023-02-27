### Use Taxonium to visualize trees with errors
pip install taxoniumtools

usher_to_taxonium \
    --input ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    --output ~/data/taxonium/realtreeA.1-taxonium.jsonl \
    --metadata metadata_random_errors.tsv \
    --columns strain,random_errors,reversions,amplicon_dropouts




#-r 0.5 -rev 0.03 -ad 0.5
conda activate usher-env
usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usher_output -o ~/data/usher_output/subtree_errors.pb
cd ~/data/usher_output
matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O ; matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk
conda activate bte
usher_to_taxonium \
    --input ~/data/usher_output/subtree_errors_collapsed.pb \
    --output ~/data/taxonium/inferredtreeA.1-taxonium.jsonl \
    --metadata ~/data/taxonium/metadata_random_errors.tsv \
    --columns strain,random_errors,reversions,amplicon_dropouts


usher_to_taxonium \
    --input ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    --output ~/data/taxonium/realtreeB.1.1.1-taxonium.jsonl \
    --metadata metadata_errors.tsv \
    --genbank hu1.gb \
    --columns strain,random_errors,reversions,amplicon_dropouts \
    --title "Real Tree" \
    --config_json https://cov2tree.org/?config={"colorMapping":{"Yes":[255,0,0],"No":[0,255,0]}}

usher_to_taxonium \
    --input ~/data/usher_output_try/subtree_errors_collapsed.pb \
    --output ~/data/taxonium/inferredtreeB.1.1.1-taxonium.jsonl \
    --metadata ~/data/taxonium/metadata_errors.tsv \
    --genbank hu1.gb \
    --columns strain,random_errors,reversions,amplicon_dropouts \
    --title "Inferred Tree with errors (-r 0.5 -rev 0.03 -ad 0.5)"


https://taxonium.org/?config={%22colorMapping%22:{%22Yes%22:[255,0,0],%22No%22:[150,255,150]}}

https://taxonium.org/?config={%22colorMapping%22:{%22Yes%22:[0,0,255],%22No%22:[200,200,200]}}
