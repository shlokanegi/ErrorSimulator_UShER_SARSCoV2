# # Run errorSimv6.py script to add random errors to our tree
# python3 argparse_errorSimulator.v6.py \
#     -t ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
#     -ref ~/data/newrefA.1.fa \
#     -r 0 -rev 0 -ad 0 > ~/data/usherA.1_rev_output/errortree10/out.txt

# conda activate usher-env; usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usherA.1_rev_output/errortree10 -o usherA.1_rev_output/errortree10/subtree_errors.pb
# cd ~/data/usherA.1_rev_output/errortree10; matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O ; matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk

# conda activate bte; python3 ~/data/calcRFD.py \
#     -t1 ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.nwk \
#     -t2 subtree_errors_collapsed.nwk > rfd.txt


## Amplicon dropout
python3 argparse_errorSimulator.v6.py \
    -t ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    -ref ~/data/newrefB.1.1.1.fa \
    -r 0 -rev 0 -ad 10 > ~/data/usherB.1.1.1_ad_output/errortree10/out.txt

conda activate usher-env; usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usherB.1.1.1_ad_output/errortree10 -o usherB.1.1.1_ad_output/errortree10/subtree_errors.pb
cd ~/data/usherB.1.1.1_ad_output/errortree10; matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O ; matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk

conda activate bte; python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.nwk \
    -t2 subtree_errors_collapsed.nwk > rfd.txt

## Reversion
python3 argparse_errorSimulator.v6.py \
    -t ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    -ref ~/data/newrefB.1.1.1.fa \
    -r 0 -rev 1 -ad 0 > ~/data/usherB.1.1.1_rev_output/errortree1.0/out.txt

conda activate usher-env; usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usherB.1.1.1_rev_output/errortree1.0 -o usherB.1.1.1_rev_output/errortree1.0/subtree_errors.pb
cd ~/data/usherB.1.1.1_rev_output/errortree1.0; matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O ; matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk

conda activate bte; python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.nwk \
    -t2 subtree_errors_collapsed.nwk > rfd.txt


## Random Errors
python3 argparse_errorSimulator.v6.py \
    -t ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    -ref ~/data/newrefB.1.1.1.fa \
    -r 1 -rev 0 -ad 0 > ~/data/usherB.1.1.1_ran_output/errortree1.0/out.txt

conda activate usher-env; usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usherB.1.1.1_ran_output/errortree1.0 -o usherB.1.1.1_ran_output/errortree1.0/subtree_errors.pb
cd ~/data/usherB.1.1.1_ran_output/errortree1.0; matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O ; matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk

conda activate bte; python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.nwk \
    -t2 subtree_errors_collapsed.nwk > rfd.txt

### ete compare ###
ete3 compare -t final-tree.nh -r ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output.nwk