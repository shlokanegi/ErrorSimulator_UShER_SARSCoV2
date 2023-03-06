'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Execute line by line, or blocks of code as required
Purpose: Random blocks of commands helpful throughout the project
'''

## Install conda
# Miniconda Installer scripts - Linux : https://docs.conda.io/en/latest/miniconda.html#linux-installers
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_22.11.1-1-Linux-x86_64.sh
chmod -v +x Miniconda*.sh

## Verify HASH
# SHA256 hash - 00938c3534750a0e4069499baf8f4e6dc1c2e471c86a59caa0dd03f4a9269db6
echo "00938c3534750a0e4069499baf8f4e6dc1c2e471c86a59caa0dd03f4a9269db6 *Miniconda3-py310_22.11.1-1-Linux-x86_64.sh" | shasum --check

## Execute Installer
bash Miniconda3-py310_22.11.1-1-Linux-x86_64.sh



######### Perform tree extraction #########
cd data
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz

matUtils summary -i public-latest.all.masked.pb.gz -c clades.txt #list of all clades

## Extract a tree with of clade - 
matUtils extract -i public-latest.all.masked.pb.gz -c A.2.2 -t subtree.nwk    # A.2.2 - 214 samples
matUtils extract -i public-latest.all.masked.pb.gz -c A.1 -t subtree_A.1.nwk  # A.1 - 2542 samples
matUtils extract -i public-latest.all.masked.pb.gz -c B.1.1.1 -t subtree_B.1.1.1.nwk  # B.1.1.1 - 7163 samples


## Use PhastSim
mkdir phastSim_output
# phastSim --outpath phastSim_output/ --seed 7 --treeFile subtree.nwk --reference newref.fa\
#     --createNewick --createMAT --createFasta --createInfo --createPhylip \
#     --scale 0.333333333 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
#     --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
#     --eteFormat 1

# phastSim --outpath phastSim_output_scaled/ --seed 7 --treeFile subtree.nwk --reference newref.fa\
#     --createNewick --createMAT --createFasta --createInfo --createPhylip \
#     --scale 0.1 --alpha 1.0 --omegaAlpha 1.0 --codon \
#     --eteFormat 1


phastSim --outpath phastSimA.1_output_scaled/ --seed 7 --treeFile subtree.nwk --reference newref.fa\
    --createNewick --createMAT --createFasta --createInfo --createPhylip \
    --scale 0.000034 --alpha 1.0 --omegaAlpha 1.0 --codon \
    --eteFormat 1
conda activate usher-env; cd phastSimA.1_output_scaled; matUtils extract -i sars-cov-2_simulation_output.mat.pb -o sars-cov-2_simulation_output_collapsed.mat.pb -O
matUtils extract -i sars-cov-2_simulation_output_collapsed.mat.pb -t sars-cov-2_simulation_output_collapsed.nwk


phastSim --outpath phastSimB.1.1.1_output_scaled/ --seed 7 --treeFile subtree_B.1.1.1.nwk --reference newrefB.1.1.1.fa\
    --createNewick --createMAT --createFasta --createInfo --createPhylip \
    --scale 0.000034 --alpha 1.0 --omegaAlpha 1.0 --codon \
    --eteFormat 1
conda activate usher-env; cd phastSimB.1.1.1_output_scaled; matUtils extract -i sars-cov-2_simulation_output.mat.pb -o sars-cov-2_simulation_output_collapsed.mat.pb -O
matUtils extract -i sars-cov-2_simulation_output_collapsed.mat.pb -t sars-cov-2_simulation_output_collapsed.nwk

## From phastSim, we get outputs in ~/data/phastSim_output
## Our file of interest --> sars-cov-2_simulation_output.mat.pb


## Extracting VCF (subtree with errors)
subtree_errors.vcf
## Extract the PB MAT file
-----

## Construct new tree from scratch using USHER sampled (use VCF + empty tree)
# clade A.2.2
usher-sampled -v subtree_random_errors.vcf -t seed.nwk -d ~/data/usher_output/ -o usher_output/subtree_random_errors.pb
usher-sampled -v subtree_errors_reversions.vcf -t seed.nwk -d ~/data/usher_output_rev/ -o usher_output_rev/subtree_errors_reversions.pb
# clade A.1
usher-sampled -v subtree_random_errors.vcf -t seed.nwk -d ~/data/usherA.1_ran_output/errortree0 -o usherA.1_ran_output/errortree0/subtree_random_errors.pb
# clade B.1.1.1
usher-sampled -v subtree_errors.vcf -t seed.nwk -d ~/data/usher_output_try/ -o usher_output_try/subtree_errors.pb
cd ~/data/usher_output_try; matUtils extract -i subtree_errors.pb -o subtree_errors_collapsed.pb -O ; matUtils extract -i subtree_errors_collapsed.pb -t subtree_errors_collapsed.nwk


#### Visualize the subtree (before and after)
matUtils extract -i ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output.mat.pb -j auspice_json/subtree_phastSim_cladeA.1.json

matUtils extract -i ~/data/usherA.1_ran_output/subtree_random_errors.pb -j subtree_random_errors_cladeA.1.json
matUtils extract -i ~/data/usherA.1_rev_output/subtree_errors_reversions.pb -j subtree_errors_reversions_cladeA.1.json

matUtils extract -i ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb -j auspice_json/subtree_phastSim_cladeB.1.1.1.json
matUtils extract -i ~/data/usher_output_try/subtree_errors_collapsed.pb -j subtree_errors_cladeB.1.1.1.json



### Running errorSimulator.py script (Random errors, reversions and amplicon_dropouts)
cd data/BTE
python3 argparse_errorSimulator.v6.py \
    -t ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output.mat.pb \
    -ref ~/data/newrefA.1.fa \
    -r 0.2 -rev 0 -ad 0

python3 argparse_errorSimulator.v6.py \
    -t ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    -ref ~/data/newrefB.1.1.1.fa \
    -r 0.5 -rev 0.03 -ad 0.5

### Run calcRFD.py to calculate Robinson-Fould distance between 2 trees
python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output.nwk \
    -t2 ~/data/usherA.1_ran_output/errortree0.0/final-tree.nh

conda activate bte; python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimB.1.1.1_output_scaled/sars-cov-2_simulation_output_collapsed.nwk \
    -t2 ~/data/usher_output_try/subtree_errors_collapsed.nwk


##################################### Do the collapsing work ####################################
conda activate usher-env
# Collapse phastSim "real tree"
cd ~/data/phastSimA.1_output_scaled/
matUtils extract -i sars-cov-2_simulation_output.mat.pb -o sars-cov-2_simulation_output_collapsed.mat.pb -O
# Make a nwk file using the collapsed MAT protobuf
matUtils extract -i sars-cov-2_simulation_output_collapsed.mat.pb -t sars-cov-2_simulation_output_collapsed.mat.nwk
# Run errorSim script using collapsed file as input
conda activate bte
python3 argparse_errorSimulator.v6.py \
    -t ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.pb \
    -ref ~/data/newrefA.1.fa \
    -r 0.5 -rev 0 -ad 0 > temp/temp_out.txt
# Recreate the inferred tree using subtree_random_errors.vcf
conda activate usher-env
#usher-sampled -v subtree_random_errors.vcf -t seed.nwk -o ~/data/temp/subtree_random_errors_collapsed.pb -c -d ~/data/temp/
usher-sampled -v subtree_random_errors.vcf -t seed.nwk -o ~/data/temp/temp.pb -d ~/data/temp/
cd temp
matUtils extract -i temp.pb -o temp_collapsed.pb -O
matUtils extract -i temp_collapsed.pb -t temp_collapsed.nwk
# compare the 2 trees
conda activate bte
cd temp
python3 ~/data/calcRFD.py \
    -t1 ~/data/phastSimA.1_output_scaled/sars-cov-2_simulation_output_collapsed.mat.nwk \
    -t2 ~/data/temp/temp_collapsed.nwk > temp_rfd.txt






###################################### EXTRAS ####################################
# Extract a VCF from a .PB file
matUtils extract -i sars-cov-2_simulation_output.mat.pb -v phastSim_MAT.vcf
# Convert a .PB to a .nwk file
matUtils extract -i sars-cov-2_simulation_output.mat.pb -t sars-cov-2_simulation_output.nwk
# Make new directories with a common name and series of numbers
mkdir -p errortree{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
mkdir -p errortree{0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10}
mkdir -p errortree{0..10}


usher-sampled -v tree2.vcf -t seed.nwk -d ~/data/usher_tree2 -o usher_tree2/tree2.pb
matUtils extract -i ~/data/usher_tree1/tree1.pb -j tree1.json

python3 ~/data/calcRFD.py \
    -t1 ~/data/usher_tree2/final-tree.nh \
    -t2 ~/data/usher_tree1/final-tree.nh
