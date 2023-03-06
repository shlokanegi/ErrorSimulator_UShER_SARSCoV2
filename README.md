# ErrorSimulator_UShER_SARSCoV2
A simulation platform that incorporates realistic errors into simulated SARS-CoV-2 genomes. 

## Table of Contents
* [Background and Motivation](https://github.com/shlokanegi/ErrorSimulator_UShER_SARSCoV2/blob/master/README.md#background-and-motivation)
* [Types of Errors simulated by error simulator](https://github.com/shlokanegi/ErrorSimulator_UShER_SARSCoV2/blob/master/README.md#types-of-errors-simulated-by-error-simulator)
* [Methodology](https://github.com/shlokanegi/ErrorSimulator_UShER_SARSCoV2/blob/master/README.md#methodology
)
* [Result](https://github.com/shlokanegi/ErrorSimulator_UShER_SARSCoV2/blob/master/README.md#result---usher-is-most-sensitive-to-reversions)
* [Conclusion and Future Work](https://github.com/shlokanegi/ErrorSimulator_UShER_SARSCoV2/blob/master/README.md#conclusion-and-future-work)
* [Usage](https://github.com/shlokanegi/ErrorSimulator_UShER_SARSCoV2/blob/master/README.md#usage)

## Background and Motivation
Many a times some mutations in the phylogenetic tree, seem more likely the result of contamination, sequencing errors, or any other kinds of random errors, than merely selection or recombination.

![image](https://user-images.githubusercontent.com/66521525/223062071-83cb2cef-833a-4265-939c-708aa97c01a0.png)

Above is a snapshot of a region from the SARS-CoV2 phylogeny, coloured by the AA at position 501 of the Spike Protein. 
The blue samples represent the omicron variant and the yellow ones are the variants which came before Omicron.
One thing that catches the eye are these bunch of yellow samples lying amongst the blues. 
These set of samples are more similar to the blue samples than this yellow group of samples. 
Yet, they have been colored as Yellow, probably because there were some reversions/errors at a site in their genome, which led them to have an Asparagine (N) and position 501. 

Just by looking at this, we cannot be entirely sure as to whether USHER was sensitive to such errors or not, since we don’t have a ground truth to validate against.
Therefore, this idea of developing an ErrorSimulator platform came into picture, wherein by simulating errors, we could evaluate the robustness of UShER, and then if needed, optimize it so that it’s not equally sensitive to both realistic mutations and random errors and such critical misclassifications could be avoided.

## Types of Errors simulated by error simulator
<img width="1437" alt="image" src="https://user-images.githubusercontent.com/66521525/223064389-2a8425d5-681f-41ca-8be7-fb082b304056.png">

1. **Random Errors** - These can occur through sequencing errors or errors during sample preparation. For example, sample contamination, mix-up, or mislabelling can all lead to random errors in final genome sequences. <br>

2. **Reversions** - Viral replication process is very prone to errors since the viral RNA polymerase lacks a proof-reading capability. Due to this, sometimes, a point mutation that previously occurred in the viral genome is corrected back to the original nucleotide. The likelihood of a reversion occurring depends on the total number of mutations in the sample as more the mutations, more the possibility that the RNA-pol makes a backmutation error, and this is what I’ve incorporated in my code while choosing samples to add reversions on. Also, recombination events can occur during something like co-infection, where a mixture of genetic material from different strains, can lead to reversions of previously mutated sites. 

3. **Amplicon Dropout and Contamination** - In the context of SARS-CoV-2 sequencing, PCR products from one sample can potentially contaminate another sample if there is carryover during library preparation or sequencing. For example, if a sequencing instrument is not properly cleaned between runs, PCR products from a previous run may carry over and contaminate the subsequent run. This can result in the presence of reads from another sample being misattributed to the first sample, leading to incorrect variant calls or other analysis results. Hence, adding amplicon dropout errors in the error simulator becomes very crucial

## Methodology
![image](https://user-images.githubusercontent.com/66521525/223066425-1f9a5efa-102b-4204-844d-b3f0919540ed.png)

## Result - UShER is most sensitive to reversions
![image](https://user-images.githubusercontent.com/66521525/223066594-2ae22410-1143-402b-9f6b-1f5d896e1ab9.png)

1. **Amplicon Dropouts** – No observable trend since dropouts depend on a lot of factors here. Since, I am randomly picking an amplicon from the list and then picking the source and replacement leaves at random as well, a dropout event would depend on whether that there are mutations in the source and replacement leaf, for the chosen amplicon range. Hence, we don’t see a clear trend here, and it does vary with every run.

2. **Random Errors** - Might look like a slight increasing trend, but the RFD is pretty much oscillates around this average of 0.15. It stayed the same till even 200% of error was incorporated, but then it did start to increase to 0.5 and 0.6 at much higher error percentages of 500% and 1000%.  An explanation for this might be that - For random error incorporation, we are adding more “false” mutations to leaves, which would essentially just increase the branch length. That’s probably why there isn’t a lot of rearrangement in the tree when USHER recreates it. Hence, the constant 0.15 value.

3. **Reversions** - This shows as clear increase in RFD. This kinda makes sense as here, we are actually removing mutations from the mutation list, which intuitively seems more than just reducing branch lengths. Hence, more rearrangement in the tree.![image]

## Conclusion and Future Work
* Optimization of UShER to reduce sensitivity towards errors.
* Developing an error-correction/ tree-pruning algorithm, mostly tailored to handle reversions better.
* Improving simulation of amplicon dropout errors by sampling amplicons based on a probabilistic model, exploiting the knowledge of some amplicons being more likely to dropout than others.
* Using a better metric that Robinson Fould to compare trees

## USAGE
Adding errors on a MAT (Mutation Annotated Tree) using error simulator
```
conda install -c conda-forge -c bioconda bte
```
```
python3 errorSimulator.py \
    -t tree.mat.pb \
    -ref reference.fa \
    -r 0.5 -rev 0.04 -ad 0.8
```

Calculating Robinson-Fould Distance between "real tree" and "inferred tree"
```
pip install ete3
```
```
python3 calcRFD.py \
    -t1 realtree.nwk \
    -t2 inferredtree.nwk
```
