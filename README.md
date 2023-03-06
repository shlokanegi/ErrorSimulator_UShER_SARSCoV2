# ErrorSimulator_UShER_SARSCoV2
A simulation platform that incorporates realistic errors into simulated SARS-CoV-2 genomes. 

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
<img width="721" alt="image" src="https://user-images.githubusercontent.com/66521525/223063296-97c4ec4d-43fd-4eba-9aee-dc7696314418.png"> <img width="723" alt="image" src="https://user-images.githubusercontent.com/66521525/223063365-bc579d56-100f-44a3-9be7-11534e1434d1.png">



## Types of Errors incorporated
