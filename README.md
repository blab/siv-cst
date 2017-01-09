# Modern-Day SIV viral diversity generated by extensive recombination and cross-species transmission

#### Sidney M. Bell<sup>1,2</sup>, Trevor Bedford<sup>1</sup>

<sup>1</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA, <sup>2</sup>Molecular and Cellular Biology Graduate Program, University of Washington, Seattle, WA, USA

## Abstract

Cross-species transmission (CST) has led to many devastating epidemics, but is still a poorly under-stood phenomenon. HIV-1 and HIV-2 (human immunodeficiency virus 1 and 2), which have collectively caused over 35 million deaths, are the result of many CSTs from chimpanzees, gorillas, and sooty mangabeys. While the immediate history of HIV is known, there are over 45 lentiviruses that infect spe-cific species of primates, and it is not understood whether host switching is a common occurrence. We thus took a phylogenetic approach to better understand the natural history of SIV recombination and CST. We used discrete trait analysis to model host species as a discrete trait in the viral phylogeny and inferred the pairwise transmission rates between each pair of 24 primate hosts. We identify 14 novel, ancient cross-species transmission events. We also find that lentiviral lineages vary widely in their abil-ity to infect new host species: SIVcol (from colobus monkeys) is evolutionarily isolated, while SIVagms (from African green monkeys) frequently move between host subspecies. We also examine the origins of SIVcpz (the predecessor of HIV-1) in greater detail than previous studies, and find that there are still large portions of the genome with unknown origins. Ultimately, these results demonstrate that while len-tiviral CST remains a rare event in the scope of evolutionary time, these viruses have a far more exten-sive history of host switching and recombination than previously described.  

## Citation

> Bell S.M., Bedford T. 2017. Modern-Day SIV viral diversity generated by extensive recombination and cross-species transmission. In prep.

## Outline

* [SIV sequence data and host metadata](data/)  
* [Analyze recombination](recombination/)  
* [Build trees for each portion of the genome](beast/main/empiricalTrees/)  
* [Infer cross-species transmissions with discrete trait analysis](beast/main/discreteTraits/)  
* Compare with a supplemental dataset: [trees](beast/supplement/empiricalTrees) and [discrete traits](beast/supplement/discreteTraits/)  
* [Figures](figures/)

## Install

Install Python packages with:

    pip install -r requirements.txt
