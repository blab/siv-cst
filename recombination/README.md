##Recombination complicates the natural history of lentiviruses  
  
Lentiviruses notoriously recombine, which effectively means that different portions of the lentiviral genome have different evolutionary--and phylogenetic--histories. Thus, in order to reconstruct the history of cross-species transmission among primate lentiviruses, we had to first assess the extent of recombination _between_ SIV lineages. We used GARD to look for evidence of recombination breakpoints across a [modified version](./lanlCompendium15_lessHIV.fasta) of the LANL compendium alignment of HIV and SIVs.  
  
To ease computational intensity, we ran GARD on 3kb regions of the genome (with 1kb overlaps on either end). We repeated this with the windows offset in order to control for proximity to region edges. Alignment coordinates listed. All analyses used the 012234 nucleotide model and 3 bins of site variation.    
  
Window|Start|End|N breakpoints found    
---|---|---|---  
  
Consensus breakpoints:  
  
We then split the [full dataset alignments](../data/) along these breakpoints (all alignments held the compendium fixed and discarded insertions, coordinates are the same) and built maximum likelihood trees with RAxML for each segment. Tracing each taxa across segments allows us to visualize how each taxa's phylogenetic placement varies across breakpoints.  
  
![Recombination summary](https://github.com/blab/siv-cst/blob/master/figures/png/Fig1.png)
