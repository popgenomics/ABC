# simple pipeline to perform simulations for downstream ABC analysis  
## requirement:  
### python:  
Be sure that you have **python** and **pypy** in your /usr/bin directory  
The pipeline also works with only the standard **python** interpreter, but is much slower.  
Comparison using the **time** command for the simulation of **1,000** multilocus datsets of **100** loci:  
**pypy** : real = 1m34.698s   
**python** : real = 21m24.772s   
  
If you do not have **pypy** installed and still want to use this pipeline despite the extra computation time, you just have to change the first line of **mscalc.py** by replacing **#!/usr/ bin/pypy** by **#!/usr/bin/python**  
  
Finaly, **priorgen.py** imports the **numpy** library.  
  
### models: 
Four categories of models are allowed:  
SI = **S**trict **I**solation model (split of an ancestral population in two isolated daughter populations)  
AM = **A**ncient **M**igration model. Like SI, but with migration between daughter populations between the time of spli **Tsplit** and **Tam**. There is then no migration between **Tam** and the present time.  
SC = **S**econdary **C**ontact model. Like SI, but with no migration between daughter populations between the time of split **Tsplit** and **Tsc**. There is then a secondary contact with continuous migration between **Tsc** and the present time.  
IM = **I**solation with **M**igration model. Like SI but with gene flow between **Tsplit** and the present time.  

### submodels:  
Two models of effective population size  (**Ne**):  
1N = all loci share the same **Ne**  
2N = two categories of loci with a neutral **Ne** (proportion **P**) and selectively-linked loci (**1-P**) with a size equal to **bf.Ne** where **bf** is a background factor reducing the **Ne** of a locus_i similarly in the 3 populations (ancestral, pop_A and pop_B).  
  
1M = all loci share the same migration rate **M**.  
2M = two categories of loci with a neutral **M** (>0) or loci linked to species barriers (**M=0**).  
  
### bpfile:  
informations about loci are required in a file named _bpfile_. Its structure is simple:  
one column per locus  
one line per feature  
line 1: useless line made to contain global informations about the analysis (species names for example, kind of markers, mutation rate, etc ... )  
line 2: locus length **L** in nucleotides.  
line 3: number of sequences in species_A.  
line 4: number of sequences in species_B.  
line 5: populational mutation rate theta = 4.N.mu.**L** with mu is the probability for a nucleotide to be mutated in one generation.  
line 6: populational recombination rate rho. As theta, by replacing **mu** per **r** the per nucleotide recombination rate.  
  

