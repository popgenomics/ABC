# simple pipeline to perform simulations for downstream ABC analysis  
## output:  
The pipeline produces two tables:  
**ABCstat.txt** : one line per multilocus simulation, one column per computed summary statistic.  
**priorfile.txt**: one line per multilocus simulation, one column per model statistic  

## requirement:  
### coalescent simulator:  
The used coalescent simulator is **msnsam** made by Jeffrey Ross-Ibarra from original **ms** coded by Richard Hudson.  
In the subdirectory msnsam, just type **./clms** to install the simulator, and put it in your *bin/*  
  
### python:  
Be sure that you have **python** and **pypy** in your /usr/bin directory  
The pipeline also works with only the standard **python** interpreter, but is much slower.  
Comparison using the **time** command for the simulation of **1,000** multilocus datsets of **100** loci:  
**pypy** : real = 1m34.698s   
**python** : real = 21m24.772s   
  
If you do not have **pypy** installed and still want to use this pipeline despite the extra computation time, you just have to change the first line of **mscalc.py** by replacing **#!/usr/ bin/pypy** by **#!/usr/bin/python**  
  
You can link **mscalc.py** and **priorgen.py** to your /usr/bin/ directory.  

Finaly, **priorgen.py** imports the **numpy** library.  
  
## models: 
Four categories of models are allowed:  
SI = **S**trict **I**solation model (split of an ancestral population in two isolated daughter populations)  
AM = **A**ncient **M**igration model. Like SI, but with migration between daughter populations between the time of spli **Tsplit** and **Tam**. There is then no migration between **Tam** and the present time.  
SC = **S**econdary **C**ontact model. Like SI, but with no migration between daughter populations between the time of split **Tsplit** and **Tsc**. There is then a secondary contact with continuous migration between **Tsc** and the present time.  
IM = **I**solation with **M**igration model. Like SI but with gene flow between **Tsplit** and the present time.  
  
Each model SI, AM, SC and IM has its own *msnsam* command line. You can get them by simply typing *./priorgen.py**  
        **#SI**  
        msnsam tbs ${nsim} -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs  
        **#AM**  
        msnsam tbs ${nsim} -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs  
        **#IM**  
        msnsam tbs ${nsim} -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs  
        **#SC**  
        msnsam tbs ${nsim} -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs  
Where ${nsim} is an integer value corresponding to : the number of loci (100 in the example) and the number of targeted simulated multilocus datasets (1,000 in the example).  

## submodels:  
Two models of effective population size  (**Ne**):  
1N = all loci share the same **Ne**  
2N = two categories of loci with a neutral **Ne** (proportion **P**) and selectively-linked loci (**1-P**) with a size equal to **bf.Ne** where **bf** is a background factor reducing the **Ne** of a locus_i similarly in the 3 populations (ancestral, pop_A and pop_B).  
  
1M = all loci share the same migration rate **M**.  
2M = two categories of loci with a neutral **M** (>0) or loci linked to species barriers (**M=0**).  
  
## bpfile:  
informations about loci are required in a file named _bpfile_. Its structure is simple:  
one column per locus  
one line per feature  
line 1: useless line made to contain global informations about the analysis (species names for example, kind of markers, mutation rate, etc ... )  
line 2: locus length **L** in nucleotides.  
line 3: number of sequences in species_A.  
line 4: number of sequences in species_B.  
line 5: populational mutation rate theta = 4.N.mu.**L** with mu is the probability for a nucleotide to be mutated in one generation.  
line 6: populational recombination rate rho. As theta, by replacing **mu** per **r** the per nucleotide recombination rate.  
  
  
## Example:  
100 loci.  
1 million of multilocus simulations.  
model: SI with 2Ne  
**priorgen.py** and **mscalc.py** are assumed to be linked in your /usr/bin  
Because 1 million of simulations can take a while, and thanks to a cluster of CPUs, 100 set of 10,000 simulations can run in parallele as following:  
for i in $(seq 1 1 100); do  
	echo "mkdir SI_2N_${i}; cd SI_2N_${i}; cp ../bpfile .; priorgen.py SI_2N 10000 | msnsam tbs 1000000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs | mscalc.py"  
done  


	
