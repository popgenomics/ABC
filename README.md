# simple pipeline to perform simulations for downstream ABC analysis  
## requirement:  
### python:  
Be sure that you have **python** and *pypy* in your /usr/bin directory  
The pipeline also works with only the standard _python_ interpreter, but is much slower.  
Comparison using the _time_ command for the simulation of _1,000_ multilocus datsets of _100_ loci:  
_pypy_ : real = 1m34.698s   
_python_ : real = 21m24.772s   
  
If you do not have _pypy_ installed and still want to use this pipeline despite the extra computation time, you just have to change the first line of _mscalc.py_ by replacing _#!/usr/ bin/pypy_ by _#!/usr/bin/python_  
  
Finaly, _priorgen.py_ imports the _numpy_ library.  
  
### models: 
Four categories of models are allowed:  
SI = _S_trict _I_solation model (split of an ancestral population in two isolated daughter populations)  
AM = _A_ncient _M_igration model. Like SI, but with migration between daughter populations between the time of spli _Tsplit_ and _Tam_. There is then no migration between _Tam_ and the present time.  
SC = _S_econdary _C_ontact model. Like SI, but with no migration between daughter populations between the time of split _Tsplit_ and _Tsc_. There is then a secondary contact with continuous migration between _Tsc_ and the present time.  
IM = _I_solation with _M_igration model. Like SI but with gene flow between _Tsplit_ and the present time.  

### submodels:  
Two models of effective population size  (_Ne_):  
1N = all loci share the same _Ne_  
2N = two categories of loci with a neutral _Ne_ (proportion _P_) and selectively-linked loci (_1-P_) with a size equal to _bf.Ne_ where _bf_ is a background factor reducing the _Ne_ of a locus_i similarly in the 3 populations (ancestral, pop_A and pop_B).  
  
1M = all loci share the same migration rate _M_.  
2M = two categories of loci with a neutral _M_ (>0) or loci linked to species barriers (_M=0_).  
  
### bpfile:  
informations about loci are required in a file named _bpfile_. Its structure is simple:  
one column per locus  
one line per feature  
line 1: useless line made to contain global informations about the analysis (species names for example, kind of markers, mutation rate, etc ... )  
line 2: locus length _L_ in nucleotides.  
line 3: number of sequences in species_A.  
line 4: number of sequences in species_B.  
line 5: populational mutation rate theta = 4.N.mu._L_ with mu is the probability for a nucleotide to be mutated in one generation.  
line 6: populational recombination rate rho. As theta, by replacing _mu_ per _r_ the per nucleotide recombination rate.  
  

