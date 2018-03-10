#!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

import sys
from numpy.random import uniform
from numpy.random import binomial
from random import shuffle
help = "\t\033[1;31;40mTakes one model specifier and a number of multilocus simulations as argument:\033[0m\n\t\t"
help += "\n\t\t".join(["SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N", "AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N", "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N", "SI_1N", "SI_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#SI\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#AM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;30m#PAM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 0 -m 2 1 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#IM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#SC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#PSC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ej tbs 2 1 -eN tbs tbs\n\n"
help += "\t\033[1;32;40mExample: ./priorgen.py SC_2M_2N 1000\033[0m\n"

if len(sys.argv) != 3:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])
N_bound = [0, 50]
T_bound = [0, 50]
M_bound = [0, 40]
P_ntrl_M_bound = [0, 1] # ntrl_M loci with M >0; (1-ntrl_M) loci with M=0
P_ntrl_N_bound = [0.5, 1] # ntrl_N loci with N and (1-ntrl_N) with bf.N
bf_bound = [0, 1] # reduction in Ne for the (1-ntrl_N)

# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = infile.readline().strip().replace(" ", "\t").split("\t")
nsamA = infile.readline().strip().replace(" ", "\t").split("\t")
nsamB = infile.readline().strip().replace(" ", "\t").split("\t")
theta = infile.readline().strip().replace(" ", "\t").split("\t")
rho = infile.readline().strip().replace(" ", "\t").split("\t")
infile.close()

# number of loci
nLoci = len(L)

# sum of nsamA + nsamB
nsam_tot = [ int(nsamA[i]) + int(nsamB[i]) for i in range(nLoci) ]



if sys.argv[1] == "SC_1M_1N":
	# secondary contact
	# ms tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
	# nsamtot theta rho L nsamA nsamB M12 M21 N1 N2 Tsc Tsplit Tsplit Na

	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#	Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1[sim], N2[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#       Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "SC_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#       Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## number of neutral loci
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7:.5f}\t{8}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		
		# vectors of size 'nLoci' containing parameters
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1[sim], N2[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tTsc\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8}\t{9:.5f}\t{10}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], Tsc[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
        #Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## number of neutral loci
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7:.5f}\t{8}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		
		# vectors of size 'nLoci' containing parameters
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tTam\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8}\t{9:.5f}\t{10}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], Tam[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5}\t{6:.5f}\t{7}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		
		# vectors of size 'nLoci' containing parameters
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7}\t{8:.5f}\t{9}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "PSC_2M_2N":
	# Tpresent --> Tmig --> Tmig + Tiso --> 2.Tmig + Tiso --> Tsplit
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tmig = [ uniform(low = 0, high = Tsplit[i], size = 1)[0]/2.0 for i in range(nMultilocus) ] # 
	Tiso = [ Tsplit[i]/2.0-Tmig[i] for i in range(nMultilocus) ]
#	print(" ".join([ str(i) for i in Tsplit ]))
#	print(" ".join([ str(i) for i in Tmig ]))
#	print(" ".join([ str(i) for i in Tiso ]))
#	print()

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tTmig\tTiso\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9}\t{10:.5f}\t{11}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], Tmig[sim], Tiso[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
			# msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ej tbs 2 1 -eN tbs tbs
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], M12_vec[locus], M21_vec[locus], Tmig[sim], Tiso[sim]+Tmig[sim], M12_vec[locus], M21_vec[locus], 2*Tmig[sim]+Tiso[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "PAM_2M_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 0 -m 2 1 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs 
	# Tpresent --> Tiso --> Tmig + Tiso --> 2.Tiso + Tmig --> Tsplit
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tmig = [ uniform(low = 0, high = Tsplit[i], size = 1)[0]/2.0 for i in range(nMultilocus) ] # 
	Tiso = [ Tsplit[i]/2.0-Tmig[i] for i in range(nMultilocus) ]
	

	## bf = factor of local reduction in Ne. Model of "background selection"
	bf = uniform(low=bf_bound[0], high=bf_bound[1], size = nMultilocus)

	## number of neutral loci
	ntrl_N = [ binomial(nLoci, uniform(P_ntrl_N_bound[0], P_ntrl_N_bound[1])) for i in range(nMultilocus) ]
	ntrl_M12 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]
	ntrl_M21 = [ binomial(nLoci, uniform(P_ntrl_M_bound[0], P_ntrl_M_bound[1])) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tprop_ntrl_N\tbf\tTsplit\tTmig\tTiso\tM12\tprop_ntrl_M12\tM21\tprop_ntrl_M21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9}\t{10:.5f}\t{11}\n".format(N1[sim], N2[sim], Na[sim], ntrl_N[sim], bf[sim], Tsplit[sim], Tmig[sim], Tiso[sim], M12[sim], ntrl_M12[sim], M21[sim], ntrl_M21[sim])
		# vectors of size 'nLoci' containing parameters
		bf_vec = ntrl_N[sim]*[1] + (nLoci - ntrl_N[sim])*[bf[sim]]
		shuffle(bf_vec)
		N1_vec = [ N1[sim]*bf_vec[i] for i in range(nLoci) ]
		N2_vec = [ N2[sim]*bf_vec[i] for i in range(nLoci) ]
		Na_vec = [ Na[sim]*bf_vec[i] for i in range(nLoci) ]
		
		M12_vec = [M12[sim]]*ntrl_M12[sim] + [0]*(nLoci - ntrl_M12[sim])
		shuffle(M12_vec)
		M21_vec = [M21[sim]]*ntrl_M21[sim] + [0]*(nLoci - ntrl_M21[sim])
		shuffle(M21_vec)
		
		for locus in range(nLoci):
		# msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 0 -m 2 1 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs 
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tiso[sim], M12_vec[locus], M21_vec[locus], Tiso[sim]+Tmig[sim], 2*Tiso[sim]+Tmig[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

