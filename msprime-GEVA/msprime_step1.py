import msprime
import tskit
import numpy as np
import math
import pandas as pd
import os
import itertools
import random

####Assuming msprime version 1 
myMap="/home/mbmenon/AgeEst/RobertRuns/Slim_map.txt"
mu=1e-6
N=500 #this is diploid num of samples



##**********MODEL 1: NO pop change***************###


demography = msprime.Demography()
demography.add_population(name="slim",initial_size=500)

demography.add_population_parameters_change(time=449, population="slim", initial_size = 500);

recomb_map =  msprime.RateMap.read_hapmap(myMap, has_header=True)
ts = msprime.sim_ancestry(samples={"slim": N}, recombination_rate=recomb_map,
                          demography=demography, random_seed=22)

muts = msprime.sim_mutations(ts, rate=mu)



with open("SlimBurst_WithMap.vcf", "w") as vcf_file:
     muts.write_vcf(vcf_file)

#Get genotype matrix per haploid genome
A = muts.genotype_matrix()
VAR=muts.variants()
POS=[]
for variant in muts.variants():
    POS.append(variant.site.position)

df01 = pd.DataFrame(A)
df01['pos'] = POS
df01.head()

df01.to_csv(path_or_buf="SlimBurst_genotypesPerHaplotype1000Haps.txt", sep='\t', na_rep='NA',
          float_format=None, columns=None, header=True, index=False,
          mode='w', encoding=None, compression='infer', quoting=None, quotechar='"',
          line_terminator=None, chunksize=None, date_format=None, doublequote=True,
          escapechar=None, decimal='.', errors='strict', storage_options=None)


###Get TMRCA and segment positions per pair

comb = list(itertools.combinations(range(0,2*N), 2))

dfFull = pd.DataFrame()

def extractPairs(PW):
    df = pd.DataFrame()
    ts_subset = muts.simplify(PW)
    TMRCA=[]
    Beg=[]
    End=[]
    for tree in ts_subset.trees():
        Beg.append(tree.interval[0])
        End.append(tree.interval[1])
        TMRCA.append(tree.time(tree.roots[0]))
    df['TMRCA'] = TMRCA
    df['beg'] = Beg
    df['end'] = End
    df['pair'] = list(itertools.repeat(PW, len(TMRCA)))
    return df


FIRSTset = [extractPairs(i) for i in comb]
dfFull = pd.DataFrame()
cobimbe = dfFull.append(FIRSTset)


cobimbe.to_csv(path_or_buf="SlimBurstAllpairs_1000haps.csv", sep=',', na_rep='NA',
          float_format=None, columns=None, header=True, index=False,
          mode='w', encoding=None, compression='infer', quoting=None, quotechar='"',
          line_terminator=None, chunksize=None, date_format=None, doublequote=True,
          escapechar=None, decimal='.', errors='strict', storage_options=None)



##**********MODEL 2: Pop bottelneck***************###


demography = msprime.Demography()
demography.add_population(name="slim", initial_size=500)

demography.add_population_parameters_change(time=4751, population="slim", initial_size = 500);
demography.add_population_parameters_change(time=4999, population="slim", initial_size = 50);
demography.add_population_parameters_change(time=5000, population="slim", initial_size = 500);

recomb_map =  msprime.RateMap.read_hapmap("/home/mbmenon/AgeEst/RobertRuns/Slim_map.txt", has_header=True)
ts = msprime.sim_ancestry(samples={"slim": N}, recombination_rate=recomb_map,
                          demography=demography, random_seed=22)

muts = msprime.sim_mutations(ts, rate=mu)


####
with open("SlimBott_WithMap.vcf", "w") as vcf_file:
     muts.write_vcf(vcf_file)

#Get genotype matrix
A = muts.genotype_matrix()
VAR=muts.variants()
POS=[]
for variant in muts.variants():
    POS.append(variant.site.position)

df01 = pd.DataFrame(A)
df01['pos'] = POS
df01.head()

df01.to_csv(path_or_buf="SlimBott_genotypesPerHaplotype1000Haps.txt", sep='\t', na_rep='NA',
          float_format=None, columns=None, header=True, index=False,
          mode='w', encoding=None, compression='infer', quoting=None, quotechar='"',
          line_terminator=None, chunksize=None, date_format=None, doublequote=True,
          escapechar=None, decimal='.', errors='strict', storage_options=None)



###Get TMRCA and segment positions per pair using function extractPair defined above

comb = list(itertools.combinations(range(0,2*N), 2))

dfFull = pd.DataFrame()

FIRSTset = [extractPairs(i) for i in comb]

dfFull = pd.DataFrame()
cobimbe = dfFull.append(FIRSTset)


cobimbe.to_csv(path_or_buf="SlimBottAllpairs_1000haps.csv", sep=',', na_rep='NA',
          float_format=None, columns=None, header=True, index=False,
          mode='w', encoding=None, compression='infer', quoting=None, quotechar='"',
          line_terminator=None, chunksize=None, date_format=None, doublequote=True,
          escapechar=None, decimal='.', errors='strict', storage_options=None)