import numpy as np
import pandas as pd
import dask.dataframe as dd
import multiprocessing
from dask.multiprocessing import get
from multiprocessing import cpu_count
import numpy as np
import os

#Input files from msprime_step1
#Run seperately for each model
model="Burst"
TMRCAsegs="Slim"+model+"Allpairs_1000haps.csv"
genoMatrix="Slim"+model+"_genotypesPerHaplotype1000Haps.txt"
IBDout="ibd_1000haps_"+model+".txt"

##Custom functions

def IsBetween(df):
    #here assuming the file with allelefrequency is called afD_df with two cols (af & pos)
    val = afD_df['pos'].between(df['beg'],df['end'])
    return afD_df[val]


def combineSets(part):
    #parted_df is the dask partitioned dataframe containing the TMRCAsegs
    partDF = parted_df.partitions[part]
    FIRSTset = [IsBetween(i) for index, i in  partDF.iterrows()]
    FIRSTset_filter = [x for x in FIRSTset if x.shape[0].compute() > 10]
    setIBD = dd.concat(FIRSTset_filter,axis=0)
    return setIBD


###################

dask_df = dd.read_csv(TMRCAsegs)
condition = dask_df['TMRCA']>0 
dask_df = dask_df[condition]
dask_df=dask_df.drop_duplicates(subset=['TMRCA','beg','end'])

dask_df['logTMRCA'] = np.log2(dask_df['TMRCA'])


#deciding on the cutoff for IBD (either median or 100TMRCA)
#med=dask_df.quantile(0.5).compute() 
condition = dask_df['logTMRCA']<=6.54
IBD = dask_df[condition]

del dask_df


geno = dd.read_table(genoMatrix)
geno_1 = geno.drop('pos', axis=1)
geno_1=geno_1[(geno_1 != 2.0).all(axis=1)] #remove all multiallelic positions 

afD = geno_1.mean(axis=1)
afD_df = pd.DataFrame(afD,columns=['af'])
afD_df['pos'] = geno['pos']
afD_df['pos'] = afD_df['pos']

del geno
del geno_1


#convert adD_df also to DASK
afD_df= dd.from_pandas(afD_df, npartitions=2*multiprocessing.cpu_count())

###Trying to do partitions of IBD and then implement the defination for each partition
#npart = round(len(IBD)/1000)
npart=500 #this can be optimised depending on size of the dataset
parted_df = IBD.repartition(npartitions=npart)



#Now for each unique segment determine the af of all the variants that fall within it and keep if the #positions > 10
IBDsets = [combineSets(i) for i in  range(npart)]

IBDsets_comb = dd.concat(IBDsets,axis=0)
IBDsets_comb.to_csv(IBDout,sep="\t",single_file=True)