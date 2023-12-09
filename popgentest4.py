import numpy
import matplotlib
from collections import defaultdict
from matplotlib.pyplot import plot,figure,subplot,xlim,ylim,\
  subplots_adjust,legend,xlabel,ylabel
import random as random
from numpy import linspace
from matplotlib.cm import viridis,hsv,Spectral 
from matplotlib.pyplot import style
from collections import defaultdict


#Make plots have a dark background
style.use('dark_background')

s=0.1
h=0.7
pop_sizes = [1000]
generations = 1000
replicates = 10
Wpp = 0.1
#Wpq = (1-h*s)*Wpp
Wqq = (1-s)*Wpp

p = 0.3
q = 1 - p



def random_genotype(f_A1):
    if random.uniform(0,1) <= f_A1:
        sperm_allele = 'A1'
    else:
        sperm_allele = 'A2'
    
    if random.uniform(0,1) <= f_A1:
        egg_allele = 'A1'
    else:
        egg_allele = 'A2'
    
    genotype = sperm_allele+egg_allele
    return genotype


def simulate_random_mating(pop_sizes,f_A1,hi):
    genotypes = defaultdict(int)
    global s
    
    for i in range(pop_sizes):
        curr_genotype = random_genotype(f_A1)
        genotypes[curr_genotype]+=1
    ghet=defaultdict(int)
    hval=[]
    for i in range(genotypes["A1A2"]+genotypes["A2A1"]):
       t=float(numpy.random.normal(hi, 0.05, 1))
       if random.uniform(0, 1) <= (1-s*t)*Wpp :
            ghet[1]+=1
            hval.append(t)
       if len(hval)>=1:
           hi=sum(hval)/len(hval)
       if hi>1:
           hi=1
       elif hi<0:
           hi=0
    y=genotypes["A1A1"]
    c=genotypes["A2A2"]
    genotypes["A2A1"]=(0.5*ghet[1])/(genotypes["A1A1"]+ghet[1]+genotypes["A2A2"])
    genotypes["A1A2"]=(0.5*ghet[1])/(genotypes["A1A1"]+ghet[1]+genotypes["A2A2"])
    genotypes["A1A1"]=(genotypes["A1A1"])/(y+ghet[1]+c)
    genotypes["A2A2"]=(genotypes["A2A2"])/(y+ghet[1]+c)
    #print(hi)
    
    
    return genotypes,ghet,hi


def simulate_genetic_drift_with_selection(n_generations,f_A1,pop_size,Wpp,Wqq):
    """Return the frequency of f_A1 over generations"""
    generations = range(n_generations)
    allele_freqs = []
    global h
    hi=h
    for generation in generations:
        genotypes = None
        #print(pop_size)
        genotypes,ghet,hi = simulate_random_mating(pop_size,f_A1,hi)
        #print(genotypes)
        numerator = (genotypes['A1A1']*Wpp+0.5*genotypes['A1A2']+0.5*genotypes['A2A1'])
        denominator = (genotypes['A1A1']*Wpp+genotypes['A1A2']+genotypes['A2A1']+genotypes['A2A2']*Wqq)
        #print("Numerator:",numerator)
        #print("Denominator:",denominator)
        f_A1 = numerator/denominator
        allele_freqs.append(f_A1)
        #print(allele_freqs)
    return list(generations),allele_freqs,hi

#Set up a figure to graph the data
figure(1,dpi=800)
toth=0.0
for i,pop_size in enumerate(pop_sizes):
    print("Simulating pop_size:%i"%pop_size)
    #Make a subplot with one column and one row per population size. 
    nrows = len(pop_sizes)
    ncols = 1
    subplot(nrows,ncols,i+1)
    ylabel("f_A1\n n = %i" %pop_size)
    for replicate in range(replicates):
        xs,ys,H= simulate_genetic_drift_with_selection(generations,p,pop_size,Wpp,Wqq)
        plot(xs,ys,'-',label = "f(A1) pop = %i"%pop_size,linewidth=0.5)
        #Set x-axis limits
        xlim(0,generations)
        #Set y-axis limits
        ylim (0.0,1.0)
        #print(h)
        toth=toth+H
print("FINAL IS:",toth/replicates)
#legend()
xlabel("Generations")
#Leave more space between subplots so they don't overlap
subplots_adjust(hspace=1.00)
#matplotlib.pyplot.show()
