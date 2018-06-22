"""
Faisal Nsour
SSIE 519 Semester Project
Fall 2016

"""
import random

def coinflip(p = 0.5):
    return p > random.random()
    
def crossover_standard(pcrossover, pmutation, mutation_list, parent1, parent2):
    """
    uses standard crossover
    parameters:
        pcrossover - prob. of crossover
        parent1, parent2 - arrays of parent chromosomes
    returns: 2-d array representing newchild1, newchild2
    """
    #rand crossover index
    xidx = random.randint(1, len(parent1) - 1)
    child1, child2 = parent1, parent2
    if coinflip(pcrossover):
        child1 = mutate_swap(pmutation, parent1[:xidx] + parent2[xidx:])
        child2 = mutate_swap(pmutation, parent2[:xidx] + parent1[xidx:])
    return child1, child2

def crossover_by_cluster(pcrossover, pmutation, parent1, parent2):
    """
    Performs crossover specific to locus-based adjacency chromosomes as described
     in Shi et al.
    """        
    p1_idx = random.randint(0, len(parent1) - 1)
    p2_idx = random.randint(0, len(parent2) - 1)
    child1, child2 = parent1, parent2
    if coinflip(pcrossover):
        p1_cluster = collect_lba_cluster(parent1, p1_idx)
        p2_cluster = collect_lba_cluster(parent2, p2_idx)
        child1 = mutate_lba(pmutation, overwrite_lba_cluster(parent1, p2_cluster))
        child2 = mutate_lba(pmutation, overwrite_lba_cluster(parent2, p1_cluster))
    return child1, child2

def overwrite_lba_cluster(lba, cluster):
    """
    LBA = locus-based adjacency list
    takes an LBA and overwrites positions identified in a cluster from another LBA
    parameters:
        -lba
        -cluster: list of edges; eg, [[1,2], [2, 1]]
    """
    _lba = lba[:]
    for e in cluster:
        _lba[e[0]] = e[1]
    return _lba

def collect_lba_cluster(lba, start_idx):
    """
          #go to px_idx; read values until no more indexes to follow; copy into other chromo
        # stop when the value is eq to idx or value is eq to idx already visited

    """
    idx = start_idx
    cluster = [] #list of edges (idx, node)
    while True:
        cluster.append([idx, lba[idx]])
        if idx == lba[idx] or lba[idx] in [e[0] for e in cluster]:
            break
        else:
            idx = lba[idx] #update index    
    return cluster
        
def mutate_lba(p, chromosome):
    """
    """
    for gene, allele in enumerate(chromosome):
        if coinflip(p):
            swap_gene = random.randrange(0, len(chromosome))
            chromosome[gene] = swap_gene
    return chromosome
    
def mutate_swap(p, chromosome):
    """
    Mutates by swapping a gene's allele with another randomly selected gene
    
    parameters:
        p - prob. of mutation
        chromosome - individual to mutate
    
    """
    for idx, allele in enumerate(chromosome):
        if coinflip(p): #mutate
            swap_idx = random.randrange(0, len(chromosome))
            temp = chromosome[swap_idx]
            chromosome[swap_idx] = chromosome[idx] 
            chromosome[idx] = temp
    return chromosome

def mutate_from_list(p, mutation_list, chromosome):
    """
        parameters:
            p - prob. of mutation
            mutation_list - values to choose from when mutating
            chromosome - individual to mutate
    """
    for idx, allele in enumerate(chromosome):
        if coinflip(p): #mutate
            chromosome[idx] = random.choice(mutation_list)
    return chromosome
    
def select(raw_pop_fit, breeding_factor):
    """
    Selects a index from a list of fitness values using roullete wheel selection
    parameters:
        pop_fit - list of fitness values
        
    return: index of selected individual
    """

    avg_raw = sum(raw_pop_fit) / float(len(raw_pop_fit))
    best_raw = max(raw_pop_fit)
    worst_raw = min(raw_pop_fit)
    scaled_pop_fitness = [scaled_fitness(i, avg_raw, best_raw, worst_raw, breeding_factor) for i in raw_pop_fit]

    sum_fit = sum(scaled_pop_fitness)
    wheel_point = random.random() * sum_fit
    #print "wheel_point", sum(pop_fit), wheel_point
    i = 0
    partial_sum = scaled_pop_fitness[0]
    #print "part sum", partial_sum
    while partial_sum < wheel_point and i < len(scaled_pop_fitness)-1:
        i += 1
        partial_sum += scaled_pop_fitness[i]
        
    #print "Selected", i
    return i
        
def scaled_fitness(fitness, avg, best, worst, breeding_factor):
    a = b = 0.0
    if (avg - worst) <= ((best - avg) / (breeding_factor - 1)):
        a = (breeding_factor - 1) / (best - avg)
        b = (best - breeding_factor*avg) / (best - avg)
    else:
        a = 1 / (avg - worst)
        b = -worst / (avg - worst)
    return (fitness * a) + b
