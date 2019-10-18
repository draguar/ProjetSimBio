# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import numpy as np
import os
import re

def target_expression(environment_file):
    """Get the target expression from given file.
    
    Parameters
    ----------
    environment_file : str
        Path and name of the environment file (giving the target relative
        expression levels)

    Returns
    -------
    target_expresstion : Numpy array
        Array of floats representing target relative expression level for each
        gene.
    """
    with open(environment_file, "r") as environment:
        target_exp = re.findall("[0-9]+\s+([0-9\.]+)", environment.read())
    return np.asarray(target_exp, dtype=float)


def expression_simulation(params_file, out_file):
    """Run  the expression  simulation with given parameters.
    
    Parameters
    ----------
    params_file : str
        Path and name of the parameters file.
    out_file : str
        Path and name of the output file.
        
    Returns
    -------
    transcript_numbers : Numpy array
        Array of ints representing the number of transcripts  for each gene,
        ordered by gene ID.
        
    Note
    ----
    Also generate a text file with the given out_file name.
    """
    
    term_command = ("python3 start_simulation.py " + params_file + " > "
                    + out_file)
    os.system(term_command)
    with open(out_file, "r") as out:
        transcript_nbs = re.findall("Transcript ID [0-9]+ : ([0-9]+)",
                                    out.read())
    return np.asarray(transcript_nbs, dtype=int)


def compute_fitness(observed_transcript_numbers, target_frequencies):
    """Compute the fitness of an individual with given gene expression pattern
    in given environment.
    
    Parameters
    ----------
    observed_transcript_numbers : Numpy array
        Array of ints representing the observed number of transcripts for each
        gene of the invidual.
    target_frequencies : Numpy array
        Array of floats representing target relative expression level for each
        gene (environment).
    
    Returns
    -------
    fitnness : float
        computed fitness of the invidual, following the formula:
            fitness = exp(-sum(|ln(observed_for_gene_i / target_for_gene_i)|))
    """
    
    observed_frequencies = (observed_transcript_numbers
                            / np.sum(observed_transcript_numbers))
    ln_freqs = np.log(observed_frequencies / target_frequencies)
    return np.exp(-np.sum(np.abs(ln_freqs)))

        
def parse_namefile_ini(u):
    """ parse an [INPUTS] line u from a .ini file to get the seeked file location  """
    u = u.split(" ")[2] # get the address of the file of interest
    u = u.rstrip() # removes the "\n" character at the end of the string
    return u


def pos_out_genes(file_ini):
    """" returns a list of open position intervals in the genome where there are not any genes
    with file_ini is the initialization file """
    f = open(file_ini, 'r')
    f_lines = f.readlines() # list of the file_ini lines
    file_tss = parse_namefile_ini(f_lines[2]) # TSS file address
    file_tts = parse_namefile_ini(f_lines[3])# TTS file address
    file_barr = parse_namefile_ini(f_lines[4]) # prot.dat file address
    
    ## list of transcription start position
    ftss = open(file_tss, 'r')
    ftss_lines = ftss.readlines() # list of the file_tss lines
    start = []
    for x in ftss_lines[1:]:
        #print(x.split("\t"))
        start.append(int(x.split("\t")[2]))
    print(start)
    
    ## list of transcription end position
    ftts = open(file_tts, 'r')
    ftts_lines = ftts.readlines() # list of the file_tss lines
    end = []
    for x in ftts_lines[1:]:
        #print(x.split("\t"))
        end.append(int(x.split("\t")[2]))
    print(end)
    
    ## list of topological barriers
    fbarr = open(file_barr, 'r')
    fbarr_lines = fbarr.readlines() # list of the file_tss lines
    barr = []
    for x in fbarr_lines[1:]:
        #print(x.split("\t"))
        barr.append(int(x.split("\t")[1]))
    print(barr)
    
    f.close()
    ftss.close()
    ftts.close()
    fbarr.close()
    
    ### list of open position intervals in the genome where there is not any gene
    n = len(barr)
    out = []
    for i in range(n-1):
        # if we assume that the topological barriers position will always be smaller than the start and end position of its corresponding gene
        interval_a = [barr[i], min(start[i], end[i])] # interval between a topological barrier position and the min(start, end) position of the gene
        out.append(interval_a)
        interval_b = [ max(start[i], end[i]), barr[i+1]] # interval between the max(start, end) position of the gene and the next topological barrier position
        out.append(interval_b)
    out.append([barr[n-1], min(start[n-1], end[n-1])])
    out.append([max(start[n-1], end[n-1]), barr[0]]) # because the genome is circular : interval between the min(start, end) position of the last gene and the first topological barrier
    
    return start, end, barr, out

def sample(out, Ngen):
    """ samples a location from a given list of intervals
    
    Parameters
    ----------
    out : Python list
        listof the open position intervals in the genome where there are not any genes   
    Ngen : int
        the length of the genome
        
    Returns
    -------
    mut_pos : int
        sampled location of the mutation
    
    """
    ### samples the interval where the exact location of the mutation
    ## for each interval, calculates the probability to sample from it
    N = 0 # total length of the intervals
    intlen = [] # list of the cummulative length of each intervals
    for x in out:
        a = x[0] # left bound of the interval x
        b = x[1] # right boud of x
        if a>b:
            print("yes")
            # then x is the last interval of the plasmid with b after the first position of the genome
            xlen = abs(Ngen - a + b - 1)
            #intlen.append(xlen)
            N += xlen
            intlen.append(N)
        else:
            xlen = b-a-1 # length of the interval
            #intlen.append(xlen)
            N += xlen
            intlen.append(N)
            print(xlen)
    probas = np.array(intlen)/N # list of the cummulative probabilities to sample in each interval
    ## samples the interval
    p = np.random.uniform(0,1) # draw a random number between 0 and 1
    sint = min(np.where(p<probas)[0]) # index of the sampled interval
    
    ### samples the exact location of the mutation
    a = out[sint][0]
    b = out[sint][1]
    if a>b:
        # then the selected interval is the last interval of the plasmid with b after the first position of the genome
        mut_pos = np.random.randint(a+1, Ngen+b) # samples the location of the mutation
        if mut_pos > Ngen:
            # then the sampled position is locater after the first position of the plasmid
        mut_pos = mut_pos - Ngen
    else:
        mut_pos = np.random.randint(a+1, b) # samples the location of the mutation
    
    return mut_pos
        


def indel(u):
    """ deletes or inserts in the plasmid a unit of length u of nucleotides """



start, end, barr, out = pos_out_genes("params.ini")
TARGET_FREQS = target_expression("environment.dat")
print (TARGET_FREQS)
initial_expression = expression_simulation("params.ini", "out.txt")
print(initial_expression)
print(compute_fitness(initial_expression, TARGET_FREQS))
