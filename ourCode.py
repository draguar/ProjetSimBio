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

def genome_inversion(genome_size, genes_start_pos, genes_end_pos, barriers_pos,
                     inversion_start, inversion_end):        
    """Perform a genome innversion onn given genome at given positions..
    
    Parameters
    ----------
    genome_size : int
        Genome size in base pair.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    inversion_start : int
        Position of the beginning of the inversion in the genome.
    inversion_end : int
        Position of the end of the inversion in the genome.
    
    Returns
    -------
    genes_start_pos : Numpy array
        Updated value of genes_start_pos after inversion.
    genes_end_pos : Numpy array
        Updated value of genes_end_pos after inversion.
    barriers_pos : Numpy array
        Updated value of barriers_pos after inversion.
        
    Notes
    -----
    inversion_start and inversion_end must not be inside a gene.
    inversion_end must be greater than inversion_start.
    """
    
    new_genes_start_pos = np.copy(genes_start_pos)
    new_genes_end_pos = np.copy(genes_end_pos)
    new_barriers_pos = np.copy(barriers_pos)
    for new_array in [new_genes_start_pos, new_genes_end_pos,
                      new_barriers_pos]:
        for (k, position) in enumerate(new_array):
            if (position > inversion_start) and (position < inversion_end):
                new_array[k] = inversion_start + (inversion_end - position)
    return (new_genes_start_pos, new_genes_end_pos, new_barriers_pos)
                
    
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


start, end, barr, out = pos_out_genes("params.ini")
TARGET_FREQS = target_expression("environment.dat")
print (TARGET_FREQS)
initial_expression = expression_simulation("params.ini", "out.txt")
print(initial_expression)
print(compute_fitness(initial_expression, TARGET_FREQS))
print(genome_inversion(30000, start, end, barr, 2020, 13000))