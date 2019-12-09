# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import pandas as pd

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


def expression_simulation(params_file, out_file, gene_start_pos):
    """Run  the expression  simulation with given parameters.
    
    Parameters
    ----------
    params_file : str
        Path and name of the parameters file.
    out_file : str
        Path and name of the output file.
    gene_start_pos : Numpy array
        Array of ints representing the begining position of genes.
        
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
        file_content = out.read()
        transcript_nbs = re.findall("Transcript ID [0-9]+ : ([0-9]+)",
                                    file_content)
        transcript_nbs = np.asarray(transcript_nbs, dtype=int)
        # Transcripts are reindexed, we fetch their starting position to map
        # transcription values to the right gene.
        starts = re.findall("\n[0-9] +[0-9] +[0-9] +[\-0-9]+ +([0-9]+)",
                            file_content)
        starts = np.asarray(starts, dtype=int)
    # Reorder transcription values
    transcription_values = np.full(len(starts), np.nan)
    for k, start_pos in enumerate(gene_start_pos):
        transcript_id, = np.where(starts == start_pos)
        transcription_values[k] = transcript_nbs[transcript_id[0]]
    return transcription_values

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

def accept_mutation(previous_fitness, new_fitness, q):
    """Accept or reject the mutation, based on fitnesses comparison.
    
    A fitness increase is always accepted. A fitness loss is accepted with
    a probability exp(fitness difference / q).
    
    Parameters
    ----------
    previous_fitness : float
        Fitness of the previous generation (before mutation).
    new_fitness : float
        Fitness of the new generation (after mutation).
    q : float
        Parameters of the Monte Carlo Metropolis algorithm, controlling the
        range of accepted fitness losses.
        
    Returns
    -------
    is_accepted : bool
        True if the mutation is accepted, False else.
    """
    
    fitness_diff = new_fitness - previous_fitness
    if fitness_diff > 0:
        return True
    else:
        return (np.random.rand() < np.exp(fitness_diff/q))


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
    event_type : str
        Always equal to "inversion"
    genome_size : int
        Unchanged value of genome_size.
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
    return ("inversion", genome_size, new_genes_start_pos, new_genes_end_pos,
            new_barriers_pos)
                
    
def parse_namefile_ini(u):
    """ parse an [INPUTS] line u from a .ini file to get the seeked file location  """
    u = u.split(" ")[2] # get the address of the file of interest
    u = u.rstrip() # removes the "\n" character at the end of the string
    return u


def pos_out_from_pos_lists(genes_start_pos, genes_end_pos, barriers_pos):        
    """Generates the list of postition intervals outside barriers and genes.
    
    Parameters
    ----------
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    
    Returns
    -------
    out_positions : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier.range
    """
    
    limits = np.sort(np.hstack((genes_start_pos, genes_end_pos, barriers_pos,
                                barriers_pos)))
    out_positions = np.array([limits[-1], limits[0]])
    for k in range(1, len(limits)-1, 2):
        out_positions = np.vstack((out_positions, [limits[k], limits[k+1]]))
    return out_positions


def pos_out_genes(file_ini):
    """" returns a list of open position intervals in the genome where there are not any genes
    with file_ini is the initialization file """
    f = open(file_ini, 'r')
    f_lines = f.readlines() # list of the file_ini lines
    file_gff = parse_namefile_ini(f_lines[1]) # GFFF file address
    file_tss = parse_namefile_ini(f_lines[2]) # TSS file address
    file_tts = parse_namefile_ini(f_lines[3])# TTS file address
    file_barr = parse_namefile_ini(f_lines[4]) # prot.dat file address
    
    ## genome size
    fgff = open(file_gff, 'r')
    fgff_content = fgff.read() # list of the file_gff lines
    Ngen = re.findall("##sequence-region .* 1 ([0-9]+)", fgff_content)[0]
    print(Ngen)
                      
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
    fgff.close()
    ftss.close()
    ftts.close()
    fbarr.close()
    
    ### list of open position intervals in the genome where there is not any gene
    out = pos_out_from_pos_lists(start, end, barr)
    return start, end, barr, out, int(Ngen)

def sample(out, Ngen, u):
    """ samples a location from a given list of intervals that satisfies the "distance to bounds condition" 
    (there is at least u nucleotides between the right and left bound and the sampled mutation position)
    
    Parameters
    ----------
    out : Python list
        list of the open position intervals in the genome where there are not any genes   
    Ngen : int
        the length of the genome
    u : int
        unit of length of nucleotides.
        
    Returns
    -------
    mut_pos : int
        sampled location of the mutation
    
    """
    space = False # boolean used to make sure there is enough space between the right and left bound and the sampled mutation position
    i=100 # max number of iterations
    
    while (not space):
        # we repeat until we find a mutation position that satisfies the distance to the bounds condition
        # or until we have repeated this loop a 100 times
        #print(i)
        ### samples the interval where the exact location of the mutation
        ## for each interval, calculates the probability to sample from it
        N = 0 # total length of the intervals
        intlen = [] # list of the cummulative length of each intervals
        for x in out:
            a = x[0] # left bound of the interval xout_positions[i][1]
            b = x[1] # right boud of x
            if a>b:
                # then x is the last interval of the plasmid with b after the first position of the genome
                xlen = abs(Ngen - a + b - 1)
                N += xlen
                intlen.append(N)
            else:
                xlen = b-a-1 # length of the interval
                N += xlen
                intlen.append(N)
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
                # then the sampled position is located after the first position of the plasmid
                
                if (mut_pos <= Ngen + out[sint][1] - u) and (mut_pos >= out[sint][0] + u):
                    # there is enough space between the right and left bound and the sampled mutation position
                    space = True
                    mut_pos = mut_pos - Ngen
                   
        else:
            mut_pos = np.random.randint(a+1, b) # samples the location of the mutation
            if (mut_pos <= out[sint][1] - u) and (mut_pos >= out[sint][0] + u):
                # there is enough space between the right and left bound and the sampled mutation position
                space = True
            
            
        #print(i, " ", mut_pos)
        
        i-=1
        if i==0:
            raise RuntimeError("A mutation position that that satisfies the distance to the bounds condition was not found")
    
    return mut_pos
        


def indel(u, genome_size, genes_start_pos, genes_end_pos, barriers_pos):
    """ deletes or inserts in the plasmid a unit length u of nucleotides 
    
    Parameters
    ----------
    u : int
        unit of length of nucleotides that is deleted or inserted.
    genome_size : int
        Genome size in base pair.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
        
    Returns
    -------
    event_type : str
        Equal to "insertion" or "deletion".
    new_genes_start_pos : Numpy array
        Updated value of genes_start_pos after inversion.
    new_genes_end_pos : Numpy array
        Updated value of genes_end_pos after inversion.
    new_barriers_pos : Numpy array
        Updated value of barriers_pos after inversion.
    genome_size: int
        Updated value of the genome size.
    
    """
    ### sample the indel position
    indel_pos = 27613 #sample(out_positions, genome_size)
    #print(indel_pos)
    
    ### initialization of the new positions
    new_genes_start_pos = np.copy(genes_start_pos)
    new_genes_end_pos = np.copy(genes_end_pos)
    new_barriers_pos = np.copy(barriers_pos)
    
    ### choose whether it is an insertion or a deletion
    p = np.random.uniform(0,1) # draw a random number between 0 and 1
    if p<.5:
        # it is an insertion
        #print("insert")
        event_type = "insertion"
        for i in range( len(genes_start_pos) ):
            if genes_start_pos[i] >= indel_pos :
                # update gene start positions
                new_genes_start_pos[i] += u
            if genes_end_pos[i] >= indel_pos :
                # update gene end positions
                new_genes_end_pos[i] += u
            if barriers_pos[i] >= indel_pos :
                # update barrieres positions
                new_barriers_pos[i] += u
        # update genome size
        genome_size += u
        
    else:
        # it is a deletion
        event_type = "deletion"
        #print("delete")
        
        ### check instead in the out list, create a dict of right interval bound : right interval bound - u_new (don't forget to check there is enough space to delete) then update start, end , barr accordingly
        
        ## find the intervals where the sampled mutation position is
#        right_bound = 0 # right bound of the interval where the mutation position is
#        i=0
#        for i in range(len(out_positions)):
#            print(i)
#            if out_positions[i][0] <= out_positions[i][1]:
#                # then this is not the last interval
#                if out_positions[i][0] <= indel_pos <= out_positions[i][1]:
#                    right_bound = out_positions[i][1]-1
#            else:
#                # this is the last interval and out_positions[i][1] is after the first position in the genome
#                if (out_positions[i][0] <= indel_pos <= genome_size) or (1 <= indel_pos <= out_positions[i][1]):
#                    right_bound = out_positions[i][1]-1
                    
#        if right_bound < (u+indel_pos):
#            # the sampled mutation position is located at a distance smaller than the unit length u
#            new_u = abs(right_bound - indel_pos)
#            print(new_u)
#        else:
#            # the sampled mutation position is located at a distance greater or equal than the unit length u
#            new_u = u
#            print(new_u)
        
        
#        ## find the intervals where the sampled mutatileaveson position is
#        found = dict({"start":0, "end":0, "barr":0}) # 
#        i=0
#        for i in range(len(genes_start_pos)-1):
#            print(i)
#            if genes_start_pos[i] <= indel_pos <= genes_start_pos[i+1]:
#                found["start"] = i+1
#            if genes_end_pos[i] <= indel_pos <= genes_end_pos[i+1]:
#                found["end"] = i+1
#            if barriers_pos[i] <= indel_pos <= barriers_pos[i+1]:
#                found["barr"] = i+1
#        print(found)
        ## check whether there is enough space to delete u nucleotides
        #if found["barr"] == 0:
        
#        if ( genes_start_pos[found["start"]] < (u+indel_pos) ) or ( genes_end_pos[found["end"]] < (u+indel_pos) ) or ( barriers_pos[found["barr"]] < (u+indel_pos) ):
#            # the sampled mutation position is located at a distance smaller than the unit length u
#            #pos = 
#            new_u = min(abs(genes_start_pos[found["start"]]-indel_pos-u), abs(genes_end_pos[found["end"]]-indel_pos-u), abs(barriers_pos[found["barr"]]-indel_pos-u))
#            print(new_u)
#        else:
#            # the sampled mutation position is located at a distance greater or equal than the unit length u
#            new_u = u
#            print(new_u)
        
        for i in range( len(genes_start_pos) ):        
            if genes_start_pos[i] >= indel_pos :
                # update gene start positions
                new_genes_start_pos[i] -= u
            if genes_end_pos[i] >= indel_pos :
                # update gene end positions
                new_genes_end_pos[i] -= u
            if barriers_pos[i] >= indel_pos :
                # update barrieres positions
                new_barriers_pos[i] -= u
        # update genome size
        genome_size -= u
    
    return (event_type, genome_size, new_genes_start_pos, new_genes_end_pos,
            new_barriers_pos)
        
        

    

def evolutive_event(DISCRET_STEP, inversion_proba, genome_size,
                    genes_start_pos, genes_end_pos, barriers_pos,
                    out_positions):
    """Generate an evolutive event on given genome.
    
    The event can either be a genome inversion (with proba inversion_proba),
    or an indel.
    
    Parameters
    ----------
    inversion_proba : float
        Probability for the event to be an inversion.
    genome_size : int
        Genome size in base pair.
    genes_start_pos : Numpy array
        Array of ints representing the begining position of genes.
    genes_end_pos : Numpy array
        Array of ints representing the ending position of genes.
    barriers_pos : Numpy array
        Array of ints representing the position of barriers.
    out_positions : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier.
    Returns
    -------
    genome_size : ine
        Updated value of genome_size after the event.
    genes_start_pos : Numpy array
        Updated value of genes_start_pos after the event.
    genes_end_pos : Numpy array
        Updated value of genes_end_pos after the event.
    barriers_pos : Numpy array
        Updated value of barriers_pos after the event.
    
    """
    
    event_position = sample(out_positions, genome_size, DISCRET_STEP)
    if np.random.rand() < inversion_proba:
        # The event will be an inversion
        event_position2 = sample(out_positions, genome_size, DISCRET_STEP)
        return genome_inversion(genome_size, genes_start_pos, genes_end_pos,
                                barriers_pos, min(event_position,
                                                  event_position2),
                                max(event_position, event_position2))
    else:
        # APPEL FONCTION INDEL
        return indel(DISCRET_STEP, genome_size, genes_start_pos, genes_end_pos,
                     barriers_pos)
        
    
def update_files(genome_size, genes_start_pos, genes_end_pos, barriers_pos,
                 gff_file, tss_file, tts_file, barr_file):
    """Write the initialization files for the transcription simulation.
    
    The event can either be a genome inversion (with proba inversion_proba),
    or an indel.
    
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
    gff_file : str
        Name of the .gff file (gene positions).
    tss_file : str array
        Name of the TSS file (gene start positions).
    tts_file : str
        Name of the TTS file (gene end positions).
    barr_file : str
        Name of the .dat file containing barrier positions
        
    Note
    ----
    Nothing is returned, but the files are created or updated.
    """
    
    sequence_name = gff_file[:-4]
    new_gff = open(gff_file, "w")
    new_tss = open(tss_file, "w")
    new_tts = open(tts_file, "w")
    new_barr = open(barr_file, "w")
    ### Headers
    new_gff.writelines(["##gff-version 3\n",
                        "#!gff-spec-version 1.20\n",
                        "#!processor NCBI annotwriter\n",
                        "##sequence-region " + sequence_name + " 1 "
                        + str(genome_size) + "\n",
                        sequence_name + "\tRefSeq\tregion\t1\t"
                        + str(genome_size) + "\t.\t+\t.\tID=id0;Name="
                        + sequence_name + "\n"])
    new_tss.write("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
    new_tts.write("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
    new_barr.write("prot_name\tprot_pos\n")
    ### Body
    for gene_index in range(len(genes_start_pos)):
        start = genes_start_pos[gene_index]
        end = genes_end_pos[gene_index]
        size = end - start
        if size > (genome_size / 2):
            # Gene oriented "-" but crossing the origin
            orient = "-"
        elif size > 0:
            orient = "+"
        elif size > (-genome_size / 2):
            orient = "-"
        else:
            # Gene oriented "+" but crossing the origin
            orient = "+"
        new_gff.write(sequence_name + "\tRefSeq\tgene\t" + str(start) + "\t"
                      + str(end) + "\t.\t" + orient + "\t.\tID=g1;Name=g"
                      + str(gene_index + 1) + "\n")
        new_tss.write(str(gene_index) + "\t" + orient + "\t" + str(start)
                      + "\t.2\n")
        new_tts.write(str(gene_index) + "\t" + orient + "\t" + str(end)
                      + "\t1.\n")
    for barrier in barriers_pos:
        new_barr.write("hns\t" +str(barrier) + "\n")
    ### Close
    for file in [new_gff, new_tss, new_tts, new_barr]:
        file.close()



def evolution(start, end, barr, out, genome_size, initial_expression, previous_fitness, PARAMS):
    """ Simulation of the evolution. 
        
    
    Parameters
    ----------
    start : Numpy array
        Array of ints representing the begining position of genes.
    end : Numpy array
        Array of ints representing the ending position of genes.
    barr : Numpy array
        Array of ints representing the position of barriers.
    out : Python list
        list of the open position intervals in the genome where there 
        are not any genes   
    genome_size : int
        Genome size in base pair.
    initial_expression : Numpy array
        Array of ints representing the initial number of transcripts 
        for each gene, ordered by gene ID.
    previous_fitness : float
        computed initial fitness of the invidual, following the formula:
            fitness = exp(-sum(|ln(observed_for_gene_i / target_for_gene_i)|))
    PARAMS : Python list
        list of the following parameters:
            TARGET_FREQS : Numpy array
                Array of floats representing target relative expression level 
                for each gene.
            NEXT_GEN_PARAMS : str
                Path and name of the parameter file necessary for 
                the expression_simulation function, at the next generation.
            NEXT_GEN_GFF : str
                Name of the .gff file (gene positions) at the next generation.
            NEXT_GEN_TSS : str
                Name of the TSS file (gene start positions) 
                at the next generation.
            NEXT_GEN_TTS : str
                Name of the TTS file (gene end positions)
                at the next generation.
            NEXT_GEN_BARRIERS : str
                Name of the .dat file containing barrier positions 
                at the next generation.
                
    Return
    ----
    
    """
    DISCRET_STEP = int(input("unit of length of nucleotides that is deleted or inserted: "))#60
    INVERSION_PROBA = float(input("Probability for an evolutive event to be an inversion: ")) #0.5 # Probability for an evolutive event to be an inversion.
    NB_GENERATIONS = int(input("number of generations: "))#30
    q_type = input("choose amongst the following:\n (1) q = (1 / 1000) * "
                   + "np.exp(- generation / 5) \n (2) q = Inf\n (3) q = cst\n")
    if q_type == "2":
        q = float("Inf")
    elif q_type == "3":
        q = float(input("Value of q : "))
    
    accepted_fitnesses = [previous_fitness]
    proposed_fitnesses = [previous_fitness]
    accepted_status = ["accepted"]
    all_types = ["initial"]
    generation_numbers = range(NB_GENERATIONS+1)
    for generation in generation_numbers[1:]:
        # Random evolutive event
        event_type, new_size, new_start, new_end, new_barr = (
                evolutive_event(DISCRET_STEP, INVERSION_PROBA, genome_size, start, end,
                                barr, out))
        # Update parameter files and run expression simulation.
        update_files(new_size, new_start, new_end, new_barr, PARAMS[2],
                     PARAMS[3], PARAMS[4], PARAMS[5])
        error = True
        while(error):
            try:
                new_expression = expression_simulation(PARAMS[1], "out.txt", new_start)
                new_fitness = compute_fitness(new_expression, PARAMS[0])
                error = False
            except IndexError:
                print(new_expression)
                print("ex:", [genome_size, start, end, barr])
                print("event:", event_type)
                print("new:", [new_size, new_start, new_end, new_barr])
            
        
        # Accept or reject the mutation.
        print("Generation ", end="")
        print(generation, end=":\n")
        print(event_type + " event")
        print("Fitness: ", end="")
        print(new_fitness)
        if q_type == "1":
            q = (1 / 1000) * np.exp(- generation / 5)
        
        if accept_mutation(previous_fitness, new_fitness, q):
            accepted_status.append("accepted")
            previous_fitness = new_fitness
            genome_size, start, end, barr = new_size, new_start, new_end, new_barr
            out = pos_out_from_pos_lists(start, end, barr)
        else:
            accepted_status.append("rejected")
            
        # Keep track of each event
        accepted_fitnesses.append(previous_fitness)
        proposed_fitnesses.append(new_fitness)
        all_types.append(event_type)
    
    return(accepted_fitnesses, proposed_fitnesses, accepted_status, all_types,
           generation_numbers)

INITIAL_PARAMETERS = "params.ini" # name of the file containing the initial parameter values necessary for the expression_simulation function
TARGET_FREQS = target_expression("environment.dat")
NEXT_GEN_PARAMS = "params_nextGen.ini"
NEXT_GEN_GFF = "nextGen/nextGen.gff"
NEXT_GEN_TSS = "nextGen/nextGenTSS.dat"
NEXT_GEN_TTS = "nextGen/nextGenTTS.dat"
NEXT_GEN_BARRIERS = "nextGen/nextGenProt.dat"
COLORS = {"initial" : "black", "deletion" : "red", "insertion" : "green", "inversion" : "purple"}
PARAMS = [TARGET_FREQS, NEXT_GEN_PARAMS, NEXT_GEN_GFF, NEXT_GEN_TSS, NEXT_GEN_TTS, NEXT_GEN_BARRIERS] # parameters to input in the evolution function


start, end, barr, out, size = pos_out_genes(INITIAL_PARAMETERS)
initial_expression = expression_simulation(INITIAL_PARAMETERS, "out.txt", start)
previous_fitness = compute_fitness(initial_expression, TARGET_FREQS)

OUTPUT_FILENAME = input("path and name of the output csv file: ")
(accepted_fitnesses, proposed_fitnesses, accepted_status, all_types,
 generation_numbers) = evolution(start, end, barr, out, size,
                                 initial_expression, previous_fitness, PARAMS)

plt.ylim(.9*min(proposed_fitnesses), 1.1*max(accepted_fitnesses))
plt.plot(accepted_fitnesses, linestyle="--", markersize=0, color="k", zorder=1)
plt.scatter(generation_numbers, accepted_fitnesses, alpha=1,
            c=[COLORS[event_type] for event_type in all_types], zorder=2)
plt.scatter(generation_numbers, proposed_fitnesses, marker="+",
            c=[COLORS[event_type] for event_type in all_types])

plt.show()

if OUTPUT_FILENAME:
    output_matrix = np.vstack((accepted_fitnesses,
                               proposed_fitnesses, accepted_status, all_types))
    output_df = pd.DataFrame(data=output_matrix, columns=generation_numbers,
                             index=["system fitness", "proposed fitness",
                                    "accepted", "event type"])
    output_df.to_csv(OUTPUT_FILENAME)