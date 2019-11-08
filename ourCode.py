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
    out = pos_out_from_pos_lists(start, end, barr)
    return start, end, barr, out

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
        print(i)
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
            
            
        print(i, " ", mut_pos)
        
        i-=1
        if i==0:
            raise RuntimeError("A mutation position that that satisfies the distance to the bounds condition was not found")
    
    return mut_pos
        


def indel(u, genome_size, genes_start_pos, genes_end_pos, barriers_pos, out_positions):
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
    out_positions : Numpy array
        2-D array of ints. Each line represents an open interval containing
        no gene nor barrier.
        
    Returns
    -------
    new_genes_start_pos : Numpy array
        Updated value of genes_start_pos after inversion.
    new_genes_end_pos : Numpy array
        Updated value of genes_end_pos after inversion.
    new_barriers_pos : Numpy array
        Updated value of barriers_pos after inversion.
    genome_size: int
        Updated value of the genome size
    
    """
    ### sample the indel position
    indel_pos = 27613 #sample(out_positions, genome_size)
    print(indel_pos)
    
    ### initialization of the new positions
    new_genes_start_pos = np.copy(genes_start_pos)
    new_genes_end_pos = np.copy(genes_end_pos)
    new_barriers_pos = np.copy(barriers_pos)
    
    ### choose whether it is an insertion or a deletion
    p = np.random.uniform(0,1) # draw a random number between 0 and 1
    if p<.5:
        # it is an insertion
        print("insert")
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
        print("delete")
        
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
    
    return new_genes_start_pos, new_genes_end_pos, new_barriers_pos, genome_size
        
        

    

def evolutive_event(inversion_proba, genome_size, genes_start_pos,
                    genes_end_pos, barriers_pos, out_positions):
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
    genes_start_pos : Numpy array
        Updated value of genes_start_pos after the event.
    genes_end_pos : Numpy array
        Updated value of genes_end_pos after the event.
    barriers_pos : Numpy array
        Updated value of barriers_pos after the event.
    """
    
    event_position = sample(out_positions, genome_size)
    if np.random.rand() < inversion_proba:
        # The event will be an inversion
        event_position2 = sample(out_positions)
        return genome_inversion(genome_size, genes_start_pos, genes_end_pos,
                                barriers_pos,  min(event_position,
                                                   event_position2),
                                max(event_position, event_position2))
    else:
        # APPEL FONCTION INDEL
        return genes_start_pos, genes_end_pos, barriers_pos
        
    
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


start, end, barr, out = pos_out_genes("params.ini")
TARGET_FREQS = target_expression("environment.dat")
INVERSION_PROBA = 0.5 # Probability for an evolutive event to be an inversion.
NEXT_GEN_GFF = "nextGen/nextGen.gff"
NEXT_GEN_TSS = "nextGen/nextGenTSS.dat"
NEXT_GEN_TTS = "nextGen/nextGenTTS.dat"
NEXT_GEN_BARRIERS = "nexten/nextGenProt.dat"
print (TARGET_FREQS)
#initial_expression = expression_simulation("params.ini", "out.txt")
#print(initial_expression)
#print(compute_fitness(initial_expression, TARGET_FREQS))
#print(genome_inversion(30000, start, end, barr, 2020, 13000))
#print (pos_out_from_pos_lists(start, end, barr))

indel_start, indel_end, indel_barr = indel(60, 30000, start, end, barr, out)

