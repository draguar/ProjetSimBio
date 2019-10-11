# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import numpy as np
import os

#
#os.system("python3 start_simulation.py params.ini > out.txt")
#with open("out.txt", "r") as out:
#    for line in out:
#        print(line)
        
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