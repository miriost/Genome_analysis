# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 13:21:44 2018

@author: mirio
"""

import pandas as pd
import numpy as np
import itertools

def number_of_records(file):
    #number of records in file. A record in a FASTA file is defined as a single-line header, followed by lines of sequence data.
    f = open(file, 'r').read()
    return f.count(">") 

def parse_file(file, outfile, reading_frame = 1, n =3):
    header = None
    seqs_df = pd.DataFrame(index=range(number_of_records(file)),columns=['id', 'sequence', 'len','ORF_start', 'ORF_end', 'ORF_len'])
#    seqs_lens = {}
    i=0
#    seq_list = {}
#    identifiers_list = []
    with open(file, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                if header is not None:
                    seqs_df['id'][i] = header.split(' ')[0][1:]  
                    print('~Writing to file~ line: ' + str(i) + ' |length: ' + str(seqs_df['len'][i]) + ' |id: '+ str(seqs_df['id'][i]))                              
                    fout.write(header + '\t' + '|length: '+ str(seqs_df['len'][i]) + '\n' )
                    i+=1
                header = line
                
            else:
                text = line.rstrip('\n')
                if np.isnan(seqs_df['len'][i]):
                    #sequence wasn't touched yet
                    seqs_df['len'][i] = len(text)
                    seqs_df['sequence'][i] =text
                    
                else:
                    seqs_df['len'][i] = seqs_df['len'][i] + len(text)
                    seqs_df['sequence'][i] = seqs_df['sequence'][i] + text
                    #print('line: ' + str(i) + ' current len: ' + str(seqs_df['len'][i]))
         #handle last sequence
        seqs_df['id'][i] = header.split(' ')[0][1:]   
        fout.write(header + '\t' + '|length: '+ str(seqs_df['len'][i]) + '\n' )
#==============================================================================
#         if header is not None:
#             fout.write(header + '\t' + 'length: ' + str(seqs_lens[i]))
#         
#==============================================================================
        fout.write('~~~~~~SUMMARY~~~~~~~~\n Total number of records:' + str(i) + '\n')
        print('~~~~~~SUMMARY~~~~~~~~\n Total number of records:' + str(i) + '\n')
        # Find shortest and logest sequences
        maxi = seqs_df['len'][1]
        maxi_index = [1]
        mini = seqs_df['len'][1]
        mini_index = [1]
        for key in range(2, i):
            if seqs_df['len'][key] > maxi:
                maxi = seqs_df['len'][key]
                maxi_index = [key]
            elif seqs_df['len'][key] == maxi:
                maxi_index.append[key]
            if seqs_df['len'][key] < mini:
                mini = seqs_df['len'][key]
                mini_index = [key]
            elif seqs_df['len'][key] == mini:
                mini_index.append[key]
                
        fout.write('Longest sequence: ' + str(maxi) + ' Indexes: '+ str(maxi_index)+ ' identefiers: ' + str(seqs_df['id'][maxi_index]) + '\n')
        fout.write('Shortest sequence: ' + str(mini) + ' Indexes: '+ str(mini_index)+' identefiers: ' + str(seqs_df['id'][mini_index]) + '\n')
        print('Longest sequence: ' + str(maxi) + ' Indexes: '+ str(maxi_index)+ ' identefiers: ' + str(seqs_df['id'][maxi_index]) + '\n')
        print('Shortest sequence: ' + str(mini) + ' Indexes: '+ str(mini_index)+' identefiers: ' + str(seqs_df['id'][mini_index]) + '\n')
        
        seqs_df.to_csv('/home/miri-o/Documents/fasta_files/df_file.csv')
        return seqs_df
        
def find_ORF(input_df, reading_frame=3, id = None):
    """
    recieves a sequence list and a reading frame, and answers the questions:
    what is the length of the longest ORF in the file? 
    What is the identifier of the sequence containing the longest ORF? 
    For a given sequence identifier, what is the longest ORF contained 
    in the sequence represented by that identifier? What is the starting
    position of the longest ORF in the sequence that contains it? 
    The position should indicate the character number in the sequence.
    """
    start_codon = 'ATG'
    stop_codon = ['TAA', 'TAG', 'TGA']
    start= -1
    stop = -1
    n=3
    longest_ORF = 0
    longest_ORF_id = None
    ORF_start=-1
    ORF_end=-1

    # find the londest ORF in each sequence
    
    # first, if we have the id, we need to find the sequence with this id and check only there
    if id is not None:
        longest_ORF_id = id
        rawseq = input_df['sequence'][input_df['id']==id]
        seq = split_ngrams_no_repetitions(rawseq, n, reading_frame)
        if start_codon in seq:
            for i in range(seq.count(start)-1):
                start = seq.index(start_codon)
                stop = min([seq.index(codon) if codon in seq[start:] else 1e6 for codon in stop_codon])
                if stop!=1e6 and (stop-start)*n>longest_ORF:
                    longest_ORF = (stop-start)*3
                    ORF_start, ORF_end =start*n+reading_frame+1, stop*n+reading_frame+1
                else:
                    continue

    else: #run over all sequences
        for index in range(len(input_df['sequence'])):
            rawseq = input_df['sequence'][index]
            seq = split_ngrams_no_repetitions(rawseq, n, reading_frame)
            
            if start_codon in seq:
                start = 0
                for i in range(seq.count(start_codon)):
                    start = seq.index(start_codon, start)
                    stop = min([seq.index(codon, start) if codon in seq[start:] else 1e6 for codon in stop_codon])
                    if stop==1e6:
                        break
                    elif (stop-start)*n>longest_ORF:
                        longest_ORF = (stop-start)*3
                        ORF_start, ORF_end =start*n+reading_frame+1, stop*n+reading_frame+1
                        longest_ORF_id = input_df['id'][index]
                        print('found new ORF, length: %d, start: %d, end: %d, id: %s'%(longest_ORF, ORF_start, ORF_end, longest_ORF_id))
                    start = stop
                    
                    
            #    print(str(seq[start]) + ' ' + str(input_df['sequence'][index][(start)*n+(reading_frame-1):((start)*n+(reading_frame-1)+3)])
            #print(input_df['id'][index] + ' stop = ' + str(stop) +' ' + str(seq[stop]) + ' |original stop: ' + str((stop)*n+(reading_frame-1)) + ' ' +str(input_df['sequence'][index][(stop)*n+(reading_frame-1):((stop)*n+(reading_frame-1)+3)]))

    return ORF_start, ORF_end, longest_ORF, longest_ORF_id
        
        
def split_ngrams_no_repetitions(seq, n, reading_frame):
    """
    'acccgtgtctgg', n=3, reading frame = 1: ['acc', 'cgt', 'gtc', 'tgg']
    reading frame = 2: ['ccc', 'gtg', 'tct']
    reading frame = 3: ['ccg', 'tgt', 'ctg']
    """
    a, b, c = zip(*[iter(seq)]*n), zip(*[iter(seq[1:])]*n), zip(*[iter(seq[2:])]*n)
    str_ngrams = []
    for ngrams in [a,b,c]:
        x = []
#        if reading_frame>1:
#            x.append(seq[0:(reading_frame-1)])
        for ngram in ngrams:
            x.append("".join(ngram))
        str_ngrams.append(x)
    return str_ngrams[reading_frame-1]

def count_repeats(df, n=3):        
    all_ngrams = [list(n) for n in itertools.product('ACTG', repeat = 3)]
    words = [''.join(ngram) for ngram in all_ngrams]
    #print(words)
    word_dict = {}
    for word in words:
        repcount = 0
        for seq in df['sequence']:
            repcount = repcount+seq.count(word)
        word_dict[word] = repcount
    #print(list(word_dict))             
    return word_dict
        
    
    
                
                
if __name__ == '__main__':
    df = parse_file("/home/miri-o/Documents/fasta_files/dna2.fasta", "/home/miri-o/Documents/fasta_files/dna2_out.xls")
    for i in range(3):
        ORF_start, ORF_end, longest_ORF, ORF_id = find_ORF(df, reading_frame=i+1)

        print('Longest ORF with reading frame %d is %s, start: %d, end: %d, length: %d' %(i+1, ORF_id, ORF_start, ORF_end, longest_ORF))
    repeats = count_repeats(df, n=3)
    print('most frequent repeat: '+ max(repeats, key=repeats.get) + ' appears %d times' %max(repeats.values()))
        
#==============================================================================
# 
#==============================================================================
# filename = "dna.example.fasta"
# parse_file(filename, 'dna_out.txt')
#==============================================================================
# 
#==============================================================================
