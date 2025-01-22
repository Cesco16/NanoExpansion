import pandas as pd
import numpy as np
import pysam
import re
import json
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#from bokeh.models import BoxZoomTool, ColumnDataSource, HoverTool
#from bokeh.models import PanTool, Range1d, ResetTool, WheelZoomTool
#from dominate.tags import h3, p, span, table, tbody, td, th, thead, tr

from utils import bed_ru_merge, extract_sequences, create_plot_input_files, extract_interruption_sequences, create_plot_interruption_files, split_interrupt_reads


parser = argparse.ArgumentParser(description='missing_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sample')
parser.add_argument('--repeat', default='CAG')
parser.add_argument('--interruption', default='CAA')
parser.add_argument('--path')
parser.add_argument('--ins1', default=3)
parser.add_argument('--ins2', default=1)
args = parser.parse_args()

sample = args.sample
repeat_motif = args.repeat
interrupt_motif = args.interruption
path = args.path
insertion_threshold = args.ins1
ins_thresh = args.ins2

## Global variables
#sample = 'vac06122023'#'native9411'#'native13204'#'native10498'
#repeat_motif = 'CAG' #'TAAAA'#'CAG'
#interrupt_motif = 'CAA' #'TGGAA'#'GAG'# 'GAG'

rep_length = len(repeat_motif)
int_length = len(interrupt_motif)

#path = "/home/PERSONALE/francesco.casadei20/GridIon/" #path where bam and vcf files are saved

## Import straglr tsv file
pd_straglr = pd.read_csv(path + sample + '/' + sample +'-straglr.tsv', sep="\t", header=1)
## Import stranger annotated file
pd_stranger = pd.read_csv(path + sample + '/' + sample + '_rep_plot.tsv', sep="\t", header=0)
## Merge straglr and stranger files
merged = pd.merge(pd_straglr, pd_stranger, left_on='start', right_on='POS')

## Merging BED RU info
merged = bed_ru_merge(merged, path + 'wf_str_repeats.bed')
## Extract sequences information as a json file
str_seq_json = extract_sequences(path + sample + '/' + sample + '_str_reads.bam', merged, repeat_motif)
## Extract relevant information for plot
str_identifier = create_plot_input_files(str_seq_json, path, sample)

## Create a list with all the STR fragments
## Import the created plot .csv file
plot_file = pd.read_csv(path + sample + f"/{str_identifier}_str-content.csv")
## Create the list
read_ids = plot_file['read_id'].unique()
STR = []
for r in read_ids:
    data = plot_file[plot_file['read_id'] == r]
    data = data.sort_values(by='start')
    rep_counter = 1
    int_counter = 0
    string_sequence = []
    #print(len(data[data['type']=='Repeat']))
    reps = data[data['type']=='Repeat']
    inter = data[data['type']=='Interruption']
    if data['type'].iloc[0] == 'Interruption':
        string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
        int_counter+=1
        for i in np.arange(0, len(reps)-1, 1):
        #print(data.iloc[i]['start'])
            if reps.iloc[i+1]['start'] - reps.iloc[i]['start'] > rep_length:
                string_sequence.append((reps.iloc[i]['sequence'],rep_counter*rep_length)) #*3
                rep_counter = 1
                string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
                int_counter+=1
            else:
                rep_counter+=1
    else:
        for i in np.arange(0, len(reps)-1, 1):
        #print(data.iloc[i]['start'])
            if reps.iloc[i+1]['start'] - reps.iloc[i]['start'] > rep_length:
                string_sequence.append((reps.iloc[i]['sequence'],rep_counter*rep_length)) #*3
                rep_counter = 1
                string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
                int_counter+=1
            else:
                rep_counter+=1
    string_sequence.append((reps.iloc[i]['sequence'],rep_counter*rep_length)) #*3
    #if int_counter < len(inter):
    #    string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  #len(inter.iloc[int_counter]['sequence'])))
    #print(reps.iloc[i]['end'])
    #print(inter.iloc[int_counter]['end'])
    #if reps.iloc[i]['end'] < inter.iloc[int_counter-1]['end']:
    if int_counter < len(inter):
        string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
    STR.append(string_sequence)

## Remove very small insertions (<=3)

#insertion_threshold = 3

for k in np.arange(0, len(STR),1):
    removes = []
    for j in np.arange(0, len(STR[k]), 1):
        #print(j,k)
        #print(STR[k][j][1])
        if STR[k][j][1] <= insertion_threshold:          #3 #STR[k][j][0] != MOTIF and 
            #print(j)
            #print(STR[k][j][1])
            removes.append(j)
    for i in sorted(removes, reverse=True):
        del STR[k][i]


## Sum up new compact STR (after removal of small deletions)

MOTIF = repeat_motif #'CAG'#merged['repeat_unit'][0]
compact_STR = []
for k in np.arange(0, len(STR),1):
    cstr = []
    cag = 0
    inser = 0
    inter = ''
    for j in np.arange(1, len(STR[k]), 1):
        if STR[k][j-1][0] == MOTIF:
            cag += STR[k][j-1][1]
            inser = 0
            inter = ''
            if STR[k][j][0] != MOTIF:
                cstr.append((MOTIF, cag, 'Repeat'))
        else:
            inser += STR[k][j-1][1]
            inter += STR[k][j-1][0]
            cag = 0
            if STR[k][j][0] == MOTIF:
                cstr.append((inter, inser, 'Interruption'))
    if STR[k][len(STR[k])-1][0] == MOTIF:
        cag += STR[k][len(STR[k])-1][1]
        cstr.append((MOTIF, cag, 'Repeat'))
    else:
        inser += STR[k][len(STR[k])-1][1]
        inter += STR[k][len(STR[k])-1][0]
        cstr.append((inter, inser, 'Interruption'))
    compact_STR.append(cstr)

## Create a dataframe with the STR motif

rows = []
for i in np.arange(0, len(compact_STR), 1):
    for j in np.arange(0, len(compact_STR[i]),1):
        rows.append((read_ids[i], compact_STR[i][j][2], compact_STR[i][j][0], compact_STR[i][j][1]))

        df_STR = pd.DataFrame(rows, columns=[
            "read_id",
            "type",
            "sequence",  # RU/interruption seq
            "length"  # len of each RU/interruption
        ])

## Analize the interruptions motif
##motifs = ['CAA','GAG','CCG','CGG','CTG','GCC','CAGCGG','CAGCAGCGG']

repeat_unit = interrupt_motif #'GAG'
repeat_unit_indexes = []
interruption_indexes = []
str_sequence = df_STR[df_STR['type']=='Interruption'].iloc[0]['sequence']

## Extract interruption information

interrupt_json = extract_interruption_sequences(df_STR, motif=interrupt_motif)

create_plot_interruption_files(interrupt_json, path)
split_interrupt_reads(path + sample + "/df_interrupt.csv", path)

## Create a list of lists with interruptions fragments
read_ids = plot_file['read_id'].unique()
INTR = []
for r in read_ids:#[0:3]:
    if os.path.exists(path + sample + f"/{r}_interrupt.csv"): 
        data_ = pd.read_csv(path + sample + f"/{r}_interrupt.csv")
    else:
        data_ = pd.DataFrame()
    strs = []
    intr = []
    for i in np.arange(0, len(data_),1):
        strs.append(data_['str_identifier'][i])
    for a in np.unique(strs):
        data = data_[data_['str_identifier'] == a]
        data = data.sort_values(by='start')
        rep_counter = 1
        int_counter = 0
        string_sequence = []
        reps = data[data['type']=='Repeat']
        inter = data[data['type']=='Interruption']
        #print(reps)
        #print(len(reps))
        if len(reps) > 1:
            if data['type'].iloc[0] == 'Interruption':
                string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
                int_counter+=1
                for i in np.arange(0, len(reps)-1, 1):
                    if reps.iloc[i+1]['start'] - reps.iloc[i]['start'] > int_length:
                        string_sequence.append((reps.iloc[i]['sequence'],rep_counter*int_length)) #*3
                        rep_counter = 1
                        string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
                        int_counter+=1
                    else:
                        rep_counter+=1
            else:
                for i in np.arange(0, len(reps)-1, 1):
                    if reps.iloc[i+1]['start'] - reps.iloc[i]['start'] > int_length:
                        string_sequence.append((reps.iloc[i]['sequence'],rep_counter*int_length)) #*3
                        rep_counter = 1
                        string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
                        int_counter+=1
                    else:
                        rep_counter+=1
            string_sequence.append((reps.iloc[i]['sequence'],rep_counter*int_length)) #*3
    #if reps.iloc[i]['end'] < inter.iloc[int_counter-1]['end']:
            if int_counter < len(inter):
                string_sequence.append((inter.iloc[int_counter]['sequence'], 
                                  len(inter.iloc[int_counter]['sequence'])))
            intr.append(string_sequence)
            #print(intr)
        else:
            #print(len(data))
            if len(data) == 1:
                #print(a)
                string_sequence.append((data.iloc[0]['sequence'], 
                                  len(data.iloc[0]['sequence'])))
                intr.append(string_sequence)
            elif len(data) > 1:
                    #print(a)
                    string_sequence.append((data.iloc[0]['sequence'], 
                                  len(data.iloc[0]['sequence'])))
                    string_sequence.append((data.iloc[1]['sequence'], 
                                  len(data.iloc[1]['sequence'])))
                    intr.append(string_sequence)
                    #print(intr)
    INTR.append(intr)

## Removes very small insertions (<=3)

#ins_thresh = 1

for f in np.arange(0, len(INTR),1):
    for k in np.arange(0, len(INTR[f]),1):
        removes = []
        for j in np.arange(0, len(INTR[f][k]), 1):
        #print(j,k)
        #print(STR[k][j][1])
            if INTR[f][k][j][1] <= ins_thresh: #1
            #print(j)
                removes.append(j)
        for i in sorted(removes, reverse=True):
            del INTR[f][k][i]

## sum up new compact INTR
MOTIF = interrupt_motif #'GAG'#merged['repeat_unit'][0]
compact_INTR = []
for f in np.arange(0, len(INTR),1):
    compact_intr = []
    for k in np.arange(0, len(INTR[f]),1):
        cstr = []
        gag = 0
        inser = 0
        inter = ''
        for j in np.arange(1, len(INTR[f][k]), 1):
            if INTR[f][k][j-1][0] == MOTIF:
                gag += INTR[f][k][j-1][1]
                inser = 0
                inter = ''
                if INTR[f][k][j][0] != MOTIF:
                    cstr.append((MOTIF, gag, 'Interruption'))
            else:
                inser += INTR[f][k][j-1][1]
                inter += INTR[f][k][j-1][0]
                gag = 0
                if INTR[f][k][j][0] == MOTIF:
                    cstr.append((inter, inser, 'Other'))
        if INTR[f][k][len(INTR[f][k])-1][0] == MOTIF:
            gag += INTR[f][k][len(INTR[f][k])-1][1]
            cstr.append((MOTIF, gag, 'Interruption'))
        else:
            inser += INTR[f][k][len(INTR[f][k])-1][1]
            inter += INTR[f][k][len(INTR[f][k])-1][0]
            cstr.append((inter, inser, 'Other'))
        compact_intr.append(cstr)
    compact_INTR.append(compact_intr)


complete_STR = []
for i in np.arange(0, len(compact_STR), 1):
    comp_str = []
    k = 0
    for j in np.arange(0, len(compact_STR[i]), 1):
        if compact_STR[i][j][2] == 'Repeat':
            comp_str.append(compact_STR[i][j])
        elif compact_STR[i][j][2] == 'Interruption':
            for f in np.arange(0, len(compact_INTR[i][k]), 1):
                comp_str.append(compact_INTR[i][k][f])
            k += 1
    complete_STR.append(comp_str)


## Draw STR motif

for j in np.arange(0, len(complete_STR), 1):
    CSTR = pd.DataFrame(complete_STR[j], columns = ['Motif','Length','Type'])

    gene_sections_lengths = CSTR['Length']                                                                             
    section_colors = []
    for i in np.arange(0, len(gene_sections_lengths),1):
        if CSTR['Type'][i] =='Repeat':
            section_colors.append('red')
        elif CSTR['Type'][i] =='Interruption':
            section_colors.append('green')
        else:
            section_colors.append('yellow')

    draw_dna_gene(gene_sections_lengths, section_colors, path, motifs=[repeat_motif, interrupt_motif, 'Other'], ids=read_ids[j]) #total_length=CSTR['Length'].cumsum()[0]

## Extract repeats information for each read

text_file = open(path + sample + '/' + f'{str_identifier}.txt', 'w')
for k in np.arange(0, len(complete_STR),1):
    CSTR = pd.DataFrame(complete_STR[k], columns = ['Motif','Length','Type'])
    cstr = ''
    for i in np.arange(0, len(CSTR), 1):
        if CSTR.iloc[i]['Type'] != 'Other':
            l = str(int(round(CSTR.iloc[i]['Length']/rep_length,0)))
            cstr += '('+ complementary_reverse(CSTR.iloc[i]['Motif']) + r'){}'.format(l)
        else:
            l = str(int(round(CSTR.iloc[i]['Length'],0)))
            cstr += '('+ complementary_reverse(CSTR.iloc[i]['Motif']) + ')'# + r'){}'.format(l)
    print(read_ids[k], ': ', cstr)
    stt = str(read_ids[k]) + ': ' + str(cstr) + '\n'
    text_file.write(stt)
    #l = int(l)
    #l = str(l).encode("utf-64")#.decode("utf-8")
    #print(l)
    #cstr += '('+ CSTR.iloc[i]['Motif'] + r'){}'.format(l)
text_file.close()
