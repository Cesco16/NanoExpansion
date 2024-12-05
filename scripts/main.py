import pandas as pd
import numpy as np
import pysam
import re
import json
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from bokeh.models import BoxZoomTool, ColumnDataSource, HoverTool
from bokeh.models import PanTool, Range1d, ResetTool, WheelZoomTool
from dominate.tags import h3, p, span, table, tbody, td, th, thead, tr

from utils import bed_ru_merge, extract_sequences, create_plot_input_files

## Global variables
sample = 'vac06122023'#'native9411'#'native13204'#'native10498'
repeat_motif = 'CAG' #'TAAAA'#'CAG'
interrupt_motif = 'CAA' #'TGGAA'#'GAG'# 'GAG'
rep_length = len(repeat_motif)
int_length = len(interrupt_motif)

path = "/home/PERSONALE/francesco.casadei20/GridIon/" #path where bam and vcf files are saved

## Import straglr tsv file
pd_straglr = pd.read_csv(path+ sample +'/'+ sample +'_straglr.tsv', sep="\t", header=1)
## Import stranger annotated file
pd_stranger = pd.read_csv(path+ sample +'/'+ sample + '_rep_plot.tsv', sep="\t", header=0)
## Merge straglr and stranger files
merged = pd.merge(pd_straglr, pd_stranger, left_on='start', right_on='POS')

## Merging BED RU info
merged = bed_ru_merge(merged, '/home/PERSONALE/francesco.casadei20/GridIon/wf_str_repeats.bed')
## Extract sequences information as a json file
str_seq_json = extract_sequences('/home/PERSONALE/francesco.casadei20/GridIon/'+ sample +'/'+ sample +'_str_reads.bam', merged)
## Extract relevant information for plot
str_identifier = create_plot_input_files(str_seq_json)

## Create a list with all the STR fragments
## Import the created plot .csv file
plot_file = pd.read_csv(path+sample+f"/{str_identifier}_str-content.csv")


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

