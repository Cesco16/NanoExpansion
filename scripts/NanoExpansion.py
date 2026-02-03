import pandas as pd
import numpy as np
import pysam
import re
import json
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


from utils import bed_ru_merge, extract_sequences, create_plot_input_files, extract_interruption_sequences, create_plot_interruption_files, split_interrupt_reads, draw_dna_gene, complementary_reverse


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
insertion_threshold = int(args.ins1)
ins_thresh = int(args.ins2)

rep_length = len(repeat_motif)
int_length = len(interrupt_motif)

## Import straglr tsv file
bd = pd.read_csv(path  + '/' +  sample + '/nanoexpansion/'+ sample + '-straglr_old.bed', sep="\t", header=0)
pd_straglr = pd.read_csv(path  + '/' +  sample + '/nanoexpansion/'+ sample + '-straglr.tsv', sep="\t", header=1)
## Import stranger annotated file
pd_stranger = pd.read_csv(path + '/' + sample + '/nanoexpansion/native'  + sample + '_rep_plot.tsv', sep="\t", header=0)
pd_stranger['POS'] = pd_straglr['start']
## Merge straglr and stranger files
merged = pd.merge(pd_straglr, pd_stranger, left_on='start', right_on='POS')
## Merging BED RU info
merged = bed_ru_merge(merged, path + '/'  + 'wf_str_repeats.bed')
## Extract sequences information as a json file
str_seq_json = extract_sequences(path + '/' +  sample + '/nanoexpansion/native' + sample + '_str_reads.bam', merged, repeat_motif)
## Extract relevant information for plot
str_identifier = create_plot_input_files(str_seq_json, sample, path)
## Create a list with all the STR fragments
## Import the created plot .csv file

plot_file = pd.read_csv(path + '/' + sample + '/nanoexpansion'  +  f"/{str_identifier}_str-content.csv")
## Create the list
read_ids = plot_file['read_id'].unique()

STR = []
for r in read_ids:

    data = plot_file[plot_file['read_id'] == r]
    data = data.sort_values(by='start')

    rep_counter = 1
    int_counter = 0
    string_sequence = []

    reps = data[data['type']=='Repeat']

    inter = data[data['type']=='Interruption']

    # --- CASE: single repeat only ---
    if len(reps) == 1:
        string_sequence.append((reps.iloc[0]['sequence'], rep_length))

    # --- CASE: multiple repeats (your original logic) ---
    elif len(reps) > 1:

        if data['type'].iloc[0] == 'Interruption':
            string_sequence.append(
                (inter.iloc[int_counter]['sequence'],
                 len(inter.iloc[int_counter]['sequence']))
            )
            int_counter += 1

        for i in np.arange(0, len(reps) - 1, 1):
            if reps.iloc[i + 1]['start'] - reps.iloc[i]['start'] > rep_length:
                string_sequence.append(
                    (reps.iloc[i]['sequence'], rep_counter * rep_length)
                )
                rep_counter = 1

                if int_counter < len(inter):
                    string_sequence.append(
                        (inter.iloc[int_counter]['sequence'],
                         len(inter.iloc[int_counter]['sequence']))
                    )
                    int_counter += 1
            else:
                rep_counter += 1

    # append last repeat block
        string_sequence.append(
            (reps.iloc[-1]['sequence'], rep_counter * rep_length)
        )

# append trailing interruption if present
    if int_counter < len(inter):
        string_sequence.append(
            (inter.iloc[int_counter]['sequence'],
             len(inter.iloc[int_counter]['sequence']))
        )

    STR.append(string_sequence)

## Remove very small insertions (<=3)


for k in np.arange(0, len(STR),1):
    removes = []
    for j in np.arange(0, len(STR[k]), 1):

        if STR[k][j][1] <= insertion_threshold: 

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


if not df_STR[df_STR['type'] == 'Interruption'].empty:
    str_sequence = df_STR[df_STR['type']=='Interruption'].iloc[0]['sequence']

    ## Extract interruption information

    interrupt_json = extract_interruption_sequences(df_STR, motif=interrupt_motif)

    create_plot_interruption_files(interrupt_json, sample, path)
    split_interrupt_reads(sample, path)

    ## Create a list of lists with interruptions fragments
    read_ids = plot_file['read_id'].unique()
    INTR = []
    for r in read_ids:
        if os.path.exists(path + '/' + sample + '/nanoexpansion' + f"/{r}_interrupt.csv"):
            data_ = pd.read_csv(path + '/' + sample + '/nanoexpansion' + f"/{r}_interrupt.csv")
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
                if int_counter < len(inter):
                    string_sequence.append((inter.iloc[int_counter]['sequence'],
                                      len(inter.iloc[int_counter]['sequence'])))
                intr.append(string_sequence)

            else:

                if len(data) == 1:

                    string_sequence.append((data.iloc[0]['sequence'],
                                      len(data.iloc[0]['sequence'])))
                    intr.append(string_sequence)
                elif len(data) > 1:

                        string_sequence.append((data.iloc[0]['sequence'],
                                      len(data.iloc[0]['sequence'])))
                        string_sequence.append((data.iloc[1]['sequence'],
                                      len(data.iloc[1]['sequence'])))
                        intr.append(string_sequence)

        INTR.append(intr)

    ## Removes very small insertions (<=3)



    for f in np.arange(0, len(INTR),1):
        for k in np.arange(0, len(INTR[f]),1):
            removes = []
            for j in np.arange(0, len(INTR[f][k]), 1):

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
        
        while len(CSTR) and CSTR.iloc[0]['Type'] == 'Other':
            CSTR = CSTR.iloc[1:]

        while len(CSTR) and CSTR.iloc[-1]['Type'] == 'Other':
            CSTR = CSTR.iloc[:-1]

        CSTR = CSTR.reset_index(drop=True)
        
        gene_sections_lengths = CSTR['Length']
        section_colors = []
        for i in np.arange(0, len(gene_sections_lengths),1):
            if CSTR['Type'][i] =='Repeat':
                section_colors.append('#E69F00')
            elif CSTR['Type'][i] =='Interruption':
                section_colors.append('#0072B2')
            else:
                section_colors.append('#4D4D4D')

        draw_dna_gene(gene_sections_lengths, section_colors, str_identifier, sample, path, motifs=[repeat_motif, interrupt_motif, 'Other'], ids=read_ids[j]) #total_length=CSTR['Length'].cumsum()[0]

    ## Extract repeats information for each read

    text_file = open(path + '/' + sample + '/nanoexpansion/' + sample + '_' + f'{str_identifier}.txt', 'w')
    for k in np.arange(0, len(complete_STR),1):
        CSTR = pd.DataFrame(complete_STR[k], columns = ['Motif','Length','Type'])
        
        while len(CSTR) and CSTR.iloc[0]['Type'] == 'Other':
            CSTR = CSTR.iloc[1:]

        while len(CSTR) and CSTR.iloc[-1]['Type'] == 'Other':
            CSTR = CSTR.iloc[:-1]

        CSTR = CSTR.reset_index(drop=True)
        
        cstr = ''
        for i in range(len(CSTR)-1, -1, -1):#for i in np.arange(0, len(CSTR), 1):
            if CSTR.iloc[i]['Type'] == 'Repeat':
                l = str(int(round(CSTR.iloc[i]['Length']/rep_length,0)))
                cstr += '('+ complementary_reverse(CSTR.iloc[i]['Motif']) + r'){}'.format(l)
            elif CSTR.iloc[i]['Type'] == 'Interruption':
                l = str(int(round(CSTR.iloc[i]['Length']/int_length,0)))
                cstr += '('+ complementary_reverse(CSTR.iloc[i]['Motif']) + r'){}'.format(l)
            else:
                l = str(int(round(CSTR.iloc[i]['Length'],0)))
                cstr += '('+ complementary_reverse(CSTR.iloc[i]['Motif']) + ')'# + r'){}'.format(l)
        print(read_ids[k], ': ', cstr)
        stt = str(read_ids[k]) + ': ' + str(cstr) + '\n'
        text_file.write(stt)

    text_file.close()

else:
    complete_STR = compact_STR
    interrupt_motif = ''
    for j in np.arange(0, len(complete_STR), 1):
        CSTR = pd.DataFrame(complete_STR[j], columns = ['Motif','Length','Type'])
        
        while len(CSTR) and CSTR.iloc[0]['Type'] == 'Other':
            CSTR = CSTR.iloc[1:]

        while len(CSTR) and CSTR.iloc[-1]['Type'] == 'Other':
            CSTR = CSTR.iloc[:-1]

        CSTR = CSTR.reset_index(drop=True)

        gene_sections_lengths = CSTR['Length']
        section_colors = []
        for i in np.arange(0, len(gene_sections_lengths),1):
            if CSTR['Type'][i] =='Repeat':
                section_colors.append('#E69F00')

            else:
                section_colors.append('#4D4D4D')

        draw_dna_gene(gene_sections_lengths, section_colors, str_identifier, sample, path, motifs=[repeat_motif, interrupt_motif, 'Other'], ids=read_ids[j]) #t

    text_file = open(path + '/' + sample + '/nanoexpansion/' + sample + '_' + f'{str_identifier}.txt', 'w')
    for k in np.arange(0, len(complete_STR), 1):
        CSTR = pd.DataFrame(complete_STR[k], columns=['Motif', 'Length', 'Type'])
        
        while len(CSTR) and CSTR.iloc[0]['Type'] == 'Other':
            CSTR = CSTR.iloc[1:]

        while len(CSTR) and CSTR.iloc[-1]['Type'] == 'Other':
            CSTR = CSTR.iloc[:-1]

        CSTR = CSTR.reset_index(drop=True)
        
        cstr = ''
        for i in range(len(CSTR)-1, -1, -1):#for i in np.arange(0, len(CSTR), 1):
            if CSTR.iloc[i]['Type'] == 'Interruption':#!= 'Other':
                print(CSTR.iloc[i]['Length'])
                l = str(int(round(CSTR.iloc[i]['Length'] / int_length, 0)))
                cstr += '(' + complementary_reverse(CSTR.iloc[i]['Motif']) + r'){}'.format(l)
                print(cstr)
            elif CSTR.iloc[i]['Type'] == 'Repeat':#!= 'Other':
                l = str(int(round(CSTR.iloc[i]['Length'] / rep_length, 0)))
                cstr += '(' + complementary_reverse(CSTR.iloc[i]['Motif']) + r'){}'.format(l)
            elif CSTR.iloc[i]['Type'] == 'Other':#else:
                l = str(int(round(CSTR.iloc[i]['Length'], 0)))
                cstr += '(' + complementary_reverse(CSTR.iloc[i]['Motif']) + ')'  # + r'){}'.format(l)
        print(read_ids[k], ': ', cstr)
        stt = str(read_ids[k]) + ': ' + str(cstr) + '\n'
        text_file.write(stt)

    text_file.close()
