import pandas as pd
import numpy as np
import pysam
import re
import json
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#from bokeh.models import BoxZoomTool, ColumnDataSource, HoverTool
#from bokeh.models import PanTool, Range1d, ResetTool, WheelZoomTool
#from dominate.tags import h3, p, span, table, tbody, td, th, thead, tr

def bed_ru_merge(merged_tsv, repeats_bed):
    """Merge RU's from STR repeats BED into straglr/stranger merged TSV."""
    bed_columns = [
        'bed_chr', 'bed_start', 'bed_end', 'bed_ru', 'bed_gene', 'bed_varid'
    ]
    remove = ['bed_start', 'bed_end', 'bed_gene']
    pd_repeats_bed = pd.read_csv(
        repeats_bed,
        sep="\t",
        header=None,
        names=bed_columns
    )
    pd_repeats_bed = pd_repeats_bed.drop(columns=remove)
    merged_tsv = pd.merge(
        merged_tsv,
        pd_repeats_bed,
        left_on='VARID',
        right_on='bed_varid'
    )
    merged_tsv = merged_tsv.drop(columns='bed_varid')
    return merged_tsv

def extract_sequences(bam, merged_tsv, repeat_motif):
    """Extract STR sequences and related info from BAM and create JSON."""
    str_seq_dict = {}

    # Need to iterate through VARID's as
    # Straglr TSV contains duplicate read ID's supporting various VARID's
    unique_varids = merged_tsv['VARID'].unique()
    #print(unique_varids)
    for varid in unique_varids:
        # Gathering supporting read_ids for the VARID
        relevant_read_ids = merged_tsv.loc[
            merged_tsv['VARID'] == varid, 'read'
        ].values

        input_bam = pysam.AlignmentFile(bam, "rb")

        for record in input_bam:
            if record.query_name in relevant_read_ids:
                read = merged_tsv['read'] == record.query_name
                print(read)
                merged_varid = merged_tsv[read & (merged_tsv['VARID'] == varid)]
                print(merged_varid)
                if not np.isnan(merged_varid['STR_NORMAL_MAX'].values[0]):###
                    strand = merged_varid['strand'].values[0]
                    chrom = merged_varid['bed_chr'].values[0]
                    repeat_start = merged_varid['read_start'].values[0]
                    size = merged_varid['size'].values[0]
                    str_normal_max = int(merged_varid['STR_NORMAL_MAX'].values[0])
                    str_pathologic_min = int(merged_varid['STR_PATHOLOGIC_MIN'].values[0])

                    # RU from Repeats BED
                    repeat_unit = repeat_motif#merged_varid['bed_ru'].values[0]

                # Create STR identifier for plot title
                    disease = merged_varid['Disease'].values[0]
                    str_identifier = f"{disease} ({varid})"

                # Collate STR Summary info into dict
                    if str_identifier not in str_seq_dict:
                        str_seq_dict[str_identifier] = {
                            'chr': chrom,
                            'repeat_unit': repeat_unit,
                            'VARID': varid,
                            'str_normal_max': str_normal_max,
                            'str_pathologic_min': str_pathologic_min,
                            'observed_reads': {}
                        }

                # Extract Haplotype
                #print(merged_varid)
                #print(relevant_read_ids)
                #print(record.query_sequence)
                #print(record.query_name)
                #print(str_identifier)
                #print(record.tags)
                    if 'HP' in dict(record.tags):
                        haplotype = dict(record.tags)['HP']
                    else:
                        haplotype = "None"

                # Extracting read and STR sequences
                # -1 as size includes repeat_start pos.
                    repeat_end = (repeat_start + size) - 1
                    read_sequence = record.query_sequence
                    if read_sequence is not None:
                        if strand == "+":
                            str_sequence = read_sequence[
                                (int(repeat_start)):(int(repeat_end) + 1)
                            ]
                        if strand == "-":
                            reversed_repeat_start = int(len(read_sequence) - repeat_end)
                            reversed_repeat_end = int(len(read_sequence) - repeat_start)
                            print(reversed_repeat_start)
                            print(reversed_repeat_end)
                            str_sequence = read_sequence[
                                (reversed_repeat_start):(reversed_repeat_end + 1)
                            ]
                    else:
                        continue
                #print(str_sequence)
                
                # Detecting the RU's and interruptions
                    repeat_unit_indexes = []
                    interruption_indexes = []

                    ru_regex = re.compile(fr'{repeat_unit}') #define a regular expression with the repeated motif
                #print(ru_regex)
                    for repeat_unit in ru_regex.finditer(str_sequence): #find the regular expression into the sequence
                    #print(repeat_unit)
                    #print(repeat_unit.start())
                    #print(repeat_unit.end())
                    #print(repeat_unit.group())
                    #print(repeat_unit.end()-repeat_unit.start())
                        repeat_unit_indexes.append((
                            repeat_unit.start(),
                            repeat_unit.end(),
                            repeat_unit.group(),
                            repeat_unit.end() - repeat_unit.start()
                        ))
                #print(repeat_unit_indexes)
                # Finding interruptions at start of seq and between RU's
                    previous_end = 0
                    for start, end, seq, length in repeat_unit_indexes:
                    #print(start)
                        if start > previous_end:
                            interruption_indexes.append((
                                previous_end,  # start of interruption
                                start,  # start of RU/end of interruption
                                str_sequence[previous_end:start],  # interruption seq
                                start - previous_end  # len of interruption
                            ))
                        previous_end = end
                    #print(interruption_indexes)
                    
                # Finding interruptions at end of sequence
                # If end of last RU is smaller than the seq len
                    str_sequence_l = len(str_sequence)
                    if previous_end < str_sequence_l:
                        interruption_indexes.append((
                            previous_end,
                            str_sequence_l,
                            str_sequence[previous_end:len(str_sequence)],
                            str_sequence_l - previous_end
                        ))

                    str_seq_dict[str_identifier]['observed_reads'].update({
                        record.query_name: {
                            "str_sequence": str_sequence,
                            "str_seq_length": str_sequence_l,
                            "haplotype": haplotype,
                            "repeat_unit_indexes": repeat_unit_indexes,
                            "interruption_indexes": interruption_indexes
                        }
                    })
                #print(str_seq_dict, '\n')
        input_bam.close()

    str_seq_json = json.dumps(str_seq_dict, indent=4)
    return str_seq_json

def truncate_interruption(interruption_seq):
    """Shorten long interruption sequences for display on Bokeh HoverTool."""
    if len(interruption_seq) > 20:
        return interruption_seq#[:21] + "..."
    else:
        return interruption_seq


#savepath = "/home/PERSONALE/francesco.casadei20/GridIon/"
#savepath = '/data/re-basecalled/'#+sample+'/'

def create_plot_input_files(str_seq_json, sample, path):
    """Extract info relevant for plots from JSON and save as CSV."""
    data = json.loads(str_seq_json)
    for str_identifier, str_data in data.items():
        print('STR_id', str_identifier)
        print('STR_data', str_data)
        rows = []

        for read_id, read_details in str_data["observed_reads"].items():
            for seq_type, indexes in {
                "Repeat": read_details["repeat_unit_indexes"],
                "Interruption": read_details["interruption_indexes"]
            }.items():

                for index in indexes:
                    rows.append([
                        str_identifier,
                        str_data['chr'],
                        str_data['VARID'],
                        str_data['repeat_unit'],
                        str_data['str_normal_max'],
                        str_data['str_pathologic_min'],
                        read_id,
                        read_details["haplotype"],
                        read_details["str_seq_length"],
                        seq_type,
                        index[2],  # sequence
                        truncate_interruption(index[2]),
                        index[0],  # start pos
                        index[1],  # end pos
                        index[3]])  # length

        df = pd.DataFrame(rows, columns=[
            "str_identifier",
            "chromosome",
            "varid",
            "repeat_unit",
            "str_normal_max",
            "str_pathologic_min",
            "read_id",
            "haplotype",
            "str_seq_length",
            "type",
            "sequence",  # RU/interruption seq
            "truncated_seq",  # seq shortened for hover tool display
            "start",  # start pos in STR seq
            "end",  # end pos in STR seq
            "length"  # len of each RU/interruption
        ])

        str_identifier = str_identifier.replace(" (", "_")
        str_identifier = str_identifier.replace(")", "")
        str_identifier = str_identifier.replace(" ", "")
        df.to_csv(path + '/' + sample + '/nanoexpansion' +  f"/{str_identifier}_str-content.csv", index=False)
        return str_identifier

def extract_interruption_sequences(str_dataset, motif):
    """Extract interruption sequences and related info from STR dataset create JSON."""
    str_seq_dict = {}

    # Need to iterate through VARID's as
    # Straglr TSV contains duplicate read ID's supporting various VARID's
    unique_varids = str_dataset['read_id'].unique() #merged_tsv['VARID'].unique()

    for varid in unique_varids:
        # Gathering supporting read_ids for the VARID
        #relevant_read_ids = merged_tsv.loc[
        #    merged_tsv['VARID'] == varid, 'read'
        #].values
        a = 0
        #input_bam = pysam.AlignmentFile(bam, "rb")
        input_file = str_dataset[str_dataset['read_id'] == varid]
        for i in np.arange(0, len(input_file),1):# in input_bam:
            #print(i)
            if input_file.iloc[i]['type'] == 'Interruption':
                #print(input_file.iloc[i]['sequence'])
                size = len(input_file.iloc[i]['sequence'])
                repeat_start = 0
            #if record.query_name in relevant_read_ids:
            #    read = merged_tsv['read'] == record.query_name
            #    merged_varid = merged_tsv[read & (merged_tsv['VARID'] == varid)]
            #    #print(merged_varid)###
            #    strand = merged_varid['strand'].values[0]
            #    chrom = merged_varid['bed_chr'].values[0]
            #    repeat_start = merged_varid['read_start'].values[0]
            #    size = merged_varid['size'].values[0]
            #    str_normal_max = int(merged_varid['STR_NORMAL_MAX'].values[0])
            #    str_pathologic_min = int(merged_varid['STR_PATHOLOGIC_MIN'].values[0])

                # RU from Repeats BED
                repeat_unit = motif#'GAG'#merged_varid['bed_ru'].values[0]

                ## Create STR identifier for plot title
                #disease = merged_varid['Disease'].values[0]
                str_identifier = varid + '_' + str(a)#f"{disease} ({varid})"

                # Collate STR Summary info into dict
                if str_identifier not in str_seq_dict:
                    str_seq_dict[str_identifier] = {
                #        'chr': chrom,
                        'repeat_unit': repeat_unit,
                #        'VARID': varid,
                #        'str_normal_max': str_normal_max,
                #        'str_pathologic_min': str_pathologic_min,
                        'observed_reads': {}
                    }


                # Extracting read and STR sequences
                # -1 as size includes repeat_start pos.
                repeat_end = (repeat_start + size) - 1
                read_sequence = input_file.iloc[i]['sequence']#record.query_sequence
                #print(input_file.iloc[i])
                #print(read_sequence)
                #print(size)
                #if read_sequence is not None:
                #    if strand == "+":
                #        str_sequence = read_sequence[
                #            (repeat_start):(repeat_end + 1)
                #        ]
                #    if strand == "-":
                #        reversed_repeat_start = len(read_sequence) - repeat_end
                #        reversed_repeat_end = len(read_sequence) - repeat_start
                #        str_sequence = read_sequence[
                #            (reversed_repeat_start):(reversed_repeat_end + 1)
                #        ]
                #else:
                #    continue
                #print(str_sequence)
                str_sequence = read_sequence
                # Detecting the RU's and interruptions
                repeat_unit_indexes = []
                interruption_indexes = []

                ru_regex = re.compile(fr'{repeat_unit}') #define a regular expression with the repeated motif
                #print(ru_regex)
                for repeat_unit in ru_regex.finditer(str_sequence): #find the regular expression into the sequence
                    repeat_unit_indexes.append((
                        repeat_unit.start(),
                        repeat_unit.end(),
                        repeat_unit.group(),
                        repeat_unit.end() - repeat_unit.start()
                    ))
                #print(repeat_unit_indexes)
                # Finding interruptions at start of seq and between RU's
                previous_end = 0
                for start, end, seq, length in repeat_unit_indexes:
                    #print(start)
                    if start > previous_end:
                        interruption_indexes.append((
                            previous_end,  # start of interruption
                            start,  # start of RU/end of interruption
                            str_sequence[previous_end:start],  # interruption seq
                            start - previous_end  # len of interruption
                        ))
                    previous_end = end
                    #print(interruption_indexes)
                    
                # Finding interruptions at end of sequence
                # If end of last RU is smaller than the seq len
                str_sequence_l = len(str_sequence)
                if previous_end < str_sequence_l:
                    interruption_indexes.append((
                        previous_end,
                        str_sequence_l,
                        str_sequence[previous_end:len(str_sequence)],
                        str_sequence_l - previous_end
                    ))
                #print(str_sequence)
                #print(str_sequence_l)
                str_seq_dict[str_identifier]['observed_reads'].update({
                    varid: { #record.query_name: {
                        "str_sequence": str_sequence,
                        "str_seq_length": str_sequence_l,
                        #"haplotype": haplotype,
                        "repeat_unit_indexes": repeat_unit_indexes,
                        "interruption_indexes": interruption_indexes
                    }
                })
                #print(str_seq_dict, '\n')
                a += 1
        #input_bam.close()

    str_seq_json = json.dumps(str_seq_dict, indent=4)
    return str_seq_json

def create_plot_interruption_files(str_seq_json, sample, path):
    """Extract info relevant for plots from JSON and save as CSV."""
    data = json.loads(str_seq_json)
    a = 0
    DF = pd.DataFrame([],columns=[
            "str_identifier",
            #"chromosome",
            #"varid",
            "repeat_unit",
            #"str_normal_max",
            #"str_pathologic_min",
            "read_id",
            #"haplotype",
            "str_seq_length",
            "type",
            "sequence",  # RU/interruption seq
            "truncated_seq",  # seq shortened for hover tool display
            "start",  # start pos in STR seq
            "end",  # end pos in STR seq
            "length"  # len of each RU/interruption
        ])
    for str_identifier, str_data in data.items():
        rows = []
        for read_id, read_details in str_data["observed_reads"].items():
            #print(read_id, read_details)
            for seq_type, indexes in {
                "Repeat": read_details["repeat_unit_indexes"],
                "Interruption": read_details["interruption_indexes"]
            }.items():
                #print(seq_type, indexes)
                for index in indexes:
                    #print(index)
                    rows.append([
                        str_identifier,
                        #str_data['chr'],
                        #str_data['VARID'],
                        str_data['repeat_unit'],
                        #str_data['str_normal_max'],
                        #str_data['str_pathologic_min'],
                        read_id,
                        #read_details["haplotype"],
                        read_details["str_seq_length"],
                        seq_type,
                        index[2],  # sequence
                        truncate_interruption(index[2]),
                        index[0],  # start pos
                        index[1],  # end pos
                        index[3]])  # length

        df = pd.DataFrame(rows, columns=[
            "str_identifier",
            #"chromosome",
            #"varid",
            "repeat_unit",
            #"str_normal_max",
            #"str_pathologic_min",
            "read_id",
            #"haplotype",
            "str_seq_length",
            "type",
            "sequence",  # RU/interruption seq
            "truncated_seq",  # seq shortened for hover tool display
            "start",  # start pos in STR seq
            "end",  # end pos in STR seq
            "length"  # len of each RU/interruption
        ])

        str_identifier = str_identifier.replace(" (", "_")
        str_identifier = str_identifier.replace(")", "")
        str_identifier = str_identifier.replace(" ", "")
        #print(df)
        #df.to_csv(f"/home/PERSONALE/francesco.casadei20/GridIon/native13204/interrupt.csv", index=False)
        a += 1
        DF = pd.concat((DF, df), axis = 0)
        #print(DF)
    DF.to_csv(path+'/' + sample+"/nanoexpansion/df_interrupt.csv", index=False)

def split_interrupt_reads(sample, path):
    int_df_path = path + '/' + sample + '/nanoexpansion/df_interrupt.csv'
    int_df = pd.read_csv(int_df_path)
    read_ids = int_df['read_id'].unique()
    
    for read in read_ids:
        #print(read)
        df = int_df[int_df['read_id']==read]
        #print(df)
        df.to_csv(path + '/' + sample + '/nanoexpansion' + f"/{read}_interrupt.csv", index=False)

def complementary_reverse(dna):
    rev = dna[::-1]
    rev=rev.replace('A','t')
    rev=rev.replace('T','a')
    rev=rev.replace('C','g')
    rev=rev.replace('G','c')
    rev=rev.replace('a','A')
    rev=rev.replace('t','T')
    rev=rev.replace('c','C')
    rev=rev.replace('g','G')
    return rev

def draw_dna_gene(sections_lengths, section_colors, str_identifier, sample, path, motifs = ['CAG','GAG','Other'], total_length=None, figsize=(10, 3), ids=None):
        #fig, ax = plt.subplots()
        reg = str_identifier.split('_')[1]
    # Definisci i colori per le diverse sezioni
    #section_colors = ['blue', 'orange', 'green', 'red']
    #fig, ax = plt.subplots()
    # Calcola la larghezza totale del gene
        if total_length is None:
            total_length = sum(sections_lengths)
            
        if ids is None:
            ids = 'a'
    
    # Disegna il gene di DNA
        current_position = 0
        fig, ax = plt.subplots(figsize=figsize)
        for length, color in zip(sections_lengths, section_colors):
        #print(current_position, current_position + length)
            ax.plot([current_position, current_position + length], [0, 0], color=color, linewidth=10)
            if color == 'red':
                ax.text(x = current_position+length/2, y = 0.01, s = str(length), size=6, ha='center', va='center')
            elif color == 'green' and length > 20:
                ax.text(x = current_position+length/2, y = -0.01, s = str(length), size=6, ha='center', va='center')
            current_position += length
    #legend_elements = [Patch(color=color, label=f'Sezione {length}') for length, color in zip(sections_lengths, section_colors)]
        legend_elements = [Patch(color=color, label=label) for color, label in zip(['red','green', 'yellow'],
                                                                                   [complementary_reverse(motifs[0]),
                                                                                   complementary_reverse(motifs[1]),
                                                                                   motifs[2]])]
        ax.legend(handles=legend_elements, loc='upper right')
        #ax.text(x = current_position, y = 0, s = str(length))
    # Nascondi assi
        ax.set_xlim(0, total_length)
        ax.set_title(f'Expansion in {reg} ({ids})')
        ax.axis('off')
    # Mostra il grafico
        plt.savefig(path + '/' + sample + '/nanoexpansion/' + ids + f'_{str_identifier}.png')
        plt.show()

