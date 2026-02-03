import pandas as pd
import numpy as np
import pysam
import re
import json
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
    for varid in unique_varids:
        # Gathering supporting read_ids for the VARID
        relevant_read_ids = merged_tsv.loc[
            merged_tsv['VARID'] == varid, 'read'
        ].values

        input_bam = pysam.AlignmentFile(bam, "rb")

        for record in input_bam:
            if record.query_name in relevant_read_ids:
                read = merged_tsv['read'] == record.query_name
                merged_varid = merged_tsv[read & (merged_tsv['VARID'] == varid)]
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
                            str_sequence = read_sequence[
                                (reversed_repeat_start):(reversed_repeat_end + 1)
                            ]
                    else:
                        continue
                
                # Detecting the RU's and interruptions
                    repeat_unit_indexes = []
                    interruption_indexes = []

                    ru_regex = re.compile(fr'{repeat_unit}') #define a regular expression with the repeated motif
                    for repeat_unit in ru_regex.finditer(str_sequence): #find the regular expression into the sequence
                        repeat_unit_indexes.append((
                            repeat_unit.start(),
                            repeat_unit.end(),
                            repeat_unit.group(),
                            repeat_unit.end() - repeat_unit.start()
                        ))
                # Finding interruptions at start of seq and between RU's
                    previous_end = 0
                    for start, end, seq, length in repeat_unit_indexes:
                        if start > previous_end:
                            interruption_indexes.append((
                                previous_end,  # start of interruption
                                start,  # start of RU/end of interruption
                                str_sequence[previous_end:start],  # interruption seq
                                start - previous_end  # len of interruption
                            ))
                        previous_end = end
                    
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
        input_bam.close()

    str_seq_json = json.dumps(str_seq_dict, indent=4)
    return str_seq_json

def truncate_interruption(interruption_seq):
    """Shorten long interruption sequences for display on Bokeh HoverTool."""
    if len(interruption_seq) > 20:
        return interruption_seq#[:21] + "..."
    else:
        return interruption_seq


def create_plot_input_files(str_seq_json, sample, path):
    """Extract info relevant for plots from JSON and save as CSV."""
    data = json.loads(str_seq_json)
    for str_identifier, str_data in data.items():
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
        df.to_csv(path + '/'+ sample + '/nanoexpansion' +  f"/{str_identifier}_str-content.csv", index=False)
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
                str_sequence = read_sequence
                # Detecting the RU's and interruptions
                repeat_unit_indexes = []
                interruption_indexes = []

                ru_regex = re.compile(fr'{repeat_unit}') #define a regular expression with the repeated motif
                for repeat_unit in ru_regex.finditer(str_sequence): #find the regular expression into the sequence
                    repeat_unit_indexes.append((
                        repeat_unit.start(),
                        repeat_unit.end(),
                        repeat_unit.group(),
                        repeat_unit.end() - repeat_unit.start()
                    ))
                # Finding interruptions at start of seq and between RU's
                previous_end = 0
                for start, end, seq, length in repeat_unit_indexes:
                    if start > previous_end:
                        interruption_indexes.append((
                            previous_end,  # start of interruption
                            start,  # start of RU/end of interruption
                            str_sequence[previous_end:start],  # interruption seq
                            start - previous_end  # len of interruption
                        ))
                    previous_end = end
                    
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
                    varid: { #record.query_name: {
                        "str_sequence": str_sequence,
                        "str_seq_length": str_sequence_l,
                        #"haplotype": haplotype,
                        "repeat_unit_indexes": repeat_unit_indexes,
                        "interruption_indexes": interruption_indexes
                    }
                })
                a += 1

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
            for seq_type, indexes in {
                "Repeat": read_details["repeat_unit_indexes"],
                "Interruption": read_details["interruption_indexes"]
            }.items():
                for index in indexes:
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
        a += 1
        DF = pd.concat((DF, df), axis = 0)
    DF.to_csv(path+ '/' + sample+"/nanoexpansion/df_interrupt.csv", index=False)

def split_interrupt_reads(sample, path):
    int_df_path = path + '/' + sample + '/nanoexpansion/df_interrupt.csv'
    int_df = pd.read_csv(int_df_path)
    read_ids = int_df['read_id'].unique()
    
    for read in read_ids:
        df = int_df[int_df['read_id']==read]
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

from matplotlib.patches import Patch, FancyArrow

def draw_dna_gene(sections_lengths,
                  section_colors,
                  str_identifier,
                  sample,
                  path,
                  motifs=['CAG','GAG','Other'],
                  total_length=None,
                  figsize=(12,3),
                  ids=None,
                  show_lengths=True):

    reg = str_identifier.split('_')[1]

    # ---- Reverse orientation ----
    sections_lengths = sections_lengths[::-1]
    section_colors   = section_colors[::-1]

    if total_length is None:
        total_length = sum(sections_lengths)

    if ids is None:
        ids = "read"

    fig, ax = plt.subplots(figsize=figsize)

    current_position = 0

    # ---- Draw segments ----
    for length, color in zip(sections_lengths, section_colors):

        ax.plot(
            [current_position, current_position + length],
            [0,0],
            color=color,
            linewidth=12,
            solid_capstyle="butt"
        )

        if show_lengths and length > 20:
            ax.text(
                current_position + length/2,
                0.05,
                str(int(length)),
                ha='center',
                va='bottom',
                fontsize=11
            )

        current_position += length

    # ---- 5' -> 3' arrow ----
    #ax.annotate(
    #    "5'",
    #    xy=(0, -0.25),
    #    ha='center',
    #    fontsize=12
    #)

    #ax.annotate(
    ##    "3'",
    #    xy=(total_length, -0.25),
    #    ha='center',
    #    fontsize=12
    #)

    #ax.arrow(
    #    0, -0.12,
    #    total_length, 0,
    #    length_includes_head=True,
    #    head_width=0.03,
    #    head_length=total_length*0.02,
    #    fc="black",
    #    ec="black"
    #)

    # ---- Total length ----
    ax.text(
        total_length/2,
        0.25,
        f"Total length: {int(total_length)} bp",
        ha='center',
        fontsize=12
    )

    # ---- Legend ----
    legend_elements = [
    Patch(facecolor='#E69F00', label=complementary_reverse(motifs[0])),
    Patch(facecolor='#0072B2', label=complementary_reverse(motifs[1])),
    Patch(facecolor='#4D4D4D', label=motifs[2])
]

    ax.legend(handles=legend_elements, frameon=False, loc="upper right")

    # ---- Axes formatting ----
    ax.set_xlim(0, total_length)
    ax.set_ylim(-0.3, 0.4)
    ax.set_title(f'Expansion in {reg} ({ids})', fontsize=14)
    ax.axis("off")

    # ---- Save ----
    plt.tight_layout()
    plt.savefig(path + '/' + sample + '/nanoexpansion/' + ids + f'_{str_identifier}.png', dpi=600)
    plt.savefig(path + '/' + sample + '/nanoexpansion/' + ids + f'_{str_identifier}.svg', dpi=600)
    plt.savefig(path + '/' + sample + '/nanoexpansion/' + ids + f'_{str_identifier}.tiff', dpi=600)
    plt.show()
