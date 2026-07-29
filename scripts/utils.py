"""
utils.py
========
Utility functions for NanoExp3 (NanoExpansion).

Contains functionality for:
- Read extraction from target STR BAM locus (analyze_target_locus)
- Somatic mosaicism statistics and GMM fitting
- Profiling of alternative motifs and interruptions
- Read structural segmentation (repeat / interruption / flanking)
- Diagnostic plots generation (length+methylation, gene linear map, epigenetic waterfall)
- Clinical reports generation in Markdown/HTML format
"""

import os
import json
import markdown
import numpy as np
import pandas as pd
import pysam
import edlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
import seaborn as sns
from collections import Counter
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde, norm, mannwhitneyu
from scipy.signal import find_peaks


STR_KNOWN_MOTIFS = {
    "DM1": ["CTG", "CAG", "CCG", "CGG", "CTC", "GAG", "CTA", "TAG", "CAA", "TTG", "CGC", "GCG"],
    "DMPK": ["CTG", "CAG", "CCG", "CGG", "CTC", "GAG", "CTA", "TAG", "CAA", "TTG", "CGC", "GCG"],
    "RFC1": ["AAAAG", "AAGGG", "AAAGG", "AAGAG", "AGGGG"],
    "C9orf72": ["GGGGCC", "GGGGGC", "CGGGGC", "GGCCCC"],
    "FXN": ["GAA", "TTC", "GGA", "TCC", "GAC", "GTC"],
    "HTT": ["CAG", "CTG", "CCG", "CGG", "CAA", "TTG"],
    "DEFAULT": ["CTG", "CAG", "CCG", "CGG", "CAA", "TTG", "CTC", "GAG", "TGT", "ACA", "CAC"],
}

STR_REFERENCE_DB = {
    "C9ORF72": {
        "Disease": "ALS / FTD",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "CCCCGG",
        "Strand": "reverse",
        "Normal_Max_Repeats": 23,
        "Pathogenic_Min_Repeats": 30,
        "Coordinates_hg38": "chr9:27573436-27573546",
        "Coordinates_T2T": "chr9:27584065-27584155",
        "Epigenetic_Target": True,
        "Interruption_Motifs": ["CCGACC", "CGGG", "CT", "CTCGGG", "CGG", "TCCCTC"],
    },
    "RFC1": {
        "Disease": "CANVAS",
        "Inheritance": "Autosomal Recessive",
        "Motif_Ref": "GGGAA",
        "Strand": "forward",
        "Normal_Max_Repeats": 200,
        "Pathogenic_Min_Repeats": 400,
        "Coordinates_hg38": "chr4:39348424-39348485",
        "Coordinates_T2T": "chr4:39318077-39318136",
        "Epigenetic_Target": False,
        "Interruption_Motifs": ["GGCT", "GGC", "AGCCTA", "CCTA"],
        "Filt": 15,
    },
    "HTT": {
        "Disease": "Huntington disease",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "CAG",
        "Strand": "forward",
        "Normal_Max_Repeats": 26,
        "Pathogenic_Min_Repeats": 27,
        "Coordinates_hg38": "chr4:3074877-3074933",
        "Coordinates_T2T": "chr4:3073603-3073687",
        "Epigenetic_Target": False,
        "Interruption_Motifs": ["CAG", "CTG", "CCG", "CGG", "CAA", "TTG"],
    },
    "DMPK": {
        "Disease": "Myotonic Dystrophy type 1",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "CAG",
        "Strand": "reverse",
        "Normal_Max_Repeats": 20,
        "Pathogenic_Min_Repeats": 50,
        "Coordinates_hg38": "chr19:45770204-45770264",
        "Coordinates_T2T": "chr19:48597739-48597756",
        "Epigenetic_Target": False,
        "Interruption_Motifs": ["GAG", "CAC", "CCG", "CTC", "TGTG", "CT", "CA"],
    },
    "FGF14": {
        "Disease": "SCA27",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "GAA",
        "Strand": "reverse",
        "Normal_Max_Repeats": 179,
        "Pathogenic_Min_Repeats": 180,
        "Coordinates_hg38": "chr13:102161574-102161726",
        "Coordinates_T2T": "chr13:101377549-101377792",
        "Epigenetic_Target": False
    },
    "BEAN1": {
        "Disease": "SCA31",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "AAAAT",
        "Strand": "forward",
        "Normal_Max_Repeats": 100,
        "Pathogenic_Min_Repeats": 110,
        "Coordinates_hg38": "chr16:66490396-66490466",
        "Coordinates_T2T": "chr16:72284666-72284761",
        "Epigenetic_Target": False
    },
    "DAB1": {
        "Disease": "SCA37",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "AAAAT",
        "Strand": "reverse",
        "Normal_Max_Repeats": 30,
        "Pathogenic_Min_Repeats": 31,
        "Coordinates_hg38": "chr1:57367024-57367124",
        "Coordinates_T2T": "chr1:57245935-57245973",
        "Epigenetic_Target": False
    },
    "AR": {
        "Disease": "Kennedy disease",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "CAG",
        "Strand": "forward",
        "Normal_Max_Repeats": 34,
        "Pathogenic_Min_Repeats": 38,
        "Coordinates_hg38": "chrX:67545316-67545419",
        "Coordinates_T2T": "chrX:65975147-65975250",
        "Epigenetic_Target": False
    },
    "ATXN3": {
        "Disease": "SCA3",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "CTG",
        "Strand": "reverse",
        "Normal_Max_Repeats": 44,
        "Pathogenic_Min_Repeats": 45,
        "Coordinates_hg38": "chr14:92071009-92071060",
        "Coordinates_T2T": "chr14:86300519-86300603",
        "Epigenetic_Target": False
    },
    "BENCHMARK": {
        "Disease": "Benchmark",
        "Inheritance": "Autosomal Dominant",
        "Motif_Ref": "CAG",
        "Strand": "reverse",
        "Normal_Max_Repeats": 20,
        "Pathogenic_Min_Repeats": 50,
        "Coordinates_hg38": "chr19:45770202-45770264",
        "Coordinates_T2T": "",
        "Epigenetic_Target": False,
        "Interruption_Motifs": ["CGG", "CAA", "GAG", "CAT", "CCG", "CG", "CTG", "AG", "TG", "GTG"],
    },
}


def analyze_target_locus(bam_path, fasta_path, chrom, start, end, motif, flank_size=50):
    samfile = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)
    
    left_flank_ref = fasta.fetch(chrom, start - flank_size, start).upper()
    right_flank_ref = fasta.fetch(chrom, end, end + flank_size).upper()
    
    results = []
    motif_len = len(motif)
    max_allowed_errors = int(flank_size * 0.35)
    
    for read in samfile.fetch(chrom, start - 150, end + 150):
        if read.is_unmapped or read.is_secondary or read.is_duplicate:
            continue
            
        read_seq = read.query_sequence
        if not read_seq:
            continue
        read_seq = read_seq.upper()
        
        str_length = None
        motif_purity = None
        seq_extracted = ""
        strategy_used = "Failed"
        
        align_left = edlib.align(left_flank_ref, read_seq, mode="HW", task="locations")
        align_right = edlib.align(right_flank_ref, read_seq, mode="HW", task="locations")
        
        left_valid = bool(align_left["locations"] and align_left["editDistance"] <= max_allowed_errors)
        right_valid = bool(align_right["locations"] and align_right["editDistance"] <= max_allowed_errors)
        
        if left_valid and right_valid:
            str_start_in_read = align_left["locations"][0][1] + 1
            str_end_in_read = align_right["locations"][0][0]
            
            if str_end_in_read > str_start_in_read:
                seq_extracted = read_seq[str_start_in_read:str_end_in_read]
                strategy_used = "Flanking_Anchors"
            else:
                if align_left["editDistance"] <= align_right["editDistance"]:
                    right_valid = False
                else:
                    left_valid = False
        
        if left_valid and not right_valid:
            str_start_in_read = align_left["locations"][0][1] + 1
            seq_extracted = read_seq[str_start_in_read:]
            strategy_used = "Truncated_Left_Only"
            
        elif right_valid and not left_valid:
            str_end_in_read = align_right["locations"][0][0]
            seq_extracted = read_seq[:str_end_in_read]
            strategy_used = "Truncated_Right_Only"
            
        if seq_extracted:
            str_length = len(seq_extracted)
            kmers = [seq_extracted[i:i+motif_len] for i in range(0, len(seq_extracted) - motif_len + 1, motif_len)]
            motif_count = kmers.count(motif)
            motif_purity = (motif_count * motif_len) / str_length if str_length > 0 else 0
        
        if str_length is None:
            pairs = read.get_aligned_pairs(matches_only=False)
            left_q = next((q for q, r in pairs if r is not None and r <= start and q is not None), None)
            right_q = next((q for q, r in reversed(pairs) if r is not None and r >= end and q is not None), None)
            
            if left_q is not None:
                if right_q is not None and right_q > left_q:
                    seq_extracted = read_seq[left_q:right_q]
                    strategy_used = "BAM_Coordinate_Map"
                else:
                    seq_extracted = read_seq[left_q:]
                    strategy_used = "BAM_Coordinate_Map_Truncated"
                
                str_length = len(seq_extracted)
                kmers = [seq_extracted[i:i+motif_len] for i in range(0, len(seq_extracted) - motif_len + 1, motif_len)]
                motif_count = kmers.count(motif)
                motif_purity = (motif_count * motif_len) / str_length if str_length > 0 else 0

        if str_length is None or not seq_extracted:
            continue
            
        locus_mods = {}
        if hasattr(read, "modified_bases_forward") and read.modified_bases_forward:
            ref_positions = read.get_reference_positions(full_length=True)
            
            if ('C', 0, 'm') in read.modified_bases_forward:
                for pos_in_read, qual in read.modified_bases_forward[('C', 0, 'm')]:
                    if pos_in_read < len(ref_positions):
                        g_pos = ref_positions[pos_in_read]
                        if g_pos and start <= g_pos <= end:
                            locus_mods.setdefault(g_pos, {'m': 0.0, 'h': 0.0})['m'] = qual / 255.0
                            
            if ('C', 0, 'h') in read.modified_bases_forward:
                for pos_in_read, qual in read.modified_bases_forward[('C', 0, 'h')]:
                    if pos_in_read < len(ref_positions):
                        g_pos = ref_positions[pos_in_read]
                        if g_pos and start <= g_pos <= end:
                            locus_mods.setdefault(g_pos, {'m': 0.0, 'h': 0.0})['h'] = qual / 255.0

        if locus_mods:
            m_values = [d['m'] for d in locus_mods.values()]
            h_values = [d['h'] for d in locus_mods.values()]
            combined_values = [min(1.0, d['m'] + d['h']) for d in locus_mods.values()]
            
            mean_5mC, mean_5hmC, mean_combined = np.mean(m_values), np.mean(h_values), np.mean(combined_values)
        else:
            mean_5mC, mean_5hmC, mean_combined = 0.0, 0.0, 0.0
            
        results.append({
            "Read_ID": read.query_name,
            "Length_bp": str_length,
            "Motif_Purity": motif_purity,
            "Meth_5mC": mean_5mC,
            "Meth_5hmC": mean_5hmC,
            "Meth_Combined": mean_combined,
            "Strategy": strategy_used,
            "Sequence": seq_extracted
        })
        
    samfile.close()
    fasta.close()
    return pd.DataFrame(results)


def calculate_mosaicism(df):
    lengths = df["Length_bp"].dropna().values
    if len(lengths) < 5:
        return {"Alleles": [np.mean(lengths) if len(lengths) > 0 else 0], "Mosaicism_Index": 0.0, "KDE_X": None, "KDE_Y": None}
        
    kde = gaussian_kde(lengths, bw_method='scott')
    x_eval = np.linspace(min(lengths) - 20, max(lengths) + 50, 500)
    y_eval = kde(x_eval)
    
    peaks, _ = find_peaks(y_eval, distance=15, prominence=0.005)
    alleles = x_eval[peaks]
    mosaicism_index = np.var(lengths)
    
    return {"Alleles": alleles.tolist(), "Mosaicism_Index": mosaicism_index, "KDE_X": x_eval, "KDE_Y": y_eval}


def profile_alternative_motifs(df, primary_motif):
    k = len(primary_motif)
    all_kmer_counts = Counter()
    read_profiles = []
    
    df_valid = df[df["Sequence"].notna()].copy()
    
    for idx, row in df_valid.iterrows():
        seq = row["Sequence"]
        read_id = row["Read_ID"]
        
        kmers = [seq[i:i+k] for i in range(0, len(seq) - k + 1, k)]
        read_kmer_counts = Counter(kmers)
        all_kmer_counts.update(kmers)
        
        total_kmers = sum(read_kmer_counts.values())
        profile = {kmer: (count / total_kmers) * 100 for kmer, count in read_kmer_counts.items()}
        profile["Read_ID"] = read_id
        profile["Total_Length_bp"] = row["Length_bp"]
        read_profiles.append(profile)
        
    df_profiles = pd.DataFrame(read_profiles).fillna(0)
    kmer_columns = [col for col in df_profiles.columns if col not in ["Read_ID", "Total_Length_bp"]]
    sorted_kmers = [item[0] for item in all_kmer_counts.most_common() if item[0] in kmer_columns]
    df_profiles = df_profiles[["Read_ID", "Total_Length_bp"] + sorted_kmers]
    
    return all_kmer_counts, df_profiles


def fit_best_gmm(df, max_components=3):
    lengths = df["Length_bp"].dropna().values.reshape(-1, 1)
    if len(lengths) < 5: return None, None
    
    best_bic = np.inf
    best_gmm = None
    actual_max = min(max_components, len(lengths))
    
    for k in range(1, actual_max + 1):
        gmm = GaussianMixture(n_components=k, random_state=42, covariance_type='full', reg_covar=9.0)
        gmm.fit(lengths)
        bic = gmm.bic(lengths)
        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm
            
    peaks_summary = []
    for i in range(best_gmm.n_components):
        peaks_summary.append({
            "Mean": best_gmm.means_[i][0],
            "Std": np.sqrt(best_gmm.covariances_[i][0][0]),
            "Weight": best_gmm.weights_[i]
        })
    return best_gmm, sorted(peaks_summary, key=lambda x: x["Mean"])


def print_global_cohort_top_motifs(df, primary_motif, interrupt_motif=None):
    k = len(primary_motif)
    repeat_perms = {primary_motif[i:] + primary_motif[:i] for i in range(len(primary_motif))}
    interrupt_perms = {interrupt_motif[i:] + interrupt_motif[:i] for i in range(len(interrupt_motif))} if interrupt_motif else set()
    
    all_kmers_normalized = []
    
    for seq in df["Sequence"].dropna():
        for i in range(0, len(seq) - k + 1, k):
            kmer = seq[i:i+k]
            if kmer in repeat_perms:
                all_kmers_normalized.append(primary_motif)
            elif kmer in interrupt_perms:
                all_kmers_normalized.append(interrupt_motif)
            else:
                all_kmers_normalized.append(kmer)
                
    counts = Counter(all_kmers_normalized)
    total_kmers = sum(counts.values())
    
    print("=== GLOBAL COHORT: TOP 5 MOST FREQUENT MOTIFS (NORMALIZED PHASE) ===")
    for motif, count in counts.most_common(5):
        pct = (count / total_kmers) * 100 if total_kmers > 0 else 0
        if motif == primary_motif:
            label = "MAIN (Grouped & Normalized)"
        elif interrupt_motif and motif == interrupt_motif:
            label = "BIOLOGICAL INTERRUPTION"
        else:
            label = "VARIANT / SEQUENCING NOISE"
            
        print(f"Motif: {motif} | Total count: {count} | Percentage: {pct:.2f}% ({label})")


def detect_most_frequent_interrupt(df, primary_motif, max_error=None):
    k = len(primary_motif)
    if max_error is None:
        max_error = 0 if k <= 3 else 1
        
    repeat_perms = {primary_motif[i:] + primary_motif[:i] for i in range(k)}
    
    all_candidate_kmers = []
    for seq in df["Sequence"].dropna():
        for i in range(0, len(seq) - k + 1, k):
            kmer = seq[i:i+k]
            if len(kmer) == k:
                if is_fuzzy_match(kmer, repeat_perms, max_distance=max_error):
                    continue
                all_candidate_kmers.append(kmer)
                
    if not all_candidate_kmers:
        print("ℹ️ No biological alternative motif detected. Pure locus.")
        return None
        
    counts = Counter(all_candidate_kmers)
    most_common_interrupt, count = counts.most_common(1)[0]
    
    print(f"🎯 Real biological interruption automatically detected: {most_common_interrupt} (n={count})")
    return most_common_interrupt


def complementary_reverse(seq):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(pairs.get(base, base) for base in reversed(seq))


def plot_sequential_locus_structure(
    df,
    repeat_motif,
    path,
    interrupt_motif=None,
    max_error=0.5, #beware this error that could be modified if needed
    strand="+",
):
    """Generates a horizontal broken bar plot showing internal locus structure per read.

    Assigns individual unique colors to each interruption motif found.
    Displays complementary reverse sequence in the legend if strand is '-'.
    """
    complete_STR, read_ids = build_complete_str_structures(
        df, repeat_motif, interrupt_motif, max_error
    )

    num_reads = len(read_ids)
    fig, ax = plt.subplots(figsize=(14, max(num_reads * 0.4, 6)))

    # Process interruption motifs to build dynamic color mappings
    if isinstance(interrupt_motif, (list, tuple)):
        all_int_motifs = [m for m in interrupt_motif if m]
    elif interrupt_motif:
        all_int_motifs = [interrupt_motif]
    else:
        all_int_motifs = []

    palette = [
    "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#56B4E9",
    "#E69F00", "#0072B2", "#999999", "#332288", "#88CCEE",
    "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677",
    "#882255", "#AA4499", "#661100", "#6699CC", "#AA4466",
    "#4477AA", "#228833", "#EE6677", "#BBBBBB"
    ]
    int_color_map = {}
    for idx, m in enumerate(all_int_motifs):
        int_color_map[m] = palette[idx % len(palette)]

    color_map = {
        "Repeat": "#0072B2",  # Standard Blue
        "Other": "#4D4D4D",  # Dark Grey
        "Truncated_Repeat": "#E69F00",  # Dark Orange
        "Truncated_Other": "#999999",  # Light Grey
    }

    for y_idx, segments in enumerate(complete_STR):
        current_pos = 0
        for motif_str, length, seg_type in segments:
            if length == 0:
                continue

            if "Interruption" in seg_type:
                color = int_color_map.get(motif_str, "#009E73")
            else:
                color = color_map.get(seg_type, color_map["Other"])

            ax.broken_barh(
                [(current_pos, length)],
                (y_idx - 0.35, 0.7),
                facecolors=color,
                edgecolor="black",
                linewidth=0.3,
            )
            current_pos += length

    ax.set_yticks(range(num_reads))
    ax.set_yticklabels(read_ids, fontsize=9, fontfamily="monospace")
    ax.set_xlabel(
        "Cumulative Locus Length (bp)",
        fontsize=11,
        fontweight="bold",
        labelpad=10,
    )
    ax.set_ylabel("Read ID", fontsize=11, fontweight="bold", labelpad=10)

    # Controlla se convertire il testo per la legenda
    is_reverse = str(strand).strip() in ["-", "-1", "reverse", "Reverse"]

    display_repeat = (
        complementary_reverse(repeat_motif) if is_reverse else repeat_motif
    )
    ax.set_title(
        f"Sequential Locus Anatomy ({display_repeat})",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    ax.xaxis.grid(True, linestyle="--", alpha=0.6, color="#CCCCCC")
    ax.set_axisbelow(True)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Legenda con termini convertiti se reverse
    legend_elements = [
        Patch(
            facecolor=color_map["Repeat"],
            edgecolor="black",
            linewidth=0.5,
            label=f"Repeat ({display_repeat})",
        ),
    ]

    for m in all_int_motifs:
        display_int = complementary_reverse(m) if is_reverse else m
        legend_elements.append(
            Patch(
                facecolor=int_color_map[m],
                edgecolor="black",
                linewidth=0.5,
                label=f"Interruption ({display_int})",
            )
        )

    legend_elements.append(
        Patch(
            facecolor=color_map["Other"],
            edgecolor="black",
            linewidth=0.5,
            label="Flanking / Other DNA",
        )
    )

    ax.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        frameon=True,
        facecolor="#F9F9F9",
    )

    plt.tight_layout()
    plt.savefig(os.path.join(path, "sequential_locus_structure.png"))
    plt.close()

def plot_gmm_and_dual_methylation(df, gmm_model, peaks_summary, chrom, start, end, path):
    if df.empty:
        print("Insufficient data: Empty DataFrame.")
        return
        
    lengths = df["Length_bp"].dropna().values
    if len(lengths) == 0:
        print("No valid length coordinates found.")
        return

    x_eval = np.linspace(max(0, min(lengths) - 30), max(lengths) + 50, 1000).reshape(-1, 1)
    
    fig, axes = plt.subplots(4, 1, figsize=(12, 11), sharex=True, 
                             gridspec_kw={'height_ratios': [1.2, 1, 1, 1]})
    ax_gmm, ax_5mc, ax_5hmc, ax_comb = axes
    
    if gmm_model is not None:
        log_density = gmm_model.score_samples(x_eval)
        ax_gmm.plot(x_eval, np.exp(log_density), color="black", lw=2, label="GMM Total Profile")
    else:
        if len(lengths) > 1:
            sns.kdeplot(lengths, ax=ax_gmm, color="black", lw=1.5, linestyle="--", label="Empirical KDE (Low Cov)")
        else:
            ax_gmm.axvline(x=lengths[0], color="black", lw=1.5, linestyle="--")

    colors = ["#2E4053", "#27AE60", "#E67E22"]
    max_weight_found = 0.1

    for i, peak in enumerate(peaks_summary or []):
        mean = peak.get("Mean", peak.get("length"))
        std = peak.get("Std", 2.0)
        weight = peak.get("Weight", peak.get("percentage"))
        if weight is None:
            weight = peak.get("count", len(df)) / len(df)
        elif weight > 1.0: 
            weight = weight / 100.0
            
        if mean is None:
            continue
            
        max_weight_found = max(max_weight_found, weight)
        component_density = weight * norm.pdf(x_eval, mean, std)
        correction_factor = std * np.sqrt(2 * np.pi)
        visual_density = component_density * correction_factor
        
        color = colors[i % len(colors)]
        
        if peak.get("type") == "Single_Allele_Fallback_LowCoverage":
            allele_type = "Estimated Allele (Low Coverage)"
            label_text = f"{allele_type}:\nMedian = {mean:.1f} bp\nReads = {peak.get('count', len(df))} ({weight*100:.0f}%)"
        else:
            allele_type = "Wild-Type" if i == 0 else f"Expanded (Sub-clone {i})"
            label_text = f"{allele_type}:\nμ = {mean:.1f} bp\nσ = {std:.1f}\nReads = {weight*100:.1f}%"
        
        ax_gmm.plot(x_eval, visual_density, color=color, lw=2, linestyle="-")
        ax_gmm.fill_between(x_eval.flatten(), visual_density.flatten(), alpha=0.15, color=color)
        ax_gmm.axvline(x=mean, color=color, linestyle=":", alpha=0.7)
        ax_gmm.text(mean + 2, max(0.05, weight * 0.4), label_text, 
                    color=color, fontweight='bold', fontsize=9, 
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    ax_gmm.set_ylabel("Fraction of Supporting Reads")
    ax_gmm.set_ylim(0, max_weight_found * 1.4)
    
    title_suffix = " (Low Coverage - Fallback)" if gmm_model is None else ""
    ax_gmm.set_title(f"Allelic Resolution Profile & Epigenetic Deconvolution (5mC / 5hmC) at {chrom}:{start}-{end}{title_suffix}")

    df_sorted = df.sort_values(by="Length_bp").reset_index(drop=True)
    df_sorted['Read_Index'] = df_sorted.index

    sc_5mc = ax_5mc.scatter(df_sorted['Length_bp'], df_sorted['Read_Index'], 
                            c=df_sorted['Meth_5mC'], cmap='Reds', vmin=0.0, vmax=1.0, s=35, edgecolors='black', lw=0.2)
    ax_5mc.set_ylabel("Read Index")
    fig.colorbar(sc_5mc, ax=ax_5mc, orientation='vertical', aspect=10, pad=0.01).set_label("5mC")
    
    sc_5hmc = ax_5hmc.scatter(df_sorted['Length_bp'], df_sorted['Read_Index'], 
                              c=df_sorted['Meth_5hmC'], cmap='Blues', vmin=0.0, vmax=1.0, s=35, edgecolors='black', lw=0.2)
    ax_5hmc.set_ylabel("Read Index")
    fig.colorbar(sc_5hmc, ax=ax_5hmc, orientation='vertical', aspect=10, pad=0.01).set_label("5hmC")
    
    sc_comb = ax_comb.scatter(df_sorted['Length_bp'], df_sorted['Read_Index'], 
                              c=df_sorted['Meth_Combined'], cmap='viridis', vmin=0.0, vmax=1.0, s=35, edgecolors='black', lw=0.2)
    ax_comb.set_ylabel("Read Index")
    ax_comb.set_xlabel("STR Length (bp)")
    fig.colorbar(sc_comb, ax=ax_comb, orientation='vertical', aspect=10, pad=0.01).set_label("Total Modified")

    for ax in axes:
        ax.grid(axis='x', linestyle=':', alpha=0.5)
        if len(df_sorted) <= 10:
            ax.set_yticks(range(len(df_sorted)))
            
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'gmm_dual_methylation.png'))
    plt.close()

import os
import markdown
from pdf2image import convert_from_bytes
from weasyprint import HTML


import fitz  # PyMuPDF



def _save_html_as_png(html_content, output_png_path, dpi=150):
    """Renderizza l'HTML in PDF tramite WeasyPrint e lo converte in PNG usando PyMuPDF (senza Poppler)."""
    try:
        # 1. Converte l'HTML in PDF in memoria
        pdf_bytes = HTML(string=html_content).write_pdf()

        # 2. Apre il PDF con PyMuPDF da byte
        doc = fitz.open(stream=pdf_bytes, filetype="pdf")
        page = doc.load_page(0)  # Prende la prima pagina

        # 3. Calcola il fattore di scala per il DPI desiderato (default PDF = 72 DPI)
        zoom = dpi / 72
        mat = fitz.Matrix(zoom, zoom)
        pix = page.get_pixmap(matrix=mat)

        # 4. Salva il file PNG
        pix.save(output_png_path)
        doc.close()

        print(f"[OK] Saved PNG Report to: {output_png_path}")

    except Exception as e:
        print(f"[WARNING] Impossible saved PNG report: {e}")


def generate_advanced_diagnostic_report(
    gene_symbol,
    df_reads,
    peaks_summary,
    global_counts,
    sample_id="SAMPLE_TEST_ONT",
    path=".",
):
    ref = STR_REFERENCE_DB[gene_symbol]
    motif_len = len(ref["Motif_Ref"])

    if peaks_summary is None:
        if df_reads is not None and not df_reads.empty:
            median_length = df_reads["Length_bp"].median()
            peaks_summary = [
                {
                    "Mean": median_length,
                    "Std": 2.0,
                    "Weight": 1.0,
                    "percentage": 100.0,
                    "count": len(df_reads),
                    "type": "Single_Allele_Fallback_LowCoverage",
                }
            ]
        else:
            peaks_summary = []

    report_md = (
        f"# FULL DIAGNOSTIC REPORT: {gene_symbol} (Dual-Epigenetics Mode)\n"
    )
    report_md += f"**Sample ID:** {sample_id}\n\n"
    report_md += "### 1. GMM Allelic Profile\n| Allele | Mean (bp) | Repeats | Status |\n| :--- | :---: | :---: | :--- |\n"

    for i, peak in enumerate(peaks_summary):
        rep = round(peak["Mean"] / motif_len, 1)
        status = (
            "⚠️ PATHOGENIC"
            if rep >= ref["Pathogenic_Min_Repeats"]
            else "✅ Normal"
        )
        suffix = (
            " (Low Coverage)"
            if peak.get("type") == "Single_Allele_Fallback_LowCoverage"
            else ""
        )
        report_md += f"| Allele {i+1}{suffix} | {peak['Mean']:.1f} | **{rep}** | {status} |\n"

    m5mic = (
        df_reads["Meth_5mC"].mean()
        if (df_reads is not None and not df_reads.empty)
        else 0.0
    )
    m5hmic = (
        df_reads["Meth_5hmC"].mean()
        if (df_reads is not None and not df_reads.empty)
        else 0.0
    )
    m_comb = (
        df_reads["Meth_Combined"].mean()
        if (df_reads is not None and not df_reads.empty)
        else 0.0
    )

    report_md += f"""
---
### 2. Methylation Characterization Profile (5mC / 5hmC)
Direct base-modification tag analysis yields the following average values at the STR locus:

* **Classical Methylation (5mC):** {m5mic:.2f}
* **Transcriptional Hydroxymethylation (5hmC):** {m5hmic:.2f}
* **Combined Modification Profile (5mC + 5hmC):** **{m_comb:.2f}** (Total CpG site occupancy index)

### Clinical-Epigenetic Considerations:
"""
    if m_comb < 0.2:
        report_md += "> 🟢 **Epigenetically Active Locus:** Low total modification levels. Open chromatin state.\n"
    else:
        report_md += f"> 🟡 **Significant Modification:** Locus exhibits an occupancy rate of {m_comb*100:.1f}%.\n"
        if m5mic > m5hmic * 2:
            report_md += "> 🔴 **5mC Skew (Total Repression):** Predominance of 5mC over 5hmC indicates stable gene silencing typical of pathogenic phenotypes.\n"
        elif m5hmic >= m5mic:
            report_md += "> 🔵 **5hmC Enrichment (Dynamic State):** High 5hmC levels suggest active demethylation or fine transcriptional regulation. Locus is not fully silenced.\n"

    html_filename = f"ADVANCED_REPORT_{sample_id}_{gene_symbol}.html"
    png_filename = f"ADVANCED_REPORT_{sample_id}_{gene_symbol}.png"

    html_output_path = os.path.join(path, html_filename)
    png_output_path = os.path.join(path, png_filename)

    html_body = markdown.markdown(
        report_md, extensions=["tables", "fenced_code"]
    )

    html_template = f"""<!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Advanced Report - {sample_id}</title>
        <style>
            body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; line-height: 1.6; margin: 40px auto; max-width: 900px; color: #2c3e50; background-color: #fbfeff; }}
            h1 {{ color: #2c3e50; border-bottom: 3px solid #34495e; padding-bottom: 12px; margin-bottom: 20px; }}
            h3 {{ color: #2980b9; margin-top: 30px; border-bottom: 1px solid #ecf0f1; padding-bottom: 6px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; background: #ffffff; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }}
            th, td {{ border: 1px solid #e2e8f0; padding: 12px 15px; text-align: left; }}
            th {{ background-color: #f8fafc; font-weight: bold; color: #475569; }}
            tr:nth-child(even) {{ background-color: #f8fafc; }}
            blockquote {{ background: #f8fafc; border-left: 4px solid #3b82f6; margin: 1.5em 0; padding: 15px 20px; border-radius: 0 6px 6px 0; }}
            blockquote p {{ margin: 0; }}
            hr {{ border: 0; border-top: 2px solid #e2e8f0; margin: 40px 0; }}
        </style>
    </head>
    <body>{html_body}</body>
    </html>"""

    # Salvataggio HTML
    with open(html_output_path, "w", encoding="utf-8") as f:
        f.write(html_template)
    print(f"[OK] Saved Advanced HTML Report to: {html_output_path}")

    # Salvataggio PNG
    _save_html_as_png(html_template, png_output_path)

    return report_md


def generate_clinical_report_with_structures(
    gene_symbol,
    df_reads,
    peaks_summary,
    allele_formulas,
    sample_id,
    path,
    database,
):
    ref = database.get(
        gene_symbol,
        {
            "Disease": "Unknown",
            "Inheritance": "N/A",
            "Normal_Max_Repeats": 30,
            "Pathogenic_Min_Repeats": 60,
        },
    )
    motif_len = len(ref.get("Motif_Ref", "CAG"))

    report_md = f"""# MOLECULAR DIAGNOSTIC REPORT: STR EXPANSE
**Sample ID:** {sample_id}  |  **Algorithm:** GMM Deconvolution & Structure Tracker
---
## 1. Genetic Locus Overview
* **Target Gene:** `{gene_symbol}`
* **Associated Condition:** {ref['Disease']}
* **Inheritance:** {ref['Inheritance']}

## 2. High-Resolution Allelic Profiling (GMM + STR Formula)
The combined statistical model with interruption tracking resolved the following alleles:

| Allele | Mean Length (bp) | Std Dev ($\\sigma$) | Estimated Repeat Units | Structural Formula (5' $\\rightarrow$ 3') | Clinical Status |
| :--- | :---: | :---: | :---: | :--- | :--- |
"""

    for i, peak in enumerate(peaks_summary):
        rep_units = round(peak["Mean"] / motif_len, 1)
        formula = allele_formulas.get(i, "N/A (Complex variance)")

        if rep_units >= ref["Pathogenic_Min_Repeats"]:
            status = "⚠️ **PATHOGENIC**"
        elif rep_units > ref["Normal_Max_Repeats"]:
            status = "🟡 **GREY ZONE / PRE-MUTATION**"
        else:
            status = "✅ Normal"

        report_md += f"| **Allele {i+1}** | {peak['Mean']:.1f} | {peak['Std']:.1f} | **{rep_units}** | `{formula}` | {status} |\n"

    m5c = (
        df_reads["Meth_5mC"].mean() if "Meth_5mC" in df_reads.columns else 0
    )
    hm5c = (
        df_reads["Meth_5hmC"].mean() if "Meth_5hmC" in df_reads.columns else 0
    )

    report_md += f"""
## 3. Epigenetic Support Profile (Dorado Native Calling)
* **Average 5mC Level (Silencing):** {m5c:.2f}
* **Average 5hmC Level (Active State):** {hm5c:.2f}
"""

    dma_block = analyze_differential_methylation(df_reads, peaks_summary)
    report_md += dma_block

    wf_file = plot_epigenetic_waterfall(df_reads, gene_symbol, sample_id, path)
    jm_file = plot_epigenetic_joint_map(df_reads, peaks_summary, gene_symbol, sample_id, path)

    report_md += "\n## 5. Integrated Graphical Appendix\n"

    if wf_file:
        wf_name = os.path.basename(wf_file)
        report_md += f'#### A. Single-Read Epigenetic Waterfall\n<img src="{wf_name}" width="100%" style="max-width:850px; margin:15px 0; border-radius:6px; box-shadow: 0 1px 4px rgba(0,0,0,0.15);"><br>\n'
    if jm_file:
        jm_name = os.path.basename(jm_file)
        report_md += f'#### B. Joint Allelic Density Topography\n<img src="{jm_name}" width="100%" style="max-width:650px; margin:15px 0; border-radius:6px; box-shadow: 0 1px 4px rgba(0,0,0,0.15);"><br>\n'

    html_filename = f"CLINICAL_REPORT_{sample_id}_{gene_symbol}.html"
    png_filename = f"CLINICAL_REPORT_{sample_id}_{gene_symbol}.png"

    html_output_path = os.path.join(path, html_filename)
    png_output_path = os.path.join(path, png_filename)

    html_body = markdown.markdown(
        report_md, extensions=["tables", "fenced_code"]
    )

    html_template = f"""<!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Clinical Report - {sample_id} ({gene_symbol})</title>
        <style>
            body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; line-height: 1.6; margin: 40px auto; max-width: 900px; color: #2c3e50; background-color: #fbfeff; }}
            h1 {{ color: #2c3e50; border-bottom: 3px solid #34495e; padding-bottom: 12px; margin-bottom: 20px; }}
            h2 {{ color: #2c3e50; margin-top: 35px; border-bottom: 2px solid #bdc3c7; padding-bottom: 8px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; background: #ffffff; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }}
            th, td {{ border: 1px solid #e2e8f0; padding: 12px 15px; text-align: left; }}
            th {{ background-color: #f8fafc; font-weight: bold; color: #475569; }}
            tr:nth-child(even) {{ background-color: #f8fafc; }}
            code {{ background-color: #f1f5f9; padding: 2px 5px; border-radius: 4px; font-family: monospace; font-size: 0.95em; }}
            hr {{ border: 0; border-top: 2px solid #e2e8f0; margin: 40px 0; }}
            blockquote {{ background: #f8fafc; border-left: 4px solid #3b82f6; margin: 1.5em 0; padding: 15px 20px; border-radius: 0 6px 6px 0; }}
        </style>
    </head>
    <body>{html_body}</body>
    </html>"""

    # Salvataggio HTML
    with open(html_output_path, "w", encoding="utf-8") as f:
        f.write(html_template)
    print(f"[OK] Saved Clinical HTML Report to: {html_output_path}")

    # Salvataggio PNG
    _save_html_as_png(html_template, png_output_path)

    return report_md

def parse_str_structure_string(complete_str_list, read_index, strand, rep_length=3, int_length=3):
    segments = complete_str_list[read_index]
    CSTR = pd.DataFrame(segments, columns=['Motif', 'Length', 'Type'])
    
    while len(CSTR) and CSTR.iloc[0]['Type'] == 'Other':
        CSTR = CSTR.iloc[1:]
    while len(CSTR) and CSTR.iloc[-1]['Type'] == 'Other':
        CSTR = CSTR.iloc[:-1]
    CSTR = CSTR.reset_index(drop=True)
    
    if CSTR.empty:
        return "Unmappable locus", CSTR

    rows_to_drop = []
    for idx in range(len(CSTR)):
        row = CSTR.iloc[idx]
        if row['Type'] in ('Other', 'Truncated_Interruption', 'Truncated_Repeat') and row['Length'] <= 2:
            orphan_seq = row['Motif']
            orphan_len = row['Length']
            
            merged = False
            for j in range(idx - 1, -1, -1):
                if CSTR.iloc[j]['Type'] == 'Interruption':
                    CSTR.at[j, 'Motif'] += orphan_seq
                    CSTR.at[j, 'Length'] += orphan_len
                    rows_to_drop.append(idx)
                    merged = True
                    break
            
            if not merged:
                for j in range(idx + 1, len(CSTR)):
                    if CSTR.iloc[j]['Type'] == 'Interruption':
                        CSTR.at[j, 'Motif'] = orphan_seq + CSTR.iloc[j]['Motif']
                        CSTR.at[j, 'Length'] += orphan_len
                        rows_to_drop.append(idx)
                        break

    if rows_to_drop:
        CSTR = CSTR.drop(rows_to_drop).reset_index(drop=True)

    raw_pieces = []
    for i in range(len(CSTR)):
        raw_motif = CSTR.iloc[i]['Motif']
        length_raw = CSTR.iloc[i]['Length']
        seg_type = CSTR.iloc[i]['Type']
        
        if seg_type in ('Repeat', 'Interruption', 'Truncated_Repeat', 'Truncated_Interruption'):
            unit_len = len(raw_motif) 
            units = int(round(length_raw / unit_len, 0))
            if units > 0:
                raw_pieces.append([raw_motif, units])
        else:
            raw_pieces.append([raw_motif, 1])

    processed_pieces = []
    for motif, count in raw_pieces:
        if (len(motif) < 3 and count <= 2) or count == 1:
            flat_seq = motif * count
            if processed_pieces and isinstance(processed_pieces[-1], str):
                processed_pieces[-1] += flat_seq
            else:
                processed_pieces.append(flat_seq)
        else:
            processed_pieces.append([motif, count])

    final_pieces = []
    for item in processed_pieces:
        if isinstance(item, str):
            if final_pieces and isinstance(final_pieces[-1], str):
                final_pieces[-1] += item
            else:
                final_pieces.append(item)
        else:
            motif, count = item
            if final_pieces and isinstance(final_pieces[-1], list) and final_pieces[-1][0] == motif:
                final_pieces[-1][1] += count
            else:
                final_pieces.append([motif, count])

    if strand == 'reverse':
        converted_pieces = []
        for item in reversed(final_pieces):
            if isinstance(item, str):
                converted_pieces.append(complementary_reverse(item))
            else:
                motif, count = item
                converted_pieces.append([complementary_reverse(motif), count])
        final_pieces = converted_pieces

    cstr_formula = ""
    for item in final_pieces:
        if isinstance(item, str):
            cstr_formula += f"({item})"
        else:
            motif, count = item
            cstr_formula += f"({motif}){count}"
    return cstr_formula, CSTR


def merge_short_other_segments(segments, min_other_len=2):
    if not segments:
        return segments

    cleaned = [list(seg) for seg in segments]

    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(cleaned):
            motif, length, seg_type = cleaned[i]
            is_other = seg_type in ("Other", "Truncated_Other")

            if is_other and length <= max(1, min_other_len) and len(cleaned) > 1:
                if i > 0 and "Interruption" in cleaned[i - 1][2]:
                    cleaned[i - 1][0] += motif
                    cleaned[i - 1][1] += length
                    del cleaned[i]
                    changed = True
                    break
                elif i < len(cleaned) - 1 and "Interruption" in cleaned[i + 1][2]:
                    cleaned[i + 1][0] = motif + cleaned[i + 1][0]
                    cleaned[i + 1][1] += length
                    del cleaned[i]
                    changed = True
                    break
                else:
                    if i > 0:
                        cleaned[i - 1][1] += length
                    else:
                        cleaned[i + 1][1] += length
                    del cleaned[i]
                    changed = True
                    break
            i += 1

        j = 0
        while j < len(cleaned) - 1:
            motif_a, len_a, type_a = cleaned[j]
            motif_b, len_b, type_b = cleaned[j + 1]
            if motif_a == motif_b and type_a == type_b:
                cleaned[j][1] += len_b
                del cleaned[j + 1]
                changed = True
            else:
                j += 1

    return cleaned


def smooth_flickering_blocks(segments, max_block_len):
    if len(segments) < 2:
        return segments

    structural_types = ("Repeat", "Truncated_Repeat", "Interruption", "Truncated_Interruption")
    cleaned = [list(s) for s in segments]

    changed = True
    while changed:
        changed = False

        candidates = []
        for idx, (_, length, seg_type) in enumerate(cleaned):
            if seg_type in structural_types and length <= max_block_len:
                if seg_type in ("Interruption", "Truncated_Interruption"):
                    has_left_repeat = idx > 0 and "Repeat" in cleaned[idx - 1][2]
                    has_right_repeat = idx < len(cleaned) - 1 and "Repeat" in cleaned[idx + 1][2]
                    if has_left_repeat and has_right_repeat:
                        continue
                candidates.append(idx)
                
        if not candidates:
            break

        idx = min(candidates, key=lambda j: cleaned[j][1])
        length = cleaned[idx][1]

        left_len = cleaned[idx - 1][1] if idx > 0 else -1
        right_len = cleaned[idx + 1][1] if idx < len(cleaned) - 1 else -1

        if left_len == -1 and right_len == -1:
            break

        if right_len > left_len:
            cleaned[idx + 1][1] += length
        else:
            cleaned[idx - 1][1] += length
        del cleaned[idx]
        changed = True

        j = 0
        while j < len(cleaned) - 1:
            motif_a, len_a, type_a = cleaned[j]
            motif_b, len_b, type_b = cleaned[j + 1]
            if motif_a == motif_b and type_a == type_b:
                cleaned[j][1] += len_b
                del cleaned[j + 1]
                changed = True
            else:
                j += 1

    return cleaned


def build_complete_str_structures(df, repeat_motif, interrupt_motif=None, max_error=1,
                                   min_other_len=None, max_flicker_block_len=None,
                                   min_confidence=0.6, confidence_decay=0.25,
                                   hysteresis_tolerance=1):
    complete_STR = []
    read_ids = []

    if min_other_len is None:
        min_other_len = len(repeat_motif)

    k_rep = len(repeat_motif)

    if not interrupt_motif:
        interrupt_motifs = []
    elif isinstance(interrupt_motif, str):
        interrupt_motifs = [interrupt_motif]
    else:
        interrupt_motifs = list(interrupt_motif)

    motifs_by_len = {}
    for m in interrupt_motifs:
        l = len(m)
        motifs_by_len.setdefault(l, [])
        if m not in motifs_by_len[l]:
            motifs_by_len[l].append(m)

    interrupt_perms_by_len = {}
    for l, ml in motifs_by_len.items():
        perms = set()
        for m in ml:
            perms |= {m[i:] + m[:i] for i in range(len(m))}
        interrupt_perms_by_len[l] = perms

    if max_flicker_block_len is None:
        max_flicker_block_len = 3 * k_rep

    err_rep = max_error
    err_int = max_error

    repeat_perms = {repeat_motif[i:] + repeat_motif[:i] for i in range(k_rep)}
    
    for _, row in df.iterrows():
        seq = row["Sequence"]
        read_id = row["Read_ID"]
        status = row.get("Status", "Complete") 
        
        if pd.isna(seq) or not seq:
            complete_STR.append([["None", 0, "Other"]])
            read_ids.append(read_id)
            continue
            
        segments = []
        i = 0
        current_type = None
        current_motif = ""
        current_len = 0

        conf_rep = 1.0
        conf_int = 1.0

        repeat_tag = "Repeat" if status == "Complete" else "Truncated_Repeat"
        interrupt_tag = "Interruption" if status == "Complete" else "Truncated_Interruption"
        
        while i < len(seq):
            window_rep = seq[i:i+k_rep]
            dist_rep = (min(levenshtein_dist(window_rep, p) for p in repeat_perms)
                        if len(window_rep) == k_rep else None)
            rep_exact = dist_rep == 0
            rep_fuzzy_possible = dist_rep is not None and dist_rep <= err_rep
            rep_ok = rep_exact or (rep_fuzzy_possible and conf_rep >= min_confidence)

            best_int_dist = None
            best_k_int = 0
            best_int_motif_label = ""
            
            for l, perms in interrupt_perms_by_len.items():
                window_int = seq[i:i+l]
                if len(window_int) == l:
                    d_int = min(levenshtein_dist(window_int, p) for p in perms)
                    if best_int_dist is None or d_int < best_int_dist:
                        best_int_dist = d_int
                        best_k_int = l
                        
                        for original_m in motifs_by_len[l]:
                            orig_perms = {original_m[x:] + original_m[:x] for x in range(len(original_m))}
                            if any(levenshtein_dist(window_int, p) == d_int for p in orig_perms):
                                best_int_motif_label = original_m
                                break

            int_exact = best_int_dist == 0 if best_int_dist is not None else False
            int_fuzzy_possible = best_int_dist is not None and best_int_dist <= err_int
            int_ok = int_exact or (int_fuzzy_possible and conf_int >= min_confidence)

            if rep_ok and int_ok:
                if current_type == repeat_tag:
                    choice = "int" if best_int_dist < (dist_rep - hysteresis_tolerance) else "rep"
                elif current_type == interrupt_tag:
                    choice = "rep" if dist_rep < (best_int_dist - hysteresis_tolerance) else "int"
                else:
                    choice = "int" if best_int_dist < dist_rep else "rep"
            elif int_ok:
                choice = "int"
            elif rep_ok:
                choice = "rep"
            else:
                choice = "other"

            strength_rep = 1.0 if rep_exact else (0.4 if rep_fuzzy_possible else 0.0)
            conf_rep = (1 - confidence_decay) * conf_rep + confidence_decay * strength_rep
            
            if best_k_int > 0:
                strength_int = 1.0 if int_exact else (0.4 if int_fuzzy_possible else 0.0)
                conf_int = (1 - confidence_decay) * conf_int + confidence_decay * strength_int

            if choice == "int":
                if current_type and current_type != interrupt_tag:
                    segments.append([current_motif, current_len, current_type])
                    current_motif = ""
                    current_len = 0
                elif current_type == interrupt_tag and current_motif != best_int_motif_label:
                    segments.append([current_motif, current_len, current_type])
                    current_motif = ""
                    current_len = 0
                    
                current_type = interrupt_tag
                current_motif = best_int_motif_label
                current_len += best_k_int
                i += best_k_int

            elif choice == "rep":
                if current_type and current_type != repeat_tag:
                    segments.append([current_motif, current_len, current_type])
                    current_motif = ""
                    current_len = 0
                current_type = repeat_tag
                current_motif = repeat_motif 
                current_len += k_rep
                i += k_rep

            else:
                other_tag = "Other" if status == "Complete" else "Truncated_Other"
                if current_type and current_type != other_tag:
                    segments.append([current_motif, current_len, current_type])
                    current_motif = ""
                    current_len = 0
                current_type = other_tag
                current_motif += seq[i]
                current_len += 1
                i += 1
                
        if current_len > 0:
            segments.append([current_motif, current_len, current_type])

        segments = merge_short_other_segments(segments, min_other_len=min_other_len)
        segments = smooth_flickering_blocks(segments, max_block_len=max_flicker_block_len)

        complete_STR.append(segments)
        read_ids.append(read_id)
        
    return complete_STR, read_ids


def draw_dna_gene(sections_lengths, section_colors, str_identifier, sample, path, 
                  motifs_list=None, total_length=None, figsize=(12, 3), 
                  ids=None, show_lengths=True):
    """
    Generates a linear diagram for a single read.
    Creates dynamic legend entries for each unique motif present rather than concatenating them.
    """
    reg = str_identifier.split('_')[1] if '_' in str_identifier else str_identifier

    sections_lengths = list(sections_lengths)[::-1]
    section_colors = list(section_colors)[::-1]

    if motifs_list is not None:
        motifs_list = list(motifs_list)[::-1]

    if total_length is None:
        total_length = sum(sections_lengths)

    if ids is None: ids = "read"

    fig, ax = plt.subplots(figsize=figsize)
    current_position = 0

    for length, color in zip(sections_lengths, section_colors):
        ax.plot(
            [current_position, current_position + length], [0, 0],
            color=color, linewidth=12, solid_capstyle="butt"
        )
        if show_lengths and length > 20:
            ax.text(
                current_position + length / 2, 0.05, str(int(length)),
                ha='center', va='bottom', fontsize=11
            )
        current_position += length

    ax.text(total_length / 2, 0.25, f"Total length: {int(total_length)} bp", ha='center', fontsize=12)

    # Dynamic legend mapping for all distinct motifs in this read
    legend_elements = []
    seen_motifs = set()

    if motifs_list:
        for m, color in zip(motifs_list, section_colors):
            if m not in seen_motifs:
                seen_motifs.add(m)
                label_text = complementary_reverse(m) if len(m) <= 10 else "Other"
                legend_elements.append(Patch(facecolor=color, label=label_text))
    else:
        legend_elements = [
            Patch(facecolor='#E69F00', label="Repeat"),
            Patch(facecolor='#0072B2', label="Interruption"),
            Patch(facecolor='#4D4D4D', label="Other")
        ]

    ax.legend(handles=legend_elements, frameon=False, loc="upper right")

    ax.set_xlim(0, total_length)
    ax.set_ylim(-0.3, 0.4)
    ax.set_title(f'Expansion Anatomy in {reg} ({ids})', fontsize=12, fontweight='bold')
    ax.axis("off")

    base_path = os.path.join(path, f"{ids}_{str_identifier}")
    plt.tight_layout()
    plt.savefig(f"{base_path}.png", dpi=300)
    plt.savefig(f"{base_path}.svg", dpi=300)
    plt.close()


def levenshtein_dist(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_dist(s2, s1)
    if len(s2) == 0:
        return len(s1)
    
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def is_fuzzy_match(window, perms, max_distance=1):
    if window in perms:
        return True
    for perm in perms:
        if len(window) == len(perm) and levenshtein_dist(window, perm) <= max_distance:
            return True
    return False


def analyze_differential_methylation(df_reads, peaks_summary):
    if peaks_summary is None or len(peaks_summary) < 2:
        return """
## 4. Differential Methylation Analysis (DMA)
> ℹ️ **Analysis N/A:** Single allele detected or insufficient coverage for inter-allelic comparison.
"""

    df_analysis = df_reads.copy()
    if "Meth_Combined" not in df_analysis.columns:
        df_analysis["Meth_Combined"] = df_analysis.get("Meth_5mC", 0) + df_analysis.get("Meth_5hmC", 0)
        
    means = [peaks_summary[0]["Mean"], peaks_summary[1]["Mean"]]
    df_analysis['Assigned_Allele'] = df_analysis['Length_bp'].apply(lambda x: np.argmin([abs(x - m) for m in means]))
    
    allele_0_meth = df_analysis[df_analysis['Assigned_Allele'] == 0]['Meth_Combined'].dropna()
    allele_1_meth = df_analysis[df_analysis['Assigned_Allele'] == 1]['Meth_Combined'].dropna()
    
    if len(allele_0_meth) < 3 or len(allele_1_meth) < 3:
        return """
## 4. Differential Methylation Analysis (DMA)
> ⚠️ **Insufficient Statistical Power:** Fewer than 3 reads assigned to one or both alleles.
"""

    m0 = allele_0_meth.mean()
    m1 = allele_1_meth.mean()
    delta_meth = m1 - m0
    
    stat, p_value = mannwhitneyu(allele_0_meth, allele_1_meth, alternative='two-sided')
    is_significant = p_value < 0.05
    status_str = "🔴 SIGNIFICANT" if is_significant else "🟢 NOT SIGNIFICANT"
    
    dma_md = f"""
## 4. Differential Methylation Analysis (DMA)
Evaluates whether total epigenetic modifications (`Meth_Combined`) are homogeneously distributed across alleles resolved by GMM.

* **Allele 1 Mean Methylation (Normal/Short):** {m0:.2f} (n = {len(allele_0_meth)} reads)
* **Allele 2 Mean Methylation (Expanded/Long):** {m1:.2f} (n = {len(allele_1_meth)} reads)
* **Absolute Variation (Delta):** {delta_meth:+.2f}
* **Statistical Test:** Mann-Whitney U test (Two-Sided)
* **p-value:** {p_value:.4f} ({status_str})

### Biological Interpretation:
"""
    if is_significant:
        if delta_meth > 0.15:
            dma_md += f"> ⚠️ **Allele-Specific Hypermethylation:** Allele 2 exhibits a significant increase in methylation compared to Allele 1.\n"
        elif delta_meth < -0.15:
            dma_md += f"> ⚠️ **Selective Short-Allele Hypermethylation:** Allele 1 presents significantly higher modification levels.\n"
        else:
            dma_md += f"> 🟡 **Moderate Asymmetry:** Statistically significant difference with reduced Delta ({delta_meth:+.2f}).\n"
    else:
        dma_md += f"> ✅ **Uniform Epigenetic Profile:** No significant variations between alleles.\n"
        
    return dma_md


def plot_epigenetic_waterfall(df_reads, gene_symbol, sample_id, out_dir):
    if df_reads is None or df_reads.empty:
        return None
    
    df_sorted = df_reads.sort_values(by="Length_bp").reset_index(drop=True)
    if "Meth_Combined" not in df_sorted.columns:
        df_sorted["Meth_Combined"] = df_sorted.get("Meth_5mC", 0) + df_sorted.get("Meth_5hmC", 0)
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6), sharey=True, 
                                   gridspec_kw={'width_ratios': [3, 1]})
    
    cmap = plt.colormaps.get_cmap('coolwarm')
    colors = cmap(df_sorted["Meth_Combined"])
    
    ax1.barh(df_sorted.index, df_sorted["Length_bp"], color=colors, edgecolor='none', height=0.85)
    ax1.set_xlabel("Read Length (bp)", fontsize=10, fontweight='bold')
    ax1.set_ylabel("Sequenced Reads (Sorted by length)", fontsize=10, fontweight='bold')
    ax1.grid(axis='x', linestyle='--', alpha=0.4)
    ax1.set_title(f"Locus {gene_symbol} - Epigenetic Waterfall", fontsize=11, pad=12, fontweight='bold')
    
    ax2.scatter(df_sorted["Meth_Combined"], df_sorted.index, color='#2c3e50', alpha=0.6, edgecolors='none', s=35)
    ax2.set_xlabel("Methylation (0-1)", fontsize=10, fontweight='bold')
    ax2.set_xlim(-0.05, 1.05)
    ax2.grid(axis='x', linestyle='--', alpha=0.4)
    ax2.axvline(x=0.2, color='green', linestyle=':', alpha=0.5)
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    cbar = fig.colorbar(sm, ax=ax2, orientation='vertical', pad=0.15)
    cbar.set_label("Meth_Combined (5mC + 5hmC)", fontsize=9, weight='bold')
    
    plt.tight_layout()
    filename = f"EPIGENETIC_WATERFALL_{sample_id}_{gene_symbol}.png"
    plot_path = os.path.join(out_dir, filename)
    plt.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close()
    return filename


def plot_epigenetic_joint_map(df_reads, peaks_summary, gene_symbol, sample_id, out_dir):
    # 1. Controlli base sugli input
    if df_reads is None or df_reads.empty or peaks_summary is None or len(peaks_summary) < 2:
        return None
        
    df_plot = df_reads.copy()
    if "Meth_Combined" not in df_plot.columns:
        df_plot["Meth_Combined"] = df_plot.get("Meth_5mC", 0) + df_plot.get("Meth_5hmC", 0)

    print(df_plot["Meth_Combined"])
    if (df_plot["Length_bp"].nunique() < 2 or df_plot["Meth_Combined"].nunique() < 2):
        print(f"⚠️ Insufficient variability for KDE plot ({gene_symbol})")
        return None
        
    means = [peaks_summary[0]["Mean"], peaks_summary[1]["Mean"]]
    df_plot['Allele'] = df_plot['Length_bp'].apply(lambda x: f"Allele {np.argmin([abs(x - m) for m in means]) + 1}")
    
    # Palette colori coerente
    palette = {"Allele 1": "#0072B2", "Allele 2": "#E69F00"}

    # 3. Jointplot
    g = sns.jointplot(
        data=df_plot, x="Length_bp", y="Meth_Combined", hue="Allele",
        kind="scatter", palette=palette,
        alpha=0.6, s=50, height=7,
        marginal_kws=dict(fill=True, common_norm=False)
    )
    
    # 4. Aggiunta opzionale del KDE
    try:
        sns.kdeplot(
            data=df_plot, x="Length_bp", y="Meth_Combined", hue="Allele", 
            ax=g.ax_joint, alpha=0.3, levels=4, palette=palette,
            warn_singular=False
        )
    except Exception as e:
        print(f"[INFO] Skipping KDE overlay due to data distribution: {e}")

    # 5. Formattazione e Titoli
    g.ax_joint.set_xlabel("Read Length (bp)", fontweight='bold')
    g.ax_joint.set_ylabel("Total Methylation", fontweight='bold')
    g.fig.suptitle(f"Epigenetic Topography: {gene_symbol} (Sample: {sample_id})", y=1.02, fontweight='bold', fontsize=11)
    
    # 6. Salvataggio e Chiusura
    filename = f"EPIGENETIC_JOINT_MAP_{sample_id}_{gene_symbol}.png"
    plot_path = os.path.join(out_dir, filename)
    
    plt.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close('all')  # Chiude tutte le figure per evitare memory leak
    
    return plot_path


import os
import markdown
from weasyprint import HTML
import fitz  # PyMuPDF
from PIL import Image
import numpy as np
from scipy.stats import mannwhitneyu


# Deve corrispondere al background-color impostato nel CSS del template HTML
_PAGE_BG_COLOR = (251, 254, 255)  # #fbfeff


def _trim_trailing_whitespace(img, bg_color=_PAGE_BG_COLOR, tolerance=10, padding=30):

    arr = np.array(img.convert("RGB")).astype(int)
    diff = np.abs(arr - np.array(bg_color)).sum(axis=2)
    content_rows = np.where((diff > tolerance).any(axis=1))[0]

    if content_rows.size == 0:
        return img

    last_content_row = int(content_rows.max())
    crop_height = min(img.height, last_content_row + padding)
    return img.crop((0, 0, img.width, crop_height))


def _save_html_as_png(html_content, output_png_path, dpi=150):

    try:
        pdf_bytes = HTML(string=html_content).write_pdf()
        doc = fitz.open(stream=pdf_bytes, filetype="pdf")

        zoom = dpi / 72
        mat = fitz.Matrix(zoom, zoom)

        page_images = []
        for page_index in range(doc.page_count):
            page = doc.load_page(page_index)
            pix = page.get_pixmap(matrix=mat)
            img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
            img = _trim_trailing_whitespace(img)
            page_images.append(img)
        doc.close()

        if not page_images:
            print("[WARNING] Nessuna pagina da renderizzare in PNG.")
            return

        if len(page_images) == 1:
            final_img = page_images[0]
        else:
            total_width = max(img.width for img in page_images)
            total_height = sum(img.height for img in page_images)
            final_img = Image.new("RGB", (total_width, total_height), "white")
            y_offset = 0
            for img in page_images:
                final_img.paste(img, (0, y_offset))
                y_offset += img.height

        final_img.save(output_png_path)
        print(
            f"[OK] Saved PNG Report ({len(page_images)} pagina/e, nulla tagliato) "
            f"to: {output_png_path}"
        )

    except Exception as e:
        print(f"[WARNING] Impossible saved PNG report: {e}")


def complementary_reverse(seq):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(pairs.get(base, base) for base in reversed(seq))


def _save_html_as_pdf(html_content, output_pdf_path):
    try:
        HTML(string=html_content).write_pdf(output_pdf_path)
        print(f"[OK] Saved PDF Report to: {output_pdf_path}")
    except Exception as e:
        print(f"[WARNING] Impossible saved PDF report: {e}")


def _assign_reads_to_nearest_allele(df_reads, peaks_summary):

    if df_reads is None or df_reads.empty or not peaks_summary or "Length_bp" not in df_reads.columns:
        return None

    means = [peak["Mean"] for peak in peaks_summary]
    df_plot = df_reads.copy()
    df_plot["Allele_Assigned"] = df_plot["Length_bp"].apply(
        lambda x: int(np.argmin([abs(x - m) for m in means])) + 1
    )
    return df_plot


def analyze_differential_methylation(df_reads, peaks_summary):

    section_md = "\n---\n\n## 4. Differential Methylation Analysis (DMA)\n"
    section_md += (
        "Evaluates whether total epigenetic modifications (`Meth_Combined`) are "
        "homogeneously distributed across alleles resolved by GMM.\n\n"
    )

    if not peaks_summary or len(peaks_summary) < 2:
        section_md += "> ℹ️ **Not applicable:** DMA requires at least 2 resolved alleles.\n"
        return section_md

    if df_reads is None or df_reads.empty or "Meth_Combined" not in df_reads.columns:
        section_md += "> ℹ️ **Not applicable:** `Meth_Combined` read-level data not available.\n"
        return section_md

    # Confronta l'allele piu' corto (Normal/Short) con il piu' lungo (Expanded/Long)
    sorted_peaks = sorted(peaks_summary, key=lambda p: p["Mean"])
    short_peak, long_peak = sorted_peaks[0], sorted_peaks[-1]
    means = [short_peak["Mean"], long_peak["Mean"]]

    df_plot = df_reads.copy()
    df_plot["Allele"] = df_plot["Length_bp"].apply(
        lambda x: f"Allele {np.argmin([abs(x - m) for m in means]) + 1}"
    )

    group1 = df_plot.loc[df_plot["Allele"] == "Allele 1", "Meth_Combined"].dropna()
    group2 = df_plot.loc[df_plot["Allele"] == "Allele 2", "Meth_Combined"].dropna()

    if len(group1) == 0 or len(group2) == 0:
        section_md += "> ℹ️ **Not applicable:** at least one allele group has no assigned reads.\n"
        return section_md

    mean1, mean2 = float(group1.mean()), float(group2.mean())
    delta = mean2 - mean1

    try:
        _, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
    except ValueError:
        section_md += (
            "> ℹ️ **Not applicable:** statistical test could not be computed "
            "(identical distributions or insufficient data).\n"
        )
        return section_md

    significant = p_value < 0.05
    sig_icon = "🔴" if significant else "🟢"
    sig_label = "SIGNIFICANT" if significant else "NOT SIGNIFICANT"

    section_md += f"* **Allele 1 Mean Methylation (Normal/Short):** {mean1:.2f} (n = {len(group1)} reads)\n"
    section_md += f"* **Allele 2 Mean Methylation (Expanded/Long):** {mean2:.2f} (n = {len(group2)} reads)\n"
    section_md += f"* **Absolute Variation (Delta):** {delta:+.2f}\n"
    section_md += "* **Statistical Test:** Mann-Whitney U test (Two-Sided)\n"
    section_md += f"* **p-value:** {p_value:.4f} ({sig_icon} {sig_label})\n\n"

    section_md += "### Biological Interpretation:\n\n"
    if significant:
        if delta > 0:
            section_md += (
                "> ⚠️ **Differential Epigenetic Silencing:** The expanded/long allele shows "
                "significantly higher methylation than the normal/short allele, consistent "
                "with allele-specific silencing typical of repeat expansion disorders.\n"
            )
        else:
            section_md += (
                "> ⚠️ **Differential Epigenetic Silencing:** The normal/short allele shows "
                "significantly higher methylation than the expanded/long allele.\n"
            )
    else:
        section_md += (
            "> ✅ **Uniform Epigenetic Profile:** No significant variations between alleles.\n"
        )

    return section_md


def generate_full_diagnostic_report(
    gene_symbol,
    df_reads,
    peaks_summary,
    allele_formulas=None,
    sample_id="SAMPLE_TEST_ONT",
    path=".",
    database=None,
    page_width_px=1150,
    page_height_px=6000,
):

    allele_formulas = allele_formulas or {}

    ref_db = database if database is not None else STR_REFERENCE_DB  # noqa: F821
    ref = ref_db.get(
        gene_symbol,
        {
            "Disease": "Unknown",
            "Inheritance": "N/A",
            "Motif_Ref": "CAG",
            "Normal_Max_Repeats": 30,
            "Pathogenic_Min_Repeats": 60,
        },
    )
    motif_len = len(ref.get("Motif_Ref", "CAG"))

    # Fallback se peaks_summary manca / bassa coverage (logica presa dal report avanzato)
    if not peaks_summary:
        if df_reads is not None and not df_reads.empty:
            median_length = df_reads["Length_bp"].median()
            peaks_summary = [
                {
                    "Mean": median_length,
                    "Std": 2.0,
                    "Weight": 1.0,
                    "percentage": 100.0,
                    "count": len(df_reads),
                    "type": "Single_Allele_Fallback_LowCoverage",
                }
            ]
        else:
            peaks_summary = []

    # ---------- 1. Genetic Locus Overview ----------
    report_md = f"# FULL NANOEXPANSION REPORT: {gene_symbol}\n"
    report_md += (
        f"**Sample ID:** {sample_id}  |  **Algorithm:** GMM Deconvolution & "
        f"Structure Tracker + Epigenetic Profiling\n\n"
    )
    report_md += "---\n\n"
    report_md += "## 1. Genetic Locus Overview\n"
    report_md += f"* **Target Gene:** `{gene_symbol}`\n"
    report_md += f"* **Associated Condition:** {ref.get('Disease', 'Unknown')}\n"
    report_md += f"* **Inheritance:** {ref.get('Inheritance', 'N/A')}\n"
    strand = ref.get("Strand", "forward")
    motif_ref = ref.get("Motif_Ref", "N/A")
    displayed_motif = (
        complementary_reverse(motif_ref) if strand == "reverse" and motif_ref != "N/A" else motif_ref
    )
    #strand_note = "" if strand == "reverse" else ""
    report_md += f"* **Reference Motif:** `{displayed_motif}` ({motif_len} bp)\n"
    report_md += f"* **Normal Max Repeats:** {ref.get('Normal_Max_Repeats', 'N/A')}\n"
    report_md += f"* **Pathogenic Min Repeats:** {ref.get('Pathogenic_Min_Repeats', 'N/A')}\n\n"

    # ---------- 2. Allele Profile (GMM) ----------
    report_md += "## 2. High-Resolution Allelic Profiling (GMM + STR Formula)\n"
    report_md += (
        "| Allele | Mean (bp) | Std Dev ($\\sigma$) | Repeats | "
        "Structural Formula (5' $\\rightarrow$ 3') | Coverage | Status |\n"
        "| :--- | :---: | :---: | :---: | :--- | :---: | :--- |\n"
    )

    total_reads = len(df_reads) if (df_reads is not None and not df_reads.empty) else None
    df_assigned = _assign_reads_to_nearest_allele(df_reads, peaks_summary)

    for i, peak in enumerate(peaks_summary):
        rep_units = round(peak["Mean"] / motif_len, 1)
        formula = allele_formulas.get(i, "N/A (Complex variance)")
        suffix = (
            " (Low Coverage)"
            if peak.get("type") == "Single_Allele_Fallback_LowCoverage"
            else ""
        )

        if df_assigned is not None:
            n_reads = int((df_assigned["Allele_Assigned"] == i + 1).sum())
            pct = (n_reads / total_reads * 100) if total_reads else None
        else:
            n_reads = peak.get("count", "N/A")
            pct = peak.get("percentage")
            if pct is None and total_reads and isinstance(n_reads, (int, float)):
                pct = (n_reads / total_reads) * 100
        coverage_str = f"{n_reads} reads" + (f" ({pct:.1f}%)" if pct is not None else "")

        if rep_units >= ref.get("Pathogenic_Min_Repeats", 60):
            status = "⚠️ **PATHOGENIC**"
        elif rep_units > ref.get("Normal_Max_Repeats", 30):
            status = "🟡 **GREY ZONE / PRE-MUTATION**"
        else:
            status = "✅ Normal"

        report_md += (
            f"| **Allele {i+1}{suffix}** | {peak['Mean']:.1f} | {peak['Std']:.1f} | "
            f"**{rep_units}** | `{formula}` | {coverage_str} | {status} |\n"
        )

    # ---------- 3. Methylation----------
    has_reads = df_reads is not None and not df_reads.empty
    m5mic = df_reads["Meth_5mC"].mean() if (has_reads and "Meth_5mC" in df_reads.columns) else 0.0
    m5hmic = df_reads["Meth_5hmC"].mean() if (has_reads and "Meth_5hmC" in df_reads.columns) else 0.0
    m_comb = (
        df_reads["Meth_Combined"].mean()
        if (has_reads and "Meth_Combined" in df_reads.columns)
        else (m5mic + m5hmic)
    )

    #report_md += "\n---\n\n### 3. Methylation Characterization Profile (5mC / 5hmC)\n"
    report_md += "\n---\n\n## 3. Methylation Characterization Profile (5mC / 5hmC)\n"
    report_md += "Direct base-modification tag analysis yields the following average values at the STR locus:\n\n"
    report_md += f"* **Classical Methylation (5mC):** {m5mic:.2f}\n"
    report_md += f"* **Transcriptional Hydroxymethylation (5hmC):** {m5hmic:.2f}\n"
    report_md += f"* **Combined Modification Profile (5mC + 5hmC):** **{m_comb:.2f}** (Total CpG site occupancy index)\n\n"
    report_md += "### Clinical-Epigenetic Considerations:\n\n"

    if m_comb < 0.2:
        report_md += "> 🟢 **Epigenetically Active Locus:** Low total modification levels. Open chromatin state.\n"
    else:
        report_md += f"> 🟡 **Significant Modification:** Locus exhibits an occupancy rate of {m_comb*100:.1f}%.\n"
        if m5mic > m5hmic * 2:
            report_md += (
                "> 🔴 **5mC Skew (Total Repression):** Predominance of 5mC over 5hmC "
                "indicates stable gene silencing typical of pathogenic phenotypes.\n"
            )
        elif m5hmic >= m5mic:
            report_md += (
                "> 🔵 **5hmC Enrichment (Dynamic State):** High 5hmC levels suggest active "
                "demethylation or fine transcriptional regulation. Locus is not fully silenced.\n"
            )

    # ---------- 4. Differential Methylation Analysis (DMA) ----------
    report_md += analyze_differential_methylation(df_reads, peaks_summary)

    # ---------- Rendering HTML ----------
    html_body = markdown.markdown(report_md, extensions=["tables", "fenced_code"])
    bg_hex = "#{:02x}{:02x}{:02x}".format(*_PAGE_BG_COLOR)

    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Full Diagnostic Report - {sample_id} ({gene_symbol})</title>
    <style>
        /* Pagina larga e molto alta: evita sia il taglio orizzontale delle
           tabelle sia lo split forzato su piu' pagine A4 (causa dei blocchi
           bianchi tra una sezione e l'altra). L'eccesso di altezza viene
           ritagliato in automatico dopo il rendering (vedi _trim_trailing_whitespace). */
        @page {{
            size: {page_width_px}px {page_height_px}px;
            margin: 0;
        }}
        body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; line-height: 1.6; margin: 40px auto; max-width: {page_width_px - 150}px; color: #2c3e50; background-color: {bg_hex}; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #34495e; padding-bottom: 12px; margin-bottom: 20px; }}
        h2 {{ color: #2c3e50; margin-top: 35px; border-bottom: 2px solid #bdc3c7; padding-bottom: 8px; }}
        h3 {{ color: #2980b9; margin-top: 30px; border-bottom: 1px solid #ecf0f1; padding-bottom: 6px; }}
        table {{ table-layout: fixed; border-collapse: collapse; width: 100%; margin: 20px 0; background: #ffffff; box-shadow: 0 1px 3px rgba(0,0,0,0.05); font-size: 0.92em; }}
        th, td {{ border: 1px solid #e2e8f0; padding: 12px 15px; text-align: left; word-break: break-word; overflow-wrap: break-word; }}
        th {{ background-color: #f8fafc; font-weight: bold; color: #475569; }}
        tr:nth-child(even) {{ background-color: #f8fafc; }}
        tr {{ page-break-inside: avoid; }}
        code {{ background-color: #f1f5f9; padding: 2px 5px; border-radius: 4px; font-family: monospace; font-size: 0.95em; word-break: break-all; white-space: normal; }}
        blockquote {{ background: #f8fafc; border-left: 4px solid #3b82f6; margin: 1.5em 0; padding: 15px 20px; border-radius: 0 6px 6px 0; page-break-inside: avoid; }}
        blockquote p {{ margin: 0; }}
        hr {{ border: 0; border-top: 2px solid #e2e8f0; margin: 40px 0; }}
    </style>
</head>
<body>{html_body}</body>
</html>"""

    html_filename = f"FULL_REPORT_{sample_id}_{gene_symbol}.html"
    png_filename = f"FULL_REPORT_{sample_id}_{gene_symbol}.png"
    pdf_filename = f"FULL_REPORT_{sample_id}_{gene_symbol}.pdf"

    html_output_path = os.path.join(path, html_filename)
    png_output_path = os.path.join(path, png_filename)
    pdf_output_path = os.path.join(path, pdf_filename)

    with open(html_output_path, "w", encoding="utf-8") as f:
        f.write(html_template)
    print(f"[OK] Saved Full HTML Report to: {html_output_path}")

    _save_html_as_png(html_template, png_output_path)
    _save_html_as_pdf(html_template, pdf_output_path)

    return report_md
