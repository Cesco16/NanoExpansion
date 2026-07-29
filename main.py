#!/usr/bin/env python3
"""
main.py
=======
NanoExpansion (NanoExp3) - Command-line pipeline for characterization
of short tandem repeat (STR) expansions from nanopore (long-read) data.

Given an aligned BAM file and a FASTA reference, the pipeline:
  1. extracts, for each read covering the target locus, the length and
     motif composition of the repeated region (with robust recovery
     for reads truncated upstream/downstream of the anchor);
  2. estimates the number and size of alleles via Gaussian Mixture
     Model (GMM) fitting on observed lengths and computes the somatic mosaicism index;
  3. profiles alternative motifs/k-mers (interruptions) present at the locus;
  4. structurally segments each read into Repeat / Interruption / Other blocks,
     producing a synthetic clinical formula (e.g., (CAG)45(CAA)2);
  5. generates diagnostic plots (length distribution, linear allele map,
     methylation if available) and a clinical report in Markdown format,
     along with a summary text file for the entire read cohort.

Typical Usage:
    python main.py --bam sample.bam --fasta hg38.fa --gene DMPK \
        --sample 13204 --outdir output/13204/nanoexpansion/

For the full list of supported genes out of the box, see
utils.STR_REFERENCE_DB (extendable by adding new entries to the dictionary).
"""

import argparse
import os
import sys

from utils import (
    STR_REFERENCE_DB,
    analyze_target_locus,
    calculate_mosaicism,
    profile_alternative_motifs,
    fit_best_gmm,
    print_global_cohort_top_motifs,
    detect_most_frequent_interrupt,
    plot_sequential_locus_structure,
    plot_gmm_and_dual_methylation,
    generate_advanced_diagnostic_report,
    parse_str_structure_string,
    generate_clinical_report_with_structures,
    build_complete_str_structures,
    draw_dna_gene,
    generate_full_diagnostic_report
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="NanoExpansion: Characterization of STR expansions from nanopore data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", required=True, help="Aligned and indexed BAM file (.bai present).")
    parser.add_argument("--fasta", required=True, help="Reference genome in FASTA format (indexed, .fai present).")
    parser.add_argument(
        "--gene",
        required=True,
        help=f"Target gene/locus. Supported genes in internal DB: {', '.join(STR_REFERENCE_DB.keys())}. "
             "If the gene is not listed, use --chrom/--start/--end/--motif to define a custom locus.",
    )
    parser.add_argument("--sample", required=True, help="Sample identifier (used for output naming).")
    parser.add_argument("--outdir", required=True, help="Output directory for plots, reports, and text files.")

    # Manual locus override (useful for genes not present in STR_REFERENCE_DB or non-hg38 assemblies)
    parser.add_argument("--chrom", default=None, help="Chromosome of the locus (manual override).")
    parser.add_argument("--start", type=int, default=None, help="Start coordinate of the locus, 0-based (manual override).")
    parser.add_argument("--end", type=int, default=None, help="End coordinate of the locus (manual override).")
    parser.add_argument("--motif", default=None, help="Reference repeated motif, e.g., CAG (manual override).")
    parser.add_argument("--strand", default=None, choices=["forward", "reverse"], help="Strand of the locus (manual override).")

    # Optional clinical metadata for custom locus
    parser.add_argument("--disease", default=None, help="Associated disease for custom locus (optional).")
    parser.add_argument("--inheritance", default=None, help="Inheritance mode for custom locus (optional).")
    parser.add_argument("--normal-max-repeats", type=float, default=None, help="Maximum normal repeat threshold (optional).")
    parser.add_argument("--pathogenic-min-repeats", type=float, default=None, help="Minimum pathogenic repeat threshold (optional).")
    parser.add_argument("--epigenetic-target", action="store_true", help="Flag indicating if the locus is an epigenetic target.")

    # Structural segmentation algorithm parameters
    parser.add_argument(
        "--interrupt-motifs",
        nargs="*",
        default=None,
        help="Space-separated list of expected interruption motifs. "
             "If omitted, uses 'Interruption_Motifs' from STR_REFERENCE_DB if present.",
    )
    parser.add_argument("--filt", type=int, default=20, help="Maximum length of 'flicker' blocks to remove during smoothing.")
    parser.add_argument("--min-other-len", type=int, default=50, help="Minimum length of 'Other' segments to prevent merging.")
    parser.add_argument("--hysteresis", type=float, default=0.5, help="Hysteresis tolerance for structural read segmentation.")
    parser.add_argument("--max-gmm-components", type=int, default=3, help="Maximum GMM components (alleles) to test.")

    return parser.parse_args()


def resolve_locus(args):
    """
    Determines locus chromosome, start, end, motif, and strand from STR_REFERENCE_DB or CLI overrides.
    """
    gene_key = args.gene.upper()
    gene_info = dict(STR_REFERENCE_DB.get(gene_key, {}))
    is_known_gene = bool(gene_info)

    if not is_known_gene and (args.chrom is None or args.start is None or args.end is None or args.motif is None):
        raise KeyError(
            f"Gene '{args.gene}' is not in internal STR_REFERENCE_DB and missing manual coordinates "
            "(--chrom/--start/--end/--motif)."
        )

    if args.chrom and args.start is not None and args.end is not None and args.motif:
        chrom, start, end, motif = args.chrom, args.start, args.end, args.motif
    else:
        fasta_lower = os.path.basename(args.fasta).lower()
        if "hs1" in fasta_lower or "t2t" in fasta_lower:
            raw_coords = gene_info.get("Coordinates_T2T") or gene_info.get("Coordinates_hg38")
        else:
            raw_coords = gene_info.get("Coordinates_hg38") or gene_info.get("Coordinates_T2T")

        if not raw_coords:
            raise KeyError(f"No coordinates available for gene '{args.gene}' in requested assembly.")

        chrom, positions = raw_coords.split(":")
        start_str, end_str = positions.split("-")
        start, end = int(start_str), int(end_str)
        motif = gene_info.get("Motif_Ref")

    strand = args.strand or gene_info.get("Strand", "forward")

    if not is_known_gene:
        if args.normal_max_repeats is None or args.pathogenic_min_repeats is None:
            print(
                f"[WARNING] Gene '{gene_key}' is not in internal DB: no validated clinical threshold available. "
                "Defaulting to placeholders (Normal_Max_Repeats=30, Pathogenic_Min_Repeats=60).",
                file=sys.stderr,
            )

        new_entry = {
            "Disease": args.disease or "N/A (custom locus)",
            "Inheritance": args.inheritance or "N/A",
            "Motif_Ref": motif,
            "Strand": strand,
            "Normal_Max_Repeats": args.normal_max_repeats if args.normal_max_repeats is not None else 30,
            "Pathogenic_Min_Repeats": args.pathogenic_min_repeats if args.pathogenic_min_repeats is not None else 60,
            "Coordinates_hg38": f"{chrom}:{start}-{end}",
            "Coordinates_T2T": "",
            "Epigenetic_Target": args.epigenetic_target,
            "Interruption_Motifs": args.interrupt_motifs if args.interrupt_motifs else None,
        }
        STR_REFERENCE_DB[gene_key] = new_entry
        gene_info = dict(new_entry)
        print(f"[INFO] Registered new entry for '{gene_key}' in STR_REFERENCE_DB: {new_entry}")

    gene_info.setdefault("Strand", strand)
    gene_info.setdefault("Disease", gene_info.get("Disease", "N/A"))
    gene_info.setdefault("Motif_Ref", motif)

    return chrom, start, end, motif, strand, gene_key, gene_info


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    chrom, start, end, motif, strand, gene_key, gene_info = resolve_locus(args)

    print(f"=== LOADED TARGET LOCUS: {gene_key} ===")
    print(f"Disease: {gene_info.get('Disease', 'N/A')}")
    print(f"CHROM  = '{chrom}'")
    print(f"START  = {start}")
    print(f"END    = {end}")
    print(f"MOTIF  = '{motif}'")
    print(f"STRAND = '{strand}'")
    print("=========================================\n")

    _GENERIC_INTERRUPT_MOTIFS_FALLBACK = ["CGG", "CAA", "GAG", "CAT", "CCG", "CG", "CTG", "AG", "TG", "GTG"]
    if args.interrupt_motifs is not None:
        interrupt_motif = args.interrupt_motifs
        print(f"Interruption motifs: {interrupt_motif or '(auto-detection)'} [source: CLI]")
    elif gene_info.get("Interruption_Motifs"):
        interrupt_motif = list(gene_info["Interruption_Motifs"])
        print(f"Interruption motifs: {interrupt_motif} [source: STR_REFERENCE_DB['{gene_key}']]")
    else:
        interrupt_motif = list(_GENERIC_INTERRUPT_MOTIFS_FALLBACK)
        print(f"Interruption motifs: {interrupt_motif} [source: generic fallback]")

    # --- 1. Read Extraction ------------------------------------------------
    df_results = analyze_target_locus(args.bam, args.fasta, chrom, start, end, motif)
    if df_results is None or df_results.empty:
        print("No valid reads found on target locus. Exiting.", file=sys.stderr)
        sys.exit(1)

    # --- 2. Mosaicism & GMM Fitting ----------------------------------------
    stats = calculate_mosaicism(df_results)
    print("--- ANALYSIS RESULTS ---")
    print(f"Total reads analyzed at locus: {len(df_results)}")
    print(f"Detected alleles (bp): {[round(a, 1) for a in stats['Alleles']]}")
    print(f"Somatic Mosaicism Index: {stats['Mosaicism_Index']:.2f}")

    gmm_model, peaks_summary = fit_best_gmm(df_results, max_components=args.max_gmm_components)

    # --- 3. Alternative Motif Profiling ------------------------------------
    global_counts, df_read_profiles = profile_alternative_motifs(df_results, motif)
    print_global_cohort_top_motifs(df_results, motif, interrupt_motif=None)
    if not interrupt_motif:
        interrupt_motif = detect_most_frequent_interrupt(df_results, motif)

    print("=== GLOBAL COHORT: TOP 5 MOST FREQUENT MOTIFS AT LOCUS ===")
    total_kmers_global = sum(global_counts.values())
    for motif_found, count in global_counts.most_common(5):
        percentage = (count / total_kmers_global) * 100 if total_kmers_global else 0
        status = "MAIN" if motif_found == motif else "ALTERNATIVE/INTERRUPT"
        print(f"Motif: {motif_found} | Total count: {count} | Percentage: {percentage:.2f}% ({status})")

    # --- 4. Length/Methylation Summary Plot --------------------------------
    plot_gmm_and_dual_methylation(df_results, gmm_model, peaks_summary, chrom, start, end, args.outdir)

    # --- 5. Markdown Advanced Diagnostic Report ----------------------------
    advanced_report = generate_advanced_diagnostic_report(
        gene_key, df_results, peaks_summary, global_counts, sample_id=args.sample, path=args.outdir
    )
    print(advanced_report)

    # --- 6. Structural Read Segmentation ----------------------------------
    rep_length = len(motif)

    if isinstance(interrupt_motif, (list, tuple)):
        candidate_lengths = [len(m) for m in interrupt_motif if m]
        int_length = max(candidate_lengths) if candidate_lengths else 3
    elif interrupt_motif:
        int_length = len(interrupt_motif)
    else:
        int_length = 3

    print("Parsing and anatomical read segmentation in progress...")
    complete_STR, read_ids = build_complete_str_structures(
        df_results,
        motif,
        interrupt_motif=interrupt_motif,
        max_flicker_block_len=args.filt,
        min_other_len=args.min_other_len,
        hysteresis_tolerance=args.hysteresis,
    )

    if peaks_summary is None:
        if not df_results.empty:
            median_length = df_results["Length_bp"].median()
            peaks_summary = [{
                "Mean": median_length,
                "Std": 2.0,
                "Weight": 1.0,
                "percentage": 100.0,
                "count": len(df_results),
                "type": "Single_Allele_Fallback_LowCoverage",
            }]
        else:
            peaks_summary = []

    # --- 7. Clinical Formula and Linear Map for each Allele ---------------
    allele_formulas_for_report = {}
    
    # Palette definition for draw_dna_gene
    palette = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]
    all_int_motifs = interrupt_motif if isinstance(interrupt_motif, (list, tuple)) else ([interrupt_motif] if interrupt_motif else [])

    for idx, peak in enumerate(peaks_summary):
        closest_read_idx = (df_results["Length_bp"] - peak["Mean"]).abs().idxmin()
        read_name = df_results.loc[closest_read_idx, "Read_ID"]

        formula_str, clean_df = parse_str_structure_string(
            complete_STR, closest_read_idx, rep_length=rep_length, int_length=int_length, strand=strand
        )
        allele_formulas_for_report[idx] = formula_str

        # Dynamic color mapping for draw_dna_gene
        section_colors = []
        for _, row in clean_df.iterrows():
            stype = row["Type"]
            smotif = row["Motif"]
            if "Repeat" in stype:
                section_colors.append("#E69F00")  # Orange: Repeat
            elif "Interruption" in stype:
                if smotif in all_int_motifs:
                    m_idx = all_int_motifs.index(smotif)
                    section_colors.append(palette[m_idx % len(palette)])
                else:
                    section_colors.append("#0072B2")
            else:
                section_colors.append("#4D4D4D")  # Grey: Flanking/Other

        print(f"--> Generating linear map for Allele {idx + 1} (Read: {read_name})")
        draw_dna_gene(
            sections_lengths=clean_df["Length"],
            section_colors=section_colors,
            str_identifier=gene_key,
            sample=args.sample,
            path=args.outdir,
            motifs_list=clean_df["Motif"].tolist(),
            ids=f"Allele{idx + 1}_{read_name}",
        )

    # --- 7. Clinical Formula and Linear Map for EACH READ -----------------
    read_formulas_for_report = {}

    # Palette definition for draw_dna_gene
    palette = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9"]
    all_int_motifs = (
        interrupt_motif
        if isinstance(interrupt_motif, (list, tuple))
        else ([interrupt_motif] if interrupt_motif else [])
    )

    # Itera direttamente su ciascuna read presente nel DataFrame
    for read_idx, row_read in df_results.iterrows():
        read_name = row_read["Read_ID"]

        # Calcola la formula e pulisce la struttura per la read corrente (usando read_idx)
        formula_str, clean_df = parse_str_structure_string(
            complete_STR,
            read_idx,
            rep_length=rep_length,
            int_length=int_length,
            strand=strand,
        )

        # Mappa la formula usando il Read_ID come chiave
        read_formulas_for_report[read_name] = formula_str

        # Dynamic color mapping for draw_dna_gene
        section_colors = []
        for _, row in clean_df.iterrows():
            stype = row["Type"]
            smotif = row["Motif"]
            if "Repeat" in stype:
                section_colors.append("#E69F00")  # Orange: Repeat
            elif "Interruption" in stype:
                if smotif in all_int_motifs:
                    m_idx = all_int_motifs.index(smotif)
                    section_colors.append(palette[m_idx % len(palette)])
                else:
                    section_colors.append("#0072B2")
            else:
                section_colors.append("#4D4D4D")  # Grey: Flanking/Other

        print(f"--> Generating linear map for Read: {read_name}")

        # Genera la mappa lineare per la read
        draw_dna_gene(
            sections_lengths=clean_df["Length"],
            section_colors=section_colors,
            str_identifier=gene_key,
            sample=args.sample,
            path=args.outdir,
            motifs_list=clean_df["Motif"].tolist(),
            ids=f"Read_{read_name}",
        )



    # --- 8. cohort Summary Text File --------------------------------------
    txt_file_path = os.path.join(args.outdir, f"{args.sample}_{gene_key}_summary.txt")
    with open(txt_file_path, "w") as text_file:
        for k in range(len(complete_STR)):
            formula, _ = parse_str_structure_string(
                complete_STR, k, rep_length=rep_length, int_length=int_length, strand=strand
            )
            line = f"{df_results.loc[k, 'Read_ID']}: {formula}\n"
            text_file.write(line)
    print(f"\nSummary text file successfully saved to: {txt_file_path}")

    # --- 9. Clinical Report & Sequential Structure Plot --------------------
    generate_clinical_report_with_structures(
        gene_symbol=gene_key,
        df_reads=df_results,
        peaks_summary=peaks_summary,
        allele_formulas=allele_formulas_for_report,
        sample_id=args.sample,
        path=args.outdir,
        database=STR_REFERENCE_DB,
    )

    plot_sequential_locus_structure(df_results, motif, args.outdir, interrupt_motif=interrupt_motif, strand=strand)

    print(f"\nAnalysis complete. Outputs available in: {args.outdir}")

    report_markdown = generate_full_diagnostic_report(
        gene_symbol=gene_key,
        df_reads=df_results,
        peaks_summary=peaks_summary,
        allele_formulas=allele_formulas_for_report,
        sample_id=args.sample,
        path=args.outdir,
        database=STR_REFERENCE_DB # Passa il tuo dizionario di riferimento del database STR
    )


if __name__ == "__main__":
    main()