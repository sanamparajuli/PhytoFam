// utils.nf
// Utility processes: sequence extraction and summary report

process EXTRACT_SEQUENCES {
    tag "extract_${label}"
    label 'low_cpu'
    publishDir "${params.outdir}/sequences", mode: 'copy'

    input:
    path proteome
    path id_list
    val  label

    output:
    path "${label}.fasta",         emit: sequences
    path "${label}_not_found.txt", emit: not_found

    script:
    """
    python3 << 'PYEOF'
label    = '${label}'
proteome = '${proteome}'
id_list  = '${id_list}'

import sys

with open(id_list) as fh:
    target_ids = set(line.strip() for line in fh if line.strip())

found     = set()
not_found = set(target_ids)
writing   = False

with open(proteome) as fh, \\
     open(f"{label}.fasta", "w") as out:
    for line in fh:
        if line.startswith(">"):
            seq_id  = line[1:].split()[0]
            writing = seq_id in target_ids
            if writing:
                found.add(seq_id)
                not_found.discard(seq_id)
        if writing:
            out.write(line)

with open(f"{label}_not_found.txt", "w") as out:
    for sid in sorted(not_found):
        print(sid, file=out)

if not_found:
    print(f"WARNING: {len(not_found)} IDs not found in proteome", file=sys.stderr)
    for sid in sorted(not_found):
        print(f"  Missing: {sid}", file=sys.stderr)

print(f"Extracted {len(found)}/{len(target_ids)} sequences for '{label}'",
      file=sys.stderr)

if len(found) == 0:
    print("ERROR: No sequences extracted. Check ID format.", file=sys.stderr)
    sys.exit(1)
    PYEOF
    """
}

process SUMMARY_REPORT {
    tag "summary_report"
    label 'low_cpu'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path hmmer_ids
    path confirmed_ids
    path final_ids
    path treefile

    output:
    path "pipeline_summary.txt", emit: summary
    path "final_gene_ids.txt",   emit: final_ids

    script:
    """
    python3 << 'PYEOF'
hmmer_ids_file     = '${hmmer_ids}'
confirmed_ids_file = '${confirmed_ids}'
final_ids_file     = '${final_ids}'
treefile           = '${treefile}'

import os
from datetime import datetime

def count_lines(filepath):
    try:
        with open(filepath) as fh:
            return sum(1 for line in fh if line.strip())
    except Exception:
        return 0

def read_ids(filepath):
    try:
        with open(filepath) as fh:
            return set(line.strip() for line in fh if line.strip())
    except Exception:
        return set()

n_hmmer     = count_lines(hmmer_ids_file)
n_confirmed = count_lines(confirmed_ids_file)
n_final     = count_lines(final_ids_file)
final_set   = read_ids(final_ids_file)
tree_ok     = os.path.exists(treefile) and os.path.getsize(treefile) > 10

with open("final_gene_ids.txt", "w") as out:
    for sid in sorted(final_set):
        print(sid, file=out)

def pct(a, b):
    return f"{a/b*100:.1f}%" if b > 0 else "N/A"

with open("pipeline_summary.txt", "w") as out:
    print("=" * 60, file=out)
    print("  Gene Family Identification Pipeline -- Summary", file=out)
    print("=" * 60, file=out)
    print(f"  Date : {datetime.now().strftime('%Y-%m-%d %H:%M')}", file=out)
    print("", file=out)
    print("  Step                                  Sequences", file=out)
    print("  " + "-" * 50, file=out)
    print(f"  1. HMMER (post isoform dedup)       : {n_hmmer:>6}", file=out)
    print(f"  2. InterProScan confirmed           : {n_confirmed:>6}"
          f" ({pct(n_confirmed, n_hmmer)} retained)", file=out)
    print(f"  3. Orthology assigned (no filter)   : {n_final:>6}"
          f" (= all confirmed)", file=out)
    print("", file=out)
    print(f"  Final gene family members           : {n_final}", file=out)
    print(f"  Phylogenetic tree produced          : {'Yes' if tree_ok else 'No'}", file=out)
    print("", file=out)
    print("  Note: RBH orthology assignments are in", file=out)
    print("        03_blast/orthology_table.tsv", file=out)
    print("", file=out)
    print("  Output files:", file=out)
    print("    sequences/final_candidates.fasta", file=out)
    print("    04_alignment/sequence_manifest.tsv (outgroup flags)", file=out)
    print("    06_phylogeny/phylogeny.treefile", file=out)
    print("    final_gene_ids.txt", file=out)
    print("=" * 60, file=out)

print("Pipeline complete. See pipeline_summary.txt")
    PYEOF
    """
}
