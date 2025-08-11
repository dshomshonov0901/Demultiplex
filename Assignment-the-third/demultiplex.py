#!/usr/bin/env python3


import argparse
import os
import gzip


def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp[b] for b in reversed(seq))

def load_indexes(index_file: str) -> set:
    """
    Read a one-column index file (one sequence per line)
    and return a set of valid index sequences (uppercased).
    """
    valid = set()
    with open(index_file) as f:
        for line in f:
            s = line.strip()
            if s:
                valid.add(s.upper())
    return valid

def write_read_block(out_r1, out_r2, r1_block, r2_block, index_info):
    """Append index-pair to headers and write both reads."""
    r1_block[0] = r1_block[0].rstrip("\n") + f" {index_info}\n"
    r2_block[0] = r2_block[0].rstrip("\n") + f" {index_info}\n"
    out_r1.writelines(r1_block)
    out_r2.writelines(r2_block)

def demultiplex(r1_fp, r2_fp, i1_fp, i2_fp, index_fp, out_dir):
    """
    Demultiplex reads into:
      - per-index matched files: matched_<INDEXSEQ>_R1.fq / _R2.fq
      - index_hop_R1.fq / _R2.fq for valid but mismatched pairs
      - unknown_R1.fq / _R2.fq for Ns or indexes not in the list
    Also writes:
      - pair_counts.tsv (all valid i1→i2rc pairs with category)
      - hopped_pairs.tsv (valid hopped only)
      - pair_matrix.tsv (24×24 matrix of valid pairs)
    """
    valid = load_indexes(index_fp)
    os.makedirs(out_dir, exist_ok=True)

    # lazily-open per-index matched files
    matched_files = {}
    def get_matched_handles(seq):
        if seq not in matched_files:
            base = os.path.join(out_dir, f"matched_{seq}")
            matched_files[seq] = (
                open(f"{base}_R1.fq", "w"),
                open(f"{base}_R2.fq", "w"),
            )
        return matched_files[seq]

    # shared category outputs
    hop_r1 = open(os.path.join(out_dir, "index_hop_R1.fq"), "w")
    hop_r2 = open(os.path.join(out_dir, "index_hop_R2.fq"), "w")
    unk_r1 = open(os.path.join(out_dir, "unknown_R1.fq"), "w")
    unk_r2 = open(os.path.join(out_dir, "unknown_R2.fq"), "w")

    # topline counts
    counts = {"unknown": 0, "hopped": 0}
    for seq in valid:
        counts[f"matched_{seq}"] = 0

    # pair-level counts
    pair_counts = {}         # (i1, i2rc) -> int
    hopped_pair_counts = {}  # (i1, i2rc) -> int
 # hopped-only


    with gzip.open(r1_fp, "rt") as r1, \
         gzip.open(r2_fp, "rt") as r2, \
         gzip.open(i1_fp, "rt") as i1, \
         gzip.open(i2_fp, "rt") as i2:

        r1_rl, r2_rl, i1_rl, i2_rl = r1.readline, r2.readline, i1.readline, i2.readline
        while True:
            # read 4-line blocks
            r1_block = [r1_rl(), r1_rl(), r1_rl(), r1_rl()]
            if not r1_block[0]:
                break  # EOF
            r2_block = [r2_rl(), r2_rl(), r2_rl(), r2_rl()]
            i1_block = [i1_rl(), i1_rl(), i1_rl(), i1_rl()]
            i2_block = [i2_rl(), i2_rl(), i2_rl(), i2_rl()]

            i1_seq = i1_block[1].rstrip("\n").upper()
            i2_seq = rev_comp(i2_block[1].rstrip("\n").upper())  # reverse complement i2
            pair_label = f"{i1_seq}-{i2_seq}"

            # Unknown: any N or not in whitelist
            if "N" in i1_seq or "N" in i2_seq:
                write_read_block(unk_r1, unk_r2, r1_block, r2_block, pair_label)
                counts["unknown"] += 1
                continue

            i1_valid = i1_seq in valid
            i2_valid = i2_seq in valid

            key = (i1_seq, i2_seq)

            if i1_valid and i2_valid:
                if i1_seq == i2_seq:
                    # matched
                    out_r1, out_r2 = get_matched_handles(i1_seq)
                    write_read_block(out_r1, out_r2, r1_block, r2_block, pair_label)
                    counts[f"matched_{i1_seq}"] += 1
                else:
                    # hopped
                    write_read_block(hop_r1, hop_r2, r1_block, r2_block, pair_label)
                    counts["hopped"] += 1
                    hopped_pair_counts[key] = hopped_pair_counts.get(key, 0) + 1

                # count ALL valid pairs (matched + hopped)
                pair_counts[key] = pair_counts.get(key, 0) + 1
            else:
                write_read_block(unk_r1, unk_r2, r1_block, r2_block, pair_label)
                counts["unknown"] += 1



    # close outputs
    for f1, f2 in matched_files.values():
        f1.close(); f2.close()
    hop_r1.close(); hop_r2.close()
    unk_r1.close(); unk_r2.close()

    # write pair summaries
    write_pair_summaries(out_dir, valid, pair_counts, hopped_pair_counts)

    return counts

def write_pair_summaries(out_dir, valid, pair_counts, hopped_pair_counts):
    """Dump long tables and a 24×24 matrix for valid pairs."""
    # long tables
    pc_path = os.path.join(out_dir, "pair_counts.tsv")
    hp_path = os.path.join(out_dir, "hopped_pairs.tsv")
    with open(pc_path, "w") as fh:
        fh.write("i1\ti2\tcount\tcategory\n")
        for (a,b), c in sorted(pair_counts.items(), key=lambda x: (-x[1], x[0])):
            cat = "matched" if a == b else "hopped"
            fh.write(f"{a}\t{b}\t{c}\t{cat}\n")
    with open(hp_path, "w") as fh:
        fh.write("i1\ti2\tcount\n")
        for (a,b), c in sorted(hopped_pair_counts.items(), key=lambda x: (-x[1], x[0])):
            fh.write(f"{a}\t{b}\t{c}\n")

    # matrix (valid vs valid)
    idxs = sorted(valid)
    mx_path = os.path.join(out_dir, "pair_matrix.tsv")
    with open(mx_path, "w") as fh:
        fh.write("i1\\i2\t" + "\t".join(idxs) + "\n")
        for a in idxs:
            row = [a] + [str(pair_counts.get((a,b), 0)) for b in idxs]
            fh.write("\t".join(row) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Demultiplex dual-indexed FASTQ files.")
    parser.add_argument("-r1", required=True, help="Read 1 FASTQ (R1)")
    parser.add_argument("-r2", required=True, help="Read 2 FASTQ (R4)")
    parser.add_argument("-i1", required=True, help="Index 1 FASTQ (R2)")
    parser.add_argument("-i2", required=True, help="Index 2 FASTQ (R3; will be reverse-complemented)")
    parser.add_argument("-idx", required=True, help="Index file: ONE sequence per line")
    parser.add_argument("-o", required=True, help="Output directory")
    args = parser.parse_args()

    stats = demultiplex(args.r1, args.r2, args.i1, args.i2, args.idx, args.o)

    print("\n--- Summary ---")
    for k, v in sorted(stats.items()):
        print(f"{k}: {v}")

    # quick peek at top 10 hopped pairs (handy in slurm log)
    # (file has the full list)
    # You can comment this out if you want totally quiet logs.
    # Read back hopped_pairs for convenience would require state; just hint user.
    print("\nDetails written to: pair_counts.tsv, hopped_pairs.tsv, pair_matrix.tsv")

if __name__ == "__main__":
    main()