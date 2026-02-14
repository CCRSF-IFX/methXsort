#!/usr/bin/env python3

import toolshed
from toolshed import nopen, reader, is_newer_b
import argparse
import sys
import os
from itertools import groupby
import gzip
import pysam
from multiprocessing import Pool, cpu_count
import subprocess
import shlex
import pathlib
import re

__version__ = "0.1.0"

def wrap(seq, width=60):
    return [seq[i:i+width] for i in range(0, len(seq), width)]

def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        yield header, "".join(s.strip() for s in next(faiter)).upper()
        

def convert_fasta(ref_fasta, out_fa=None):
    """
    Code originally copied from: bwa-meth and adapted for use in methxsort.py
    Convert a fasta file to a bwameth compatible format.
    """
    if out_fa is None:
        out_fa = ref_fasta + ".converted"
    msg = "c2t/g2a in %s to %s" % (ref_fasta, out_fa)
    if is_newer_b(ref_fasta, out_fa):
        sys.stderr.write("already converted: %s\n" % msg)
        return out_fa
    sys.stderr.write("converting %s\n" % msg)
    try:
        fh = open(out_fa, "w")
        for header, seq in fasta_iter(ref_fasta):
            ########### Reverse ######################
            fh.write(">r%s\n" % header)
            for line in wrap(seq.replace("G", "A")):
                fh.write(line + '\n')

            ########### Forward ######################
            fh.write(">f%s\n" % header)
            for line in wrap(seq.replace("C", "T")):
                fh.write(line + '\n')
        fh.close()
    except:
        try:
            fh.close()
        except UnboundLocalError:
            pass
        os.unlink(out_fa)
        raise
    return out_fa


def convert_reads_c2t_r1_g2a_r2(read, read2=None, out=None, out2=None, with_orig_seq=False):
    """
    Fast conversion using awk+sed+gzip for large FASTQ files.
    C->T in read (single or read1), G->A in read2 (if paired-end).
    If with_orig_seq is True, store the original sequence in the header line as an extra field using awk.
    """
    if out is None:
        out = read + ".meth"
    if read2 and out2 is None:
        out2 = read2 + ".meth"

    if with_orig_seq:
        # For read1: C->T, add original sequence to header (awk command)
        cmd1 = (
            f"zcat {read} | "
            "awk 'NR%4==1{{getline seq; print $0 \" ORIGINAL_SEQ:\" seq; print seq; getline; print $0; getline; print $0}}' | "
            "sed '2~4s/C/T/g;2~4s/c/t/g' | gzip > {out}"
        ).format(read=read, out=out)
        print(f"[convert-reads] CMD1: {cmd1}", file=sys.stdout)

        # For read2: G->A, add original sequence to header (awk command)
        if read2 and out2:
            cmd2 = (
                f"zcat {read2} | "
                "awk 'NR%4==1{{getline seq; print $0 \" ORIGINAL_SEQ:\" seq; print seq; getline; print $0; getline; print $0}}' | "
                "sed '2~4s/G/A/g;2~4s/g/a/g' | gzip > {out2}"
            ).format(read2=read2, out2=out2)
            print(f"[convert-reads] CMD2: {cmd2}", file=sys.stdout)
            # Run both commands in parallel
            p1 = subprocess.Popen(cmd1, shell=True)
            p2 = subprocess.Popen(cmd2, shell=True)
            p1.wait()
            p2.wait()
            if p1.returncode != 0 or p2.returncode != 0:
                raise subprocess.CalledProcessError(p1.returncode if p1.returncode != 0 else p2.returncode, cmd1 if p1.returncode != 0 else cmd2)
            return out, out2
        else:
            subprocess.check_call(cmd1, shell=True)
            return out
    else:
        # Default: fast sed-only implementation (no original sequence in header)
        cmd1 = (
            f"zcat {read} | "
            "sed '2~4s/C/T/g;2~4s/c/t/g' | gzip > {out}"
        ).format(read=read, out=out)
        print(f"[convert-reads] CMD1: {cmd1}", file=sys.stdout)

        if read2 and out2:
            cmd2 = (
                f"zcat {read2} | "
                "sed '2~4s/G/A/g;2~4s/g/a/g' | gzip > {out2}"
            ).format(read2=read2, out2=out2)
            print(f"[convert-reads] CMD2: {cmd2}", file=sys.stdout)
            # Run both commands in parallel
            p1 = subprocess.Popen(cmd1, shell=True)
            p2 = subprocess.Popen(cmd2, shell=True)
            p1.wait()
            p2.wait()
            if p1.returncode != 0 or p2.returncode != 0:
                raise subprocess.CalledProcessError(p1.returncode if p1.returncode != 0 else p2.returncode, cmd1 if p1.returncode != 0 else cmd2)
            return out, out2
        else:
            subprocess.check_call(cmd1, shell=True)
            return out

# Old Python implementation (commented out for reference)
# def convert_reads_c2t_r1_g2a_r2(read, read2=None, out=None, out2=None):
#     """
#     Convert C->T in read (single or read1) and G->A in read2 (if paired-end).
#     Store the original sequence in the header line as an extra field.
#     """
#     def convert_r1(seq):
#         return seq.replace("C", "T").replace("c", "t")
#     def convert_r2(seq):
#         return seq.replace("G", "A").replace("g", "a")
#     if out is None:
#         out = read + ".meth"
#     if read2 and out2 is None:
#         out2 = read2 + ".meth"

#     def open_out(filename):
#         if filename and filename.endswith('.gz'):
#             return gzip.open(filename, "wt")
#         else:
#             return open(filename, "w")

#     with nopen(read) as r1, open_out(out) as o1:
#         while True:
#             header = r1.readline()
#             if not header:
#                 break
#             seq = r1.readline()
#             plus = r1.readline()
#             qual = r1.readline()
#             # Add original sequence to header
#             header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
#             o1.write(header)
#             o1.write(convert_r1(seq))
#             o1.write(plus)
#             o1.write(qual)
#     if read2:
#         with nopen(read2) as r2, open_out(out2) as o2:
#             while True:
#                 header = r2.readline()
#                 if not header:
#                     break
#                 seq = r2.readline()
#                 plus = r2.readline()
#                 qual = r2.readline()
#                 header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
#                 o2.write(header)
#                 o2.write(convert_r2(seq))
#                 o2.write(plus)
#                 o2.write(qual)
#     return out, out2 if read2 else out

def bam_to_fastq_with_original_seq(bamfile, out1, out2=None):
    """
    Convert BAM to FASTQ, using the original sequence if present in the header (e.g., as ORIGINAL_SEQ:...), 
    otherwise use 'N' * read length. Output read 1 to out1, read 2 to out2 (if paired-end).
    """
    def open_out(filename):
        if filename and filename.endswith('.gz'):
            return gzip.open(filename, "wt")
        else:
            return open(filename, "w")

    fq1 = open_out(out1)
    fq2 = open_out(out2) if out2 else None

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            # Determine if read is read1 or read2
            if read.is_read1:
                fq = fq1
            elif read.is_read2 and fq2:
                fq = fq2
            else:
                fq = fq1
            parts = [read.query_name]
            # Try to extract original sequence from header
            orig_seq = None
            if "ORIGINAL_SEQ:" in read.query_name:
                # e.g. @name ORIGINAL_SEQ:ACGT...
                parts = read.query_name.split(" ORIGINAL_SEQ:")
                if len(parts) > 1:
                    orig_seq = parts[1].strip()
            # If not found, use N's of the same length as the read
            if not orig_seq or orig_seq == "":
                orig_seq = "N" * read.query_length
            fq.write(f"@{parts[0]}\n")
            fq.write(f"{orig_seq}\n")
            fq.write("+\n")
            qual = read.qual if read.qual else "I" * read.query_length
            if read.is_reverse:
                qual = qual[::-1]
            fq.write(f"{qual}\n")

    fq1.close()
    if fq2:
        fq2.close()

def process_chunk(records):
    out = []
    for header, seq, plus, qual in records:
        header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
        seq = seq.replace("C", "T").replace("c", "t")
        out.extend([header, seq, plus, qual])
    return out

def read_fastq_chunk(fh, chunk_size):
    records = []
    for _ in range(chunk_size):
        header = fh.readline()
        if not header:
            break
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        records.append((header, seq, plus, qual))
    return records

def process_chunk_r1(records):
    out = []
    for header, seq, plus, qual in records:
        header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
        seq = seq.replace("C", "T").replace("c", "t")
        out.extend([header, seq, plus, qual])
    return out

def process_chunk_r2(records):
    out = []
    for header, seq, plus, qual in records:
        header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
        seq = seq.replace("G", "A").replace("g", "a")
        out.extend([header, seq, plus, qual])
    return out

def parallel_convert_fastq(infile, outfile, chunk_size=100000, mode="r1"):
    opener = gzip.open if infile.endswith('.gz') else open
    process_chunk = process_chunk_r1 if mode == "r1" else process_chunk_r2
    with opener(infile, 'rt') as fh, gzip.open(outfile, 'wt') if outfile.endswith('.gz') else open(outfile, 'w') as out_fh:
        pool = Pool(cpu_count())
        while True:
            chunk = read_fastq_chunk(fh, chunk_size)
            if not chunk:
                break
            results = pool.map(process_chunk, [chunk])
            for recs in results:
                out_fh.writelines(recs)
        pool.close()
        pool.join()

def process_bam_chunk(records):
    out = []
    for read in records:
        if read.is_unmapped:
            continue
        parts = [read.query_name]
        # Try to extract original sequence from header
        orig_seq = None
        if "ORIGINAL_SEQ:" in read.query_name:
            parts = read.query_name.split(" ORIGINAL_SEQ:")
            if len(parts) > 1:
                orig_seq = parts[1].strip()
        if not orig_seq or orig_seq == "":
            sys.stderr.write(f"Warning: No original sequence found for read {read.query_name}.\n")
            sys.exit(1)
            #orig_seq = "N" * read.query_length
        fq_entry = f"@{parts[0]}\n{orig_seq}\n+\n"
        qual = read.qual if read.qual else "I" * read.query_length
        if read.is_reverse:
            qual = qual[::-1]
        fq_entry += f"{qual}\n"
        out.append((read.is_read1, fq_entry))
    return out

def read_bam_chunk(bam, chunk_size):
    records = []
    try:
        for _ in range(chunk_size):
            read = next(bam)
            records.append(read)
    except StopIteration:
        pass
    return records

# def parallel_bam_to_fastq(bamfile, out1, out2=None, chunk_size=100000):
#     def open_out(filename):
#         if filename and filename.endswith('.gz'):
#             return gzip.open(filename, "wt")
#         else:
#             return open(filename, "w")

#     fq1 = open_out(out1)
#     fq2 = open_out(out2) if out2 else None

#     with pysam.AlignmentFile(bamfile, "rb") as bam:
#         for read in bam.fetch(until_eof=True):
#             if read.is_unmapped:
#                 continue
#             # Determine if read is read1 or read2
#             if read.is_read1:
#                 fq = fq1
#             elif read.is_read2 and fq2:
#                 fq = fq2
#             else:
#                 fq = fq1
#             parts = [read.query_name]
#             # Try to extract original sequence from header
#             orig_seq = None
#             if "ORIGINAL_SEQ:" in read.query_name:
#                 parts = read.query_name.split(" ORIGINAL_SEQ:")
#                 if len(parts) > 1:
#                     orig_seq = parts[1].strip()
#             if not orig_seq or orig_seq == "":
#                 orig_seq = "N" * read.query_length
#             fq.write(f"@{parts[0]}\n")
#             fq.write(f"{orig_seq}\n")
#             fq.write("+\n")
#             qual = read.qual if read.qual else "I" * read.query_length
#             if read.is_reverse:
#                 qual = qual[::-1]
#             fq.write(f"{qual}\n")

#     fq1.close()
#     if fq2:
#         fq2.close()

def run_bbsplit(read1, read2, host, graft, out_host, out_graft, 
                bbsplit_index_build=1, bbsplit_index_path="bbsplit_index",
                bbsplit_path="bbsplit.sh", bbsplit_extra=""):
    """
    Run bbsplit.sh to split reads and require BAM output files.
    # ambiguous=best ambiguous2=split
    # AMBIGUOUS_bbsplit_graft.bam  AMBIGUOUS_bbsplit_host.bam
    """
    if not out_host.lower().endswith(".bam") or not out_graft.lower().endswith(".bam"):
        raise ValueError("Output files must have .bam suffix.")
    
    # Check for AMBIGUOUS BAM files before running bbsplit
    ambiguous_host = "AMBIGUOUS_" + os.path.basename(out_host)
    ambiguous_graft = "AMBIGUOUS_" + os.path.basename(out_graft)
    for amb_file in [ambiguous_host, ambiguous_graft]:
        if os.path.exists(amb_file):
            sys.stderr.write(f"Error: {amb_file} already exists. Please remove or rename it before proceeding.\n")
            sys.exit(1)

    cmd = (
        f"{bbsplit_path} build={bbsplit_index_build} "
        f"path={bbsplit_index_path} "
        f"in={read1} "
        f"{f'in2={read2} ' if read2 else ''}"
        f"out_{host}={out_host} "
        f"out_{graft}={out_graft} "
        f"{bbsplit_extra}"
    )
    print(f"Running: {cmd}", file=sys.stderr)
    subprocess.check_call(shlex.split(cmd))

def build_bbsplit_index(host_fa, graft_fa, host_name, graft_name, bbsplit_idx_dir="bbsplit_idx_convert", bbsplit_path="bbsplit.sh", build="1"):
    """
    Build bbsplit index using the provided converted reference FASTA files and custom names.
    """
    cmd = (
        f"{bbsplit_path} path={bbsplit_idx_dir} build={build} "
        f"ref_{host_name}={host_fa} ref_{graft_name}={graft_fa}"
    )
    print(f"[bbsplit-build] CMD: {cmd}", file=sys.stdout)
    subprocess.check_call(cmd, shell=True)

def ensure_bbsplit_index_structure(bbsplit_index_path, build="1", host=None, graft=None):
    """
    1. Ensure the bbsplit index folder structure exists:
    <bbsplit_index_path>/ref/genome/<build>
    <bbsplit_index_path>/ref/index/<build>
    2. Also ensure that host and graft names exist in namelist.txt.
    """
    base = pathlib.Path(bbsplit_index_path)
    genome_dir = base / "ref" / "genome" / str(build)
    index_dir = base / "ref" / "index" / str(build)
    namelist_file = genome_dir / "namelist.txt"

    # Check directory structure
    if not genome_dir.is_dir() or not index_dir.is_dir():
        raise FileNotFoundError(
            f"Required bbsplit index directories do not exist:\n"
            f"  {genome_dir}\n"
            f"  {index_dir}\n"
            "Please build the bbsplit index first."
        )

    # Check host and graft names in namelist.txt if provided
    if host or graft:
        if not namelist_file.is_file():
            raise FileNotFoundError(f"namelist.txt not found: {namelist_file}")
        with open(namelist_file) as f:
            names = set(line.strip() for line in f if line.strip())
        missing = [n for n in (host, graft) if n and n not in names]
        if missing:
            raise ValueError(
                f"The following names are missing in {namelist_file}: {', '.join(missing)}"
            )

def fastq_count(filename):
    """
    Fast count of reads in a FASTQ file (gzipped or plain).
    """
    import subprocess
    import os

    if filename.endswith('.gz'):
        cmd = f"zcat {filename} | wc -l"
    else:
        cmd = f"wc -l < {filename}"
    n_lines = int(subprocess.check_output(cmd, shell=True).decode().strip())
    return n_lines // 4

def stat_split(raw_fastq, host_fastq, graft_fastq, xengsort_prefix=None, sample_name=None):
    """
    Output CSV with sample name, raw read number, host read number, graft read number, percent graft.
    If xengsort_prefix is provided, also count both/neither/ambiguous reads from xengsort output.
    """
    import os
    n_raw = fastq_count(raw_fastq)
    n_host = fastq_count(host_fastq)
    n_graft = fastq_count(graft_fastq)
    percent_graft = (n_graft / n_raw * 100) if n_raw else 0

    print("Sample_name,raw_read_number,host_read_number,graft_read_number,percent_graft", end="")
    if xengsort_prefix:
        print(",both_read_number,neither_read_number,ambiguous_read_number")
    else:
        print()

    print(f"{sample_name},{n_raw},{n_host},{n_graft},{percent_graft:.2f}", end="")

    if xengsort_prefix:
        def count_xengsort_reads(suffix):
            fname = f"{xengsort_prefix}-{suffix}.1.fq.gz"
            return fastq_count(fname) if os.path.exists(fname) else "NA"

        n_both = count_xengsort_reads("both")
        n_neither = count_xengsort_reads("neither")
        n_ambiguous = count_xengsort_reads("ambiguous")
        print(f",{n_both},{n_neither},{n_ambiguous}")
    else:
        print()

def filter_fastq_by_bam(read, bamfile, out, filterbyname_path="filterbyname.sh", read2=None, out2=None):
    """
    Extract read IDs from BAM, save to a temp file, and use filterbyname.sh to extract reads from FASTQ.
    The read ID file is named based on the output FASTQ basename with .readids suffix.
    If out2 and read2 are provided, treat as paired-end and run filterbyname.sh for both read1/read2.
    Otherwise, run for single-end.
    Both --read and --out are required and cannot be None.
    """

    def get_readids_filename(out):
        base = os.path.basename(out)
        base = re.sub(r'\.(fastq|fq)(\.gz)?$', '', base)
        base_dir = os.path.dirname(out)
        return os.path.join(base_dir, base + ".readids")

    # Extract read IDs for read1 and read2
    readids1 = get_readids_filename(out)
    with open(readids1, "w") as f:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read_bam in bam.fetch(until_eof=True):
                if read_bam.is_unmapped:
                    continue
                if out2:
                    if read_bam.is_read1:
                        f.write(read_bam.query_name + "\n")
                else:
                    f.write(read_bam.query_name + "\n")
    # Run filterbyname.sh for single-end or paired-end
    if out2 and read2:
        # Paired-end
        cmd = (
            f"{filterbyname_path} in={read} in2={read2} "
            f"out={out} out2={out2} names={readids1} include=t"
        )
        print(f"[bam-to-fastq] CMD: {cmd}", file=sys.stdout)
        subprocess.check_call(cmd, shell=True)
    else:
        # Single-end
        cmd = (
            f"{filterbyname_path} in={read} out={out} names={readids1} include=t"
        )
        print(f"[bam-to-fastq] CMD: {cmd}", file=sys.stdout)
        subprocess.check_call(cmd, shell=True)

def filter_xengsort_extra(xengsort_extra, used_params):
    """
    Remove any parameters from xengsort_extra that are already provided by main args.
    If any are found, error out with a clear message.
    """
    tokens = shlex.split(xengsort_extra)
    for i, token in enumerate(tokens):
        # Check for direct match or --param=value style
        for pat in used_params:
            if token == pat or (pat.startswith("--") and token.startswith(pat + "=")):
                sys.stderr.write(
                    f"Error: Parameter '{pat}' is already provided by main arguments and should not be included in --xengsort_extra.\n"
                )
                sys.exit(1)
    return xengsort_extra

def run_xengsort_classify(args):
    """
    Run xengsort classify with argument checking for redundant parameters in xengsort_extra.
    """
    used_params = [
        "--fastq", "-q", "--pairs", "-p", "--index", "--out", 
        "-o", "--prefix", "--threads", "-T", "-j"
    ]
    filtered_extra = filter_xengsort_extra(args.xengsort_extra, used_params)
    cmd = (
        f"{args.xengsort_path} classify "
        f"--fastq {args.read} "
        f"{f'--pairs {args.read2} ' if args.read2 else ''}"
        f"--index {args.index} "
        f"--out {args.out_prefix} "
        f"--threads {args.threads} "
        f"{filtered_extra}"
    )
    print(f"[xengsort-classify] CMD: {cmd}", file=sys.stdout)
    subprocess.check_call(cmd, shell=True)

def restore_fastq_from_xengsort(read, out, read2=None, out2=None):
    """
    Restore original sequences in FASTQ files classified by xengsort.
    The input FASTQ should have the original sequence in the header as 'ORIGINAL_SEQ:...'.
    The output FASTQ will have the original sequence as the sequence line, and the header will be cleaned.
    If read2 and out2 are provided, process both files in parallel using paste+awk for speed.
    """
    import subprocess

    if read2 and out2:
        cmd1 = (
            f"zcat {read} | "
            "awk 'NR%4==1{{split($0, h, \" ORIGINAL_SEQ:\"); header=h[1]; orig_seq=h[2]}} "
            "NR%4==2{{seq=$0}} NR%4==3{{plus=$0}} NR%4==0{{qual=$0; print header \"\\n\" orig_seq \"\\n\" plus \"\\n\" qual}}' | gzip > {out}"
        ).format(out=out)
        cmd2 = (
            f"zcat {read2} | "
            "awk 'NR%4==1{{split($0, h, \" ORIGINAL_SEQ:\"); header=h[1]; orig_seq=h[2]}} "
            "NR%4==2{{seq=$0}} NR%4==3{{plus=$0}} NR%4==0{{qual=$0; print header \"\\n\" orig_seq \"\\n\" plus \"\\n\" qual}}' | gzip > {out2}"
        ).format(out2=out2)
        print(f"[restore-fastq] CMD1: {cmd1}", file=sys.stdout)
        print(f"[restore-fastq] CMD2: {cmd2}", file=sys.stdout)
        p1 = subprocess.Popen(cmd1, shell=True, executable="/bin/bash")
        p2 = subprocess.Popen(cmd2, shell=True, executable="/bin/bash")
        p1.wait()
        p2.wait()
        if p1.returncode != 0 or p2.returncode != 0:
            raise subprocess.CalledProcessError(p1.returncode if p1.returncode != 0 else p2.returncode, cmd1 if p1.returncode != 0 else cmd2)
    else:
        # Single-end: use awk
        cmd = (
            f"zcat {read} | "
            "awk 'NR%4==1{{split($0, h, \" ORIGINAL_SEQ:\"); header=h[1]; orig_seq=h[2]}} "
            "NR%4==2{{seq=$0}} NR%4==3{{plus=$0}} NR%4==0{{qual=$0; print header \"\\n\" orig_seq \"\\n\" plus \"\\n\" qual}}' | gzip > {out}"
        ).format(out=out)
        print(f"[restore-fastq] CMD: {cmd}", file=sys.stdout)
        subprocess.check_call(cmd, shell=True, executable="/bin/bash")

# def parse_xengsort_summary(logfile, outfile):
#     """
#     Parse xengsort log file and write the Classification Statistics table as CSV.
#     """
#     import csv

#     stats_started = False
#     header = []
#     row = []
#     with open(logfile) as f:
#         for line in f:
#             line = line.strip()
#             if not stats_started:
#                 if line.startswith("prefix") and "host" in line and "graft" in line:
#                     header = line.split()
#                     stats_started = True
#             elif stats_started and line and not line.startswith("|"):
#                 # This is the data row
#                 row = line.split()
#                 break  # Only the first data row is needed

#     if not header or not row:
#         sys.stderr.write("Error: Could not find Classification Statistics in log file.\n")
#         sys.exit(1)

#     # Write to CSV
#     with open(outfile, "w", newline="") as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(header)
#         writer.writerow(row)

#     print(f"Classification statistics written to {outfile}")

def parse_xengsort_summary_table(logfile, outfile):
    """
    Parse the vertical table from xengsort log and write as CSV.
    Only the section with | prefix | host | graft | ambiguous | both | neither | is parsed.
    """
    import csv

    data = {}
    with open(logfile) as f:
        for line in f:
            line = line.strip()
            if line.startswith("| prefix"):
                # Prefix line
                parts = line.split("|")
                if len(parts) > 2:
                    data["prefix"] = parts[2].strip().split("/")[-1]  # Get the last part of the prefix
            elif line.startswith("| host"):
                parts = line.split("|")
                if len(parts) > 2:
                    data["host"] = parts[2].strip()
                    data["host_pct"] = parts[3].strip(" %|")
            elif line.startswith("| graft"):
                parts = line.split("|")
                if len(parts) > 2:
                    data["graft"] = parts[2].strip()
                    data["graft_pct"] = parts[3].strip(" %|")
            elif line.startswith("| ambiguous"):
                parts = line.split("|")
                if len(parts) > 2:
                    data["ambiguous"] = parts[2].strip()
                    data["ambiguous_pct"] = parts[3].strip(" %|")
            elif line.startswith("| both"):
                parts = line.split("|")
                if len(parts) > 2:
                    data["both"] = parts[2].strip()
                    data["both_pct"] = parts[3].strip(" %|")
            elif line.startswith("| neither"):
                parts = line.split("|")
                if len(parts) > 2:
                    data["neither"] = parts[2].strip()
                    data["neither_pct"] = parts[3].strip(" %|")

    # Write to CSV
    header = ["prefix", "host", "host_pct", "graft", "graft_pct", "ambiguous", "ambiguous_pct", "both", "both_pct", "neither", "neither_pct"]
    with open(outfile, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow([data.get(col, "") for col in header])

    print(f"Classification statistics (vertical table) written to {outfile}")

if __name__ == "__main__":
