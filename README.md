# SF-methXsort

**methXsort** is a command-line toolkit for sorting bisulfite sequencing reads into host species and graft species in xenograft experiments. 

MethXsort supports both `xengsort` and `bbsplit` for sorting reads into host and graft species. We recommend `xengsort`, as it is accurate and much faster based on benchmarking results.

<img src="docs/methXsort_flowchart.png" width="400">

---

## Installation

### Requirements
- Python >= 3.12
- pysam >= 0.15.0
- xengsort >= 2.0.9

### Install from source

```bash
git clone https://github.com/CCRSF-IFX/methXsort.git
cd methXsort
pip install -e .
```

See [INSTALLATION.md](INSTALLATION.md) for more details.

---

## Usage

All commands are run via the installed command:

```bash
methXsort <subcommand> [options]
```

### Main Subcommands

#### Convert Reference Genome

Convert a reference FASTA for bisulfite mapping (C→T and G→A):

```bash
methXsort convert-ref <ref_fasta> [-o OUTPUT]
```

#### Build xengsort Index

```bash
methXsort xengsort-index --host <host.fa> --graft <graft.fa> --index <index_dir> [-n N] [--fill FILL] [--statistics STAT] [-k K] [--xengsort_path <path>] [--xengsort_extra <extra>]
```

#### Convert Reads

Convert reads for bisulfite mapping (C→T for R1, G→A for R2):

```bash
methXsort convert-reads --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] [--out <R1_out>] [--out2 <R2_out>] [--with_orig_seq]
```
- `--with_orig_seq`: Store the original sequence in the header (slower, but traceable).


#### Classify Reads with xengsort

```bash
methXsort xengsort-classify --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] --index <index_dir> --out_prefix <prefix> --threads <N> [--xengsort_path <path>] [--xengsort_extra <extra>]
```

Output: 

* {prefix}.host.1.fq.gz: host reads

* {prefix}.graft.1.fq.gz: graft reads

* {prefix}.both.1.fq.gz: reads that could originate from both

* {prefix}.neither.1.fq.gz: reads that originate from neither host nor graft

* {prefix}.ambiguous.1.fq.gz: (few) ambiguous reads that cannot be classified,


####  Split Statistics

Output CSV statistics for split reads:

```bash
methXsort stat-split --raw <raw_R1.fastq.gz> --host <host_R1.fastq.gz> --graft <graft_R1.fastq.gz>
```

#### Restore FASTQ from xengsort Output

Restore original sequences in FASTQ files classified by xengsort:

```bash
methXsort restore-fastq --read <classified_R1.fq.gz> --out <restored_R1.fq.gz> [--read2 <classified_R2.fq.gz> --out2 <restored_R2.fq.gz>]
```

---

## Example Workflow

1. **Convert reference genomes:**
    ```bash
    methXsort convert-ref mm10.fa -o mm10_converted.fa
    methXsort convert-ref hg38.fa -o hg38_converted.fa
    ```

2. **Build bbsplit and xengsort indices:**
    ```bash
    methXsort xengsort-index --host mm10_converted.fa --graft hg38_converted.fa --index xengsort_index_7B
    ```

3. **Convert reads:**
    ```bash
    methXsort convert-reads --read sample_R1.fastq.gz --read2 sample_R2.fastq.gz --with_orig_seq
    ```

4. **Run bbsplit or xengsort:**
    ```bash
    methXsort xengsort-classify --read sample_R1.meth.gz --read2 sample_R2.meth.gz --index xengsort_index_7B --out_prefix sample_xengsort --threads 8
    ```

5. **Restore original FASTQ:**
    ```bash
    methXsort restore-fastq --read sample_xengsort-graft.1.fq.gz --out sample_graft_R1_restored.fq.gz --read2 sample_xengsort-graft.2.fq.gz --out2 sample_graft_R2_restored.fq.gz
    ```

### Alternstive workflow using `bbsplit`


#### 5. Build bbsplit Index

```bash
methXsort bbsplit-index --host <host.fa> --graft <graft.fa> --host_name <host> --graft_name <graft> [--bbsplit_path <path>] [--bbsplit_index_path <dir>]
```

#### 4. Run bbsplit

Split reads into host and graft using bbsplit:

```bash
methXsort bbsplit --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] --host <host_name> --graft <graft_name> --out_host <host.bam> --out_graft <graft.bam> [--bbsplit_path <path>] [--bbsplit_extra <extra>]
```

#### 3. Filter FASTQ by BAM

Extract reads from FASTQ that are present in a BAM file (e.g., after bbsplit):

```bash
methXsort filter-fastq-by-bam --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] --bam <file.bam> --out <R1_out> [--out2 <R2_out>] [--filterbyname_path <path>]
```

---

## Notes

- For all subcommands, use `-h` or `--help` to see detailed options.
- Make sure all required external tools (`bbsplit.sh`, `filterbyname.sh`, `xengsort`) are in your `PATH` or specify their locations with the appropriate options.
- For paired-end data, always provide both `--read2` and `--out2` where required.

---

## Contact

Email: ccrsfifx@nih.gov

