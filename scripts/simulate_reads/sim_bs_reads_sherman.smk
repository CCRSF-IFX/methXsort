
graft_fa = config["graft_fa"]
host_fa = config["host_fa"]
graft_name = os.path.basename(graft_fa).split(".")[0]
host_name = os.path.basename(host_fa).split(".")[0]
outdir = os.path.abspath(config["outdir"])
read_number = config["read_pair_number"]
bbsplit_path = config["bbsplit_idx_path"]
bbsplit_idx = config["bbsplit_idx"]

# Note: methXsort should be installed via pip (pip install methXsort)
# If not installed, you can still use: python /path/to/methXsort.py
sherman_path = config["sherman_path"]
stat_accuracy_path = config.get("stat_accuracy_path", None)

rule all:
    input:
        expand(os.path.join(outdir, "fastq/bbsplit_{genome}_{pair}.fastq.gz"), 
               genome=["graft", "host"], 
               pair =["R1", "R2"]),
        os.path.join(outdir, "read_number_stat.txt"),
        os.path.join(outdir, "accuracy_stat.txt")

rule sherman: 
    input: 
        graft = graft_fa,
        host = host_fa,
    params:
        od_graft = os.path.join(outdir, "graft"),
        od_host = os.path.join(outdir, "host"),
    output:
        graft = os.path.join(outdir, "graft/sherman.log"), 
        host = os.path.join(outdir, "host/sherman.log"),
        sim_graft_R1 = os.path.join(outdir, "sim_reads/simulated_bs_graft_R1.fastq.gz"),
        sim_graft_R2 = os.path.join(outdir, "sim_reads/simulated_bs_graft_R2.fastq.gz"),
        sim_host_R1 = os.path.join(outdir, "sim_reads/simulated_bs_host_R1.fastq.gz"),
        sim_host_R2 = os.path.join(outdir, "sim_reads/simulated_bs_host_R2.fastq.gz"),
    shell:
        """
mkdir -p {params.od_graft} && ln -s {input.graft} {params.od_graft}/graft.fa
mkdir -p {params.od_host} && ln -s {input.host} {params.od_host}/host.fa
cd  {params.od_graft} && \
 {sherman_path} --conversion_rate 60 -l 150 -n {read_number} --genome_folder {params.od_graft} -pe > {output.graft} 2>&1 & 
cd  {params.od_host} && \
 {sherman_path} --conversion_rate 60 -l 150 -n {read_number} --genome_folder {params.od_host} -pe > {output.host} 2>&1 &
wait
cd ../ && 
cat {params.od_graft}/simulated_1.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1graft_/' |gzip -c > {output.sim_graft_R1}
cat {params.od_graft}/simulated_2.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1graft_/' |gzip -c > {output.sim_graft_R2}
cat {params.od_host}/simulated_1.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1host_/' |gzip -c > {output.sim_host_R1}
cat {params.od_host}/simulated_2.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1host_/' |gzip -c > {output.sim_host_R2}
"""

rule merge_graft_host:
    input:
        graft_R1 = rules.sherman.output.sim_graft_R1,
        graft_R2 = rules.sherman.output.sim_graft_R2,
        host_R1 = rules.sherman.output.sim_host_R1,
        host_R2 = rules.sherman.output.sim_host_R2,
    output:
        merged_R1 = os.path.join(outdir, "sim_reads/merged_simulated_bs_R1.fastq.gz"),
        merged_R2 = os.path.join(outdir, "sim_reads/merged_simulated_bs_R2.fastq.gz"),
        merged_cvt_R1 = os.path.join(outdir, "sim_reads/merged_simulated_bscvt_R1.fastq.gz"),
        merged_cvt_R2 = os.path.join(outdir, "sim_reads/merged_simulated_bscvt_R2.fastq.gz"),
    shell:
        """
cat {input.graft_R1} {input.host_R1} > {output.merged_R1}
cat {input.graft_R2} {input.host_R2} > {output.merged_R2}
methXsort convert-reads \
    --read {output.merged_R1} --read2 {output.merged_R2} \
    --out {output.merged_cvt_R1} --out2 {output.merged_cvt_R2}
"""

rule convert_ref:
    input:
        graft_fa = graft_fa,
        host_fa = host_fa,
    output:
        graft_fa_cvt = os.path.join(outdir, "ref_idx/graft_cvt.fa"),
        host_fa_cvt = os.path.join(outdir, "ref_idx/host_cvt.fa"),
    shell:
        """
methXsort convert-ref \
    {input.graft_fa} --out {output.graft_fa_cvt} 
methXsort convert-ref \
    {input.host_fa} --out {output.host_fa_cvt}
"""

rule bbsplit_idx:
    input:
        graft_fa = rules.convert_ref.output.graft_fa_cvt,
        host_fa = rules.convert_ref.output.host_fa_cvt,
    output:
        log = os.path.join(outdir, "ref_idx/bbsplit_index.log"),
    params:
        bbsplit_index_path = os.path.join(outdir, "ref_idx/bbsplit_index"),
        build = 1

    shell:
        """
methXsort bbsplit-index \
    --host {input.host_fa} --graft {input.graft_fa} \
    --host_name {host_name} --graft_name {graft_name} \
    --bbsplit_index_path {params.bbsplit_index_path} \
    --bbsplit_index_build 1 > {output.log} 2>&1

"""

rule bbsplit:
    input:
        bbsplit_idx_log = rules.bbsplit_idx.output.log,
        reads_R1 = rules.merge_graft_host.output.merged_cvt_R1,
        reads_R2 = rules.merge_graft_host.output.merged_cvt_R2,
    params:
        bbsplit_index_path = os.path.join(outdir, "ref_idx/bbsplit_index"),
        bbsplit_idx = bbsplit_idx,
    output:
        bbsplit_log = os.path.join(outdir, "bbsplit.log"),
        bam_graft = os.path.join(outdir, "bbsplit/bbsplit_graft.bam"),
        bam_host = os.path.join(outdir, "bbsplit/bbsplit_host.bam"),
        scafstats = os.path.join(outdir, "bbsplit/scafstats.log"),
        refstats = os.path.join(outdir, "bbsplit/refstats.log"),
    shell:
        """
mkdir -p {outdir}/bbsplit && cd {outdir}/bbsplit && \
methXsort bbsplit \
    --read {input.reads_R1} --read2 {input.reads_R2} \
    --host {host_name} --graft {graft_name} \
    --bbsplit_index_build 1 \
    --bbsplit_index_path {params.bbsplit_index_path} \
    --out_host {output.bam_host} --out_graft {output.bam_graft} \
    --bbsplit_extra "ambiguous=best ambiguous2=split scafstats={output.scafstats} refstats={output.refstats}" \
    > {output.bbsplit_log} 2>&1
"""

rule convert_bam_to_fastq:
    input:
        R1 = rules.merge_graft_host.output.merged_R1,
        R2 = rules.merge_graft_host.output.merged_R2,
        bam_graft = rules.bbsplit.output.bam_graft,
        bam_host = rules.bbsplit.output.bam_host,
    output:
        fastq_graft_R1 = os.path.join(outdir, "fastq/bbsplit_graft_R1.fastq.gz"),
        fastq_graft_R2 = os.path.join(outdir, "fastq/bbsplit_graft_R2.fastq.gz"),
        fastq_host_R1 = os.path.join(outdir, "fastq/bbsplit_host_R1.fastq.gz"),
        fastq_host_R2 = os.path.join(outdir, "fastq/bbsplit_host_R2.fastq.gz"),
    shell:
        """
methXsort filter-fastq-by-bam \
    --read {input.R1} --read2 {input.R2} \
    --bam {input.bam_graft} \
    --out {output.fastq_graft_R1} --out2 {output.fastq_graft_R2} 
methXsort filter-fastq-by-bam \
    --read {input.R1} --read2 {input.R2} \
    --bam {input.bam_host} \
    --out {output.fastq_host_R1} --out2 {output.fastq_host_R2}
"""

rule stat_read_number: 
    input: 
        raw_R1 = rules.merge_graft_host.output.merged_cvt_R1,
        fastq_graft_R1 = rules.convert_bam_to_fastq.output.fastq_graft_R1,
        fastq_host_R1 = rules.convert_bam_to_fastq.output.fastq_host_R1
    output: 
        read_number_stat = os.path.join(outdir, "read_number_stat.txt")
    params: batch = "-l nodes=1:ppn=16,mem=64g"
    shell: 
        """
methXsort stat-split --raw {input.raw_R1} \
             --graft {input.fastq_graft_R1} --host {input.fastq_host_R1} \
              > {output.read_number_stat}
"""

rule stat_accuracy:
    input:
        raw_R1 = rules.merge_graft_host.output.merged_cvt_R1,
        fastq_graft_R1 = rules.convert_bam_to_fastq.output.fastq_graft_R1,
        fastq_host_R1 = rules.convert_bam_to_fastq.output.fastq_host_R1
    output:
        accuracy_stat = os.path.join(outdir, "accuracy_stat.txt")
    params: batch = "-l nodes=1:ppn=16,mem=64g"
    shell:
        """
python {stat_accuracy_path} {read_number} \
             {input.fastq_host_R1},host  {input.fastq_graft_R1},graft \
              > {output.accuracy_stat}
"""
