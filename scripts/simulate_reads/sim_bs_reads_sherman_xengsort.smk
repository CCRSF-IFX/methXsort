
graft_fa = config["graft_fa"]
host_fa = config["host_fa"]
graft_name = os.path.basename(graft_fa).split(".")[0]
host_name = os.path.basename(host_fa).split(".")[0]
outdir = os.path.abspath(config["outdir"])
read_number = config["read_pair_number"]
xengsort_index_prefix = config.get("xengsort_index_prefix", None)

# Note: methXsort should be installed via pip (pip install methXsort)
# If not installed, you can still use: python /path/to/methXsort.py
sherman_path = config["sherman_path"]
stat_accuracy_path = config.get("stat_accuracy_path", None)

rule all:
    input:
        os.path.join(outdir, "ref_idx/xengsort_index.log"),
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
methXsort convert-reads --with_orig_seq \
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

rule xengsort_idx:
    input:
        graft_fa = rules.convert_ref.output.graft_fa_cvt,
        host_fa = rules.convert_ref.output.host_fa_cvt,
    output:
        log = os.path.join(outdir, "ref_idx/xengsort_index.log"),
    params:
        xengsort_index_prefix = os.path.join(outdir, "ref_idx/xengsort_index"),
        build = 1

    shell:
        """
methXsort xengsort-index \
    --host {input.host_fa} --graft {input.graft_fa} \
    --index {params.xengsort_index_prefix} \
    -n 7_000_000_000 --fill 0.88 --statistics summary -k 25 \
    > {output.log} 2>&1
"""

rule xengsort_classify:
    input:
        reads_R1 = rules.merge_graft_host.output.merged_cvt_R1,
        reads_R2 = rules.merge_graft_host.output.merged_cvt_R2,
    params:
        xengsort_out_prefix = os.path.join(outdir, "xengsort_classify/xengsort_classify"),
    output:
        classify_log = os.path.join(outdir, "classify.log"),
        fq_graft_R1 = os.path.join(outdir, "xengsort_classify/xengsort_classify-graft.1.fq.gz"),
        fq_graft_R2 = os.path.join(outdir, "xengsort_classify/xengsort_classify-graft.2.fq.gz"),
        fq_host_R1 = os.path.join(outdir, "xengsort_classify/xengsort_classify-host.1.fq.gz"),
        fq_host_R2 = os.path.join(outdir, "xengsort_classify/xengsort_classify-host.2.fq.gz"),
        #fq_ni_R1 = os.path.join(outdir, "xengsort_classify/xengsort_classify-neither.1.fq.gz"),
        #fq_ni_R2 = os.path.join(outdir, "xengsort_classify/xengsort_classify-neither.2.fq.gz"),
        #fq_ambiguous_R1 = os.path.join(outdir, "xengsort_classify/xengsort_classify-ambiguous.1.fq.gz"),
        #fq_ambiguous_R2 = os.path.join(outdir, "xengsort_classify/xengsort_classify-ambiguous.2.fq.gz"),
    shell:
        """
mkdir -p {outdir}/xengsort_classify && cd {outdir}/xengsort_classify && \
methXsort xengsort-classify \
    --read {input.reads_R1} --read2 {input.reads_R2} \
    --index {xengsort_index_prefix} \
    --out_prefix {params.xengsort_out_prefix} \
    --threads 33 \
    > {output.classify_log} 2>&1
"""

rule convert_restore_fastq:
    input:
        R1_graft = rules.xengsort_classify.output.fq_graft_R1,
        R2_graft = rules.xengsort_classify.output.fq_graft_R2,
        R1_host = rules.xengsort_classify.output.fq_host_R1,
        R2_host = rules.xengsort_classify.output.fq_host_R2,
    output:
        fastq_graft_R1 = os.path.join(outdir, "fastq/bbsplit_graft_R1.fastq.gz"),
        fastq_graft_R2 = os.path.join(outdir, "fastq/bbsplit_graft_R2.fastq.gz"),
        fastq_host_R1 = os.path.join(outdir, "fastq/bbsplit_host_R1.fastq.gz"),
        fastq_host_R2 = os.path.join(outdir, "fastq/bbsplit_host_R2.fastq.gz"),
    shell:
        """
methXsort restore-fastq \
    --read {input.R1_graft} --read2 {input.R2_graft} \
    --out {output.fastq_graft_R1} --out2 {output.fastq_graft_R2} 
methXsort restore-fastq \
    --read {input.R1_host} --read2 {input.R2_host} \
    --out {output.fastq_host_R1} --out2 {output.fastq_host_R2}
"""

rule stat_read_number: 
    input: 
        raw_R1 = rules.merge_graft_host.output.merged_cvt_R1,
        fastq_graft_R1 = rules.convert_restore_fastq.output.fastq_graft_R1,
        fastq_host_R1 = rules.convert_restore_fastq.output.fastq_host_R1
    output: 
        read_number_stat = os.path.join(outdir, "read_number_stat.txt")
    params: 
        batch = "-l nodes=1:ppn=16,mem=64g",
        xengsort_out_prefix = os.path.join(outdir, "xengsort_classify/xengsort_classify"),
    shell: 
        """
methXsort stat-split --raw {input.raw_R1} \
             --graft {input.fastq_graft_R1} --host {input.fastq_host_R1} \
             --xengsort_prefix {params.xengsort_out_prefix} \
             --sample_name "simulation_100k" \
              > {output.read_number_stat}
"""

rule stat_accuracy:
    input:
        raw_R1 = rules.merge_graft_host.output.merged_cvt_R1,
        fastq_graft_R1 = rules.convert_restore_fastq.output.fastq_graft_R1,
        fastq_host_R1 = rules.convert_restore_fastq.output.fastq_host_R1
    output:
        accuracy_stat = os.path.join(outdir, "accuracy_stat.txt")
    params: batch = "-l nodes=1:ppn=16,mem=64g"
    shell:
        """
python {stat_accuracy_path} {read_number} \
             {input.fastq_host_R1},host  {input.fastq_graft_R1},graft \
              > {output.accuracy_stat}
"""
