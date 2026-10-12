# Kallisto read-assignment rules.

rule build_kallisto_index:
    input:
        reference=get_reference_fasta
    output:
        index=READ_ASSIGNMENT_DIR + "/indexes/{assignment_reference}/kallisto/reference.idx"
    params:
        index_dir=READ_ASSIGNMENT_DIR + "/indexes/{assignment_reference}/kallisto",
        index_name="reference.idx",
    shell:
        """
        mkdir -p {params.index_dir:q}
        cd {params.index_dir:q}
        kallisto index -i {params.index_name:q} {input.reference:q}
        """


rule run_kallisto:
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2,
        index=rules.build_kallisto_index.output.index,
    output:
        abundance=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/{sample}/abundance.tsv",
        run_info=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/{sample}/run_info.json",
        bootstraps=temp(expand(READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kallisto/{{sample}}/bs_abundance_{bootstrap}.tsv", bootstrap=range(BOOTSTRAPS))),
        log=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/{sample}/kallisto.log",
    params:
        out_dir=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/{sample}",
        bootstraps=BOOTSTRAPS,
    threads: config["threads"]
    shell:
        """
        status=0
        kallisto quant -i {input.index:q} -o {params.out_dir:q} -b {params.bootstraps} --seed {config[seed]} -t {threads} --plaintext {input.r1:q} {input.r2:q} \
            2> {output.log:q} || status=$?

        if [ "$status" -ne 0 ]; then
            if grep -q "no reads pseudoaligned" {output.log:q}; then
                echo "[info] kallisto: no reads pseudoaligned for {wildcards.sample}; writing empty outputs."
                {PYTHON:q} -m cosmos_abundance.assignment.kallisto make-empty \
                    --abundance_output {output.abundance:q} \
                    --run_info_output {output.run_info:q} \
                    --bootstraps_output_dir {params.out_dir:q} \
                    --num_bootstraps {params.bootstraps} \
                    --read1 {input.r1:q}
            else
                echo "[error] kallisto failed for {wildcards.sample} with unexpected error:" >&2
                cat {output.log:q} >&2 || true
                exit "$status"
            fi
        fi
        """


rule combine_kallisto_bootstraps:
    input:
        expand(READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kallisto/{{sample}}/bs_abundance_{bootstrap}.tsv", bootstrap=range(BOOTSTRAPS))
    output:
        combined_bootstraps=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/{sample}/{sample}_kallisto_combined_bootstraps.tsv"
    shell:
        """
        {PYTHON:q} -m cosmos_abundance.assignment.kallisto combine-bootstraps \
            --input_files {input:q} \
            --output_file {output.combined_bootstraps:q}
        """


rule collate_kallisto_results:
    input:
        lambda wc: expand(
            READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kallisto/{sample}/{sample}_kallisto_combined_bootstraps.tsv",
            sample=get_samples(wc.sample_set),
        )
    output:
        READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/results/data/kallisto_results.all_samples.tsv"
    shell:
        """
        {PYTHON:q} -m cosmos_abundance.assignment.kallisto collate-results \
            --input_files {input:q} \
            --output_file {output:q}
        """


rule collate_kallisto_proportion_mapped:
    input:
        lambda wc: expand(
            READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kallisto/{sample}/run_info.json",
            sample=get_samples(wc.sample_set),
        )
    output:
        READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kallisto/results/data/kallisto_proportion_mapped.all_samples.tsv"
    shell:
        """
        {PYTHON:q} -m cosmos_abundance.assignment.kallisto proportion-mapped \
            --input_files {input:q} \
            --output_file {output:q}
        """
