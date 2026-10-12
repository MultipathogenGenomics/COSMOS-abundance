# KMA read-assignment rules.

rule build_kma_index:
    input:
        reference=get_reference_fasta
    output:
        multiext(
            READ_ASSIGNMENT_DIR + "/indexes/{assignment_reference}/kma/reference.kma_index",
            ".comp.b",
            ".length.b",
            ".name",
            ".seq.b",
        )
    params:
        index_prefix=READ_ASSIGNMENT_DIR + "/indexes/{assignment_reference}/kma/reference.kma_index"
    threads: config["threads"]
    shell:
        """
        kma_index -i {input.reference:q} -o {params.index_prefix:q}
        """


rule run_kma:
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2,
        index=rules.build_kma_index.output,
    output:
        default=multiext(
            READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kma/{sample}/{sample}_kma",
            ".mapstat",
            ".res",
            ".aln",
            ".fsa",
            ".frag.gz",
        ),
        vcf=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kma/{sample}/{sample}_kma.vcf.gz",
        sam=temp(READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kma/{sample}/{sample}_kma.sam"),
    params:
        out_prefix=READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kma/{sample}/{sample}_kma",
        index_prefix=READ_ASSIGNMENT_DIR + "/indexes/{assignment_reference}/kma/reference.kma_index",
        kma_options="-ef -ex_mode -sam 4 -vcf 2",
    threads: config["threads"]
    shell:
        """
        kma -ipe {input.r1:q} {input.r2:q} -o {params.out_prefix:q} -t_db {params.index_prefix:q} -t {threads} \
            {params.kma_options} > {output.sam:q}
        """


rule collate_kma_results:
    input:
        res=lambda wc: expand(
            READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kma/{sample}/{sample}_kma.res",
            sample=get_samples(wc.sample_set),
        ),
        mapstat=lambda wc: expand(
            READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kma/{sample}/{sample}_kma.mapstat",
            sample=get_samples(wc.sample_set),
        ),
    output:
        READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kma/results/data/kma_results.all_samples.tsv"
    shell:
        """
        {PYTHON:q} -m cosmos_abundance.assignment.kma collate-results \
            --input_files {input.res:q} \
            --output_file {output:q}
        """


rule collate_kma_proportion_mapped:
    input:
        mapstat=lambda wc: expand(
            READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kma/{sample}/{sample}_kma.mapstat",
            sample=get_samples(wc.sample_set),
        ),
        fragments=lambda wc: expand(
            READ_ASSIGNMENT_DIR + "/{{assignment_reference}}/{{sample_set}}/kma/{sample}/{sample}_kma.frag.gz",
            sample=get_samples(wc.sample_set),
        ),
    output:
        READ_ASSIGNMENT_DIR + "/{assignment_reference}/{sample_set}/kma/results/data/kma_proportion_mapped.all_samples.tsv"
    shell:
        """
        {PYTHON:q} -m cosmos_abundance.assignment.kma proportion-mapped \
            --input_files {input.mapstat:q} \
            --output_file {output:q}
        """
