# only variable needed to change


PROJECT_PATH='/Users/yiquan/PycharmProjects/HA1'
FASTA = PROJECT_PATH + '/result/tagfree_fasta/{SAMPLENAME}_tagfree.fa'
SAMPLENAMES, = glob_wildcards(FASTA)
REF=PROJECT_PATH + '/ref/PR8HA_ORF.fa'
#print(SAMPLENAMES)
FAS=PROJECT_PATH + '/result/tagfree_fasta/{SAMPLENAME}_tagfree.fa'
RESULT_PATH = PROJECT_PATH + '/result/{SAMPLENAME}'
RM_PRIMER_FAS = RESULT_PATH + '/{SAMPLENAME}_rm_primer.fa'
UNTRIM_FAS = RESULT_PATH + '/{SAMPLENAME}_untrim.fa'
RM_PRIMER_FQ = RESULT_PATH + '/{SAMPLENAME}_rm_primer.fq'
ASSEMBLED_FQ = RESULT_PATH + '/{SAMPLENAME}_assembled.fq'
TABLE = RESULT_PATH + '/{SAMPLENAME}_count.tsv'

rule all:
    input:
        expand(TABLE, SAMPLENAME=SAMPLENAMES),
        expand(RM_PRIMER_FAS, SAMPLENAME = SAMPLENAMES),
        expand(RM_PRIMER_FQ, SAMPLENAME = SAMPLENAMES),
        expand(ASSEMBLED_FQ, SAMPLENAME = SAMPLENAMES)
# remove primer(only for paired-primer are removed) &filter out small reads(-m 100:100(R1:R2))
rule rm_primer:
    input:
        FAS
    params:
        rename=lambda wc: "'{id}_{adapter_name}'"
    output:
        trim_o=RM_PRIMER_FAS,
        untrim_o = UNTRIM_FAS
    shell:
        '''cutadapt --pair-adapters --interleaved '''\
        '''-g "Amp1=CACTTGCAGCTGCAGATGCA;e=0.2" -g "Amp2=GGGAGCAATTGAGCTCAGTG;e=0.2" -g "Amp3=CCCCGGAAATAGCAGAAAGA;e=0.2" '''\
        '''-G "Amp1=CGAATCTTTCGAATGATGA;e=0.2" -G "Amp2=CTTGATCTCTTACTTTGGG;e=0.2" -G "Amp3=TGGCTCCAAATAGACCTCT;e=0.2" '''\
        '''--rename={params} '''\
        '''--untrimmed-output {output.untrim_o} '''
        '''-O 10 -m 100:100 -o {output.trim_o} {input}'''

rule fas2fq:
    input:
        RM_PRIMER_FAS
    output:
        RM_PRIMER_FQ
    shell:
        "seqtk seq -F '#' {input} > {output}"
#merge paired-end reads
rule flash:
    input:
        RM_PRIMER_FQ
    output:
        ASSEMBLED_FQ
    shell:
        "flash -m 100 -M 200 -I {input} -c > {output} "

rule fq2count:
    input:
        ASSEMBLED_FQ
    params:
        REF_FA=REF
    output:
        TABLE
    shell:
        'python fastq2count.py {input} {params} {output}'
