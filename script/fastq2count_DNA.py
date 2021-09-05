#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import Counter

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "---":"-"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def sum_mut(aa1,aa2):
    return sum ( aa1[i] != aa2[i] for i in range(len(aa1)) )


def call_mutid(mutpep,refseq,shift,amplicon):
  mut_id_ls = []
  assert (len(mutpep) == len(refseq))
  for n in range(len(mutpep)):
    pos = n+shift
    if refseq[n]!=mutpep[n]:
       mut_id_ls.append(amplicon+'|'+refseq[n]+str(pos)+mutpep[n])
  return mut_id_ls

def Call_Amp_ID(record_id):
    amp=str(record_id)[-4:]
    return amp


def cal_fastq_dic(fastq,ref):

    print ("reading %s" % fastq)

    Rrecords = SeqIO.parse(fastq,"fastq")
    ref = str(ref.seq)
    mut_id_ls = []
    error_read=0
    len_error=0
    shift=0
    amp_ref=''
    for record in Rrecords:
        amp= Call_Amp_ID(record.id)
        if amp == 'Amp1':
            amp_ref = ref[51:375] #amplicon 1 start and end
            shift=52
        elif amp == 'Amp2':
            amp_ref = ref[375:699] #amplicon 2 start and end
            shift=376
        elif amp == 'Amp3':
            amp_ref = ref[699:1026]  #amplicon 3 start and end
            shift=700
        else:
            error_read +=1
        seq = str(record.seq)
        if len(seq) != len(amp_ref):
            len_error+=1
        if len(seq) == len(amp_ref):
            mut = seq
            if sum_mut(mut, amp_ref) == 0:

                mut_id = amp+'|'+'WT'
                mut_id_ls.append(mut_id)
            else:
                call_mut_id = call_mutid(mut, amp_ref,shift,amp)
                mut_id = "-".join(sorted(call_mut_id, key=lambda x:int(x[6:-1])))
                mut_id_ls.append(mut_id)
    print('no amplicon tag, a total of %s reads were exclued ' %error_read)
    print('deletion or insertion, a total of %s reads were excluded' %len_error)
    dna_dict = Counter(mut_id_ls)
    return dna_dict

def write_mut_table(mut_dic,outfilename):
  outfile = open(outfilename,'w')
  outfile.write("\t".join(['Mutation', 'count'])+"\n")
  for mut in mut_dic.keys():
    Mutation = mut
    count = mut_dic[mut]
    outfile.write("\t".join(map(str,[Mutation, count]))+"\n")
  outfile.close()
  print('Written %s' %outfilename)


def main():
    if len(sys.argv) != 4:
        sys.exit('[usage] python fastq2count.py <fastq file> < reference> <mutation table filename>')
    fastqfile = sys.argv[1]
    ref = sys.argv[2]
    outfilename = sys.argv[3]
    ref = SeqIO.read(ref,"fasta")
    mutation_dic = cal_fastq_dic(fastqfile,ref)
    write_mut_table(mutation_dic, outfilename)


if __name__ == "__main__":
  main()
