'''
===================================================================================
Extract paired end UMI(5') from fastq, combine as its read name and merge the UMIs
===================================================================================

'''

import sys
import regex
import collections
from optparse import OptionParser
import os
import glob
import pandas as pd
import sys
import gzip
import itertools
from future.utils import iteritems

class Record:
  """A record representing a :term:`fastq` formatted record.
  Attributes
  ----------
  identifier : string Sequence identifier
  seq : string Sequence quals : string String representation of quality scores.
  format : string Quality score format. Can be one of ``sanger``, ``phred33``, ``phred64`` or ``solexa``.
  """

  def __init__(self, identifier, seq, quals, entry_format=None):
    self.identifier, self.seq, self.quals, entry_format = (
      identifier, seq, quals, entry_format)

  def guessFormat(self):
    '''return quality score format -
    might return several if ambiguous.'''
    RANGES = {
      'phred33': (33, 77),
      'solexa': (59, 106),
      'phred64': (64, 106),
    }
    c = [ord(x) for x in self.quals]
    mi, ma = min(c), max(c)
    r = []
    for entry_format, v in iteritems(RANGES):
      m1, m2 = v
      if mi >= m1 and ma < m2:
        r.append(entry_format)
    return r

  def __str__(self):
    return "@%s\n%s\n+\n%s" % (self.identifier, self.seq, self.quals)

class Extractor:
  ''' A functor which extracts barcodes from a read(s), filters the
  read(s) and updates the read(s). Keeps track of events in
  read_counts Counter
  '''
  def extract_5prime(self, sequence):
    return (sequence[:self.pattern_length],
            sequence[self.pattern_length:])
  def MergeRead(self,UMI1,UMI2,R1,R2):
    uim12_id = UMI1 + "_" + UMI2
    read12_seq = R1 + "_" + R2

    return uim12_id,read12_seq

  def __init__(self,pattern=None):
    self.pattern = pattern
    self.pattern_length = len(self.pattern)

  def __call__(self, read1, read2):
    umi1,new_read1=self.extract_5prime(read1.seq)
    umi2,new_read2=self.extract_5prime(read2.seq)
    uim12_id,read12_seq=self.MergeRead(umi1,umi2,new_read1,new_read2)
    return uim12_id,read12_seq



def fastqIterate(infile):
  '''iterate over contents of fastq file.'''

  def convert2string(b):
    if type(b) == str:
      return b
    else:
      return b.decode("utf-8")

  while 1:
    line1 = convert2string(infile.readline())
    if not line1:
      break
    if not line1.startswith('@'):
      print("parsing error: expected '@' in line %s" % line1)
    line2 = convert2string(infile.readline())
    line3 = convert2string(infile.readline())
    if not line3.startswith('+'):
      print("parsing error: expected '+' in line %s" % line3)
    line4 = convert2string(infile.readline())
    # incomplete entry
    if not line4:
      print("incomplete entry for %s" % line1)

    yield Record(line1[1:-1], line2[:-1], line4[:-1])

def EFread(readlist,Ecutoff,Scutoff):
  totalread = len(readlist)
  if totalread < Scutoff: return 'bad'
  if len(list(set(map(len, readlist)))) > 1: return 'bad'
  realread = ''
  for j in range(0,len(readlist[0])):
    bases = {}
    for read in readlist:
      base = read[j]
      if base in bases: bases[base] += 1
      else: bases[base] = 1
    check = 'bad'
    bs = ''
    for b in bases.keys():
      if float(bases[b])/float(totalread) >= Ecutoff:
        bs += b
    if len(bs) != 1: return 'bad'
    else:
      realread += bs
      check = 'good'
    if check == 'bad':
      return 'bad'
  return realread


def main():
  if len(sys.argv) !=5:
    sys.exit('[usage] python Dedup_UMI.py <folder> <pattern> <Ecutoff> <Scutoff>')
  folder = sys.argv[1]
  pattern = sys.argv[2]
  Ecutoff = float(sys.argv[3])
  Scutoff = int(sys.argv[4])
  infiles=glob.glob(folder+'/*_L001_R1_001.fastq.gz')
  #check number of files
  filenames=[]
  for file in infiles:
    filename=os.path.basename(file).rsplit('_L001_')[0]
    filenames.append(filename)
  print('there are total %s pared end files \n' % len(filenames))
  # set up read extractor
  #Open decompose tmp files
  ReadExtractor=Extractor(pattern=pattern)
  for n in filenames:
    outfile=open('result/'+n+'_tagfree.fa', 'w')
    print('working on...', n, 'and its mate')
    r1file=gzip.open(folder+'/'+n+'_L001_R1_001.fastq.gz','rt',encoding="utf-8")
    r2file=gzip.open(folder+'/'+n+'_L001_R2_001.fastq.gz','rt',encoding="utf-8")
    #read1s and read2s are generator for loop fastq reads
    read1s=fastqIterate(r1file)
    read2s=fastqIterate(r2file)
    progCount = 0
    displayMax = 1000000
    Rclusters = {}
    for read1 in read1s:
      read2=next(read2s)
      progCount += 1
      if progCount % displayMax == 0:
        print('Processed %s Million reads'% (progCount/displayMax))
      #extract UMI1_UMI2 as ID, read1_read2 as read
      ID,read= ReadExtractor(read1,read2)
      #filter out neither read or ID containing N
      if 'N' in read or 'N' in ID: continue

      #collect reads by barcode while iterating over reads Rclusters{ 'ID1':[read1,read2...], 'ID2':[read1,read2...]}
      if ID in Rclusters: Rclusters[ID].append(read)
      else: Rclusters[ID] = [read]

    #return barcode and its numbers statistics
    numbarcode={}
    for i in Rclusters.keys():
      num=len(Rclusters[i])
      if num in numbarcode:
        numbarcode[num]+=1
      else:
        numbarcode[num] = 1
    readsperbcstats = pd.DataFrame(sorted(numbarcode.items()),
                                   columns=['number of reads', 'number of barcodes']).set_index('number of reads')
    print(readsperbcstats.to_string())

    #loop over Rclusters and return unique ID and its reads
    for (iID,(ID, reads)) in enumerate(Rclusters.items()):
      if (iID + 1) % 2e5 == 0:
        print("Barcodes examined so far: {0}".format(
          iID + 1))
      efreeread = EFread(reads, Ecutoff, Scutoff)
      if efreeread != 'bad':
        r1 = efreeread.rsplit('_')[0]
        r2 = efreeread.rsplit('_')[1]
        outfile.write('>' + ID + "\n" + r1 + "\n")
        outfile.write('>' + ID + "\n" + r2 + "\n")
    outfile.close()



if __name__ == '__main__':
  main()

