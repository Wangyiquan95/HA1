#!/usr/bin/python
#APRIL 5, 2014 BY NICHOLAS C. WU
#TESTING EXAMPLE: python script/Fastq2ErrorFreeFasta.py -i fastq/ -o otest -F _R1_ -R _R2_ -d 4 -b 1-3 -p 4-12 -e 0.9 -s 3
#Note: This script is for paired-end reads 
#Note: Barcode and tag should be at the 5' end before the reads
#
#EAMPLE for barcode and tag parameters:
#B stands for population ID (Barcode); N stands for error-correction tag
#Example 1
#  Adaptor sequence: BBBNNNNNNNNN
#  Then, barcode and tag parameters should be: -b 1-3 -p 4-12
#
#Example 2
#  Adaptor sequence: NNNNNBBB
#  Then, barcode and tag parameters should be: -b 6-8 -p 1-5
#
#Example for identifier parameters:
#
#####################################################################
#############################   ENJOY   #############################
#####################################################################
import os
import glob
import sys
import getopt
import subprocess
import gzip
import itertools
from optparse import OptionParser
from Bio import SeqIO

#SUBROUTINE 1: DECOMPOSE READS TO MINMIZE MEMORY USAGE
def formatread(countr1,Decomp,tagstart,tagend,barstart,barend,tmp,read1i,read2i):
  #Open decompose tmp files
  prefixcomb = list(itertools.product(['A','C','T','G'],repeat=Decomp))
  dprefixs    = {}
  for prefix in prefixcomb:
    prefix = ''.join(prefix)
    dprefixs[prefix] = open(tmp+prefix,'w')
  #Demultiplex by tmp files
  for r1file in sorted(countr1):
    print('working on...', r1file, 'and its mate')
    r2file  = r1file.replace(read1i,read2i)
    r1file  = gzip.open(r1file,'rt',encoding='utf-8')
    r2file  = gzip.open(r2file,'rt',encoding='utf-8')
    countline  = 0
    for r1rec in r1file.readlines():
      r2rec = next(r2file)
      r1rec = r1rec.rstrip()
      r2rec = r2rec.rstrip()
      countline += 1
      if countline == 4: countline = 0
      elif countline%2 == 0:
        #filter out small read
        #if len(r1rec) < 200 or len(r2rec) <200: continue
        bar = r1rec[barstart:barend]+'-'+r2rec[barstart:barend]
        tag = r1rec[tagstart:tagend]+'-'+r2rec[tagstart:tagend]
        tagdecomp = tag[0:Decomp]
        readstart = max([barend,tagend])
        if 'N' not in r1rec and 'N' not in r2rec:
          dprefixs[tagdecomp].write(tag+"\t"+bar+"\t"+r1rec[tagend::]+"\t"+r2rec[tagend::]+"\n")
  #Close all decompose tmp files
  for dprefix in dprefixs.keys():
    dprefixs[dprefix].close()

#SUBROUTINE 2a: ERROR CORRECTION
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

#SUBROUTINE 2b: ERROR CORRECTION
def errorcorrect(tmp,ofile,Ecutoff,Scutoff):
  tmpfiles = sorted(glob.glob(tmp+'*'))
  countfile = 0
  ofile = open(ofile,'w')
  for tmpfile in tmpfiles:
    countfile += 1
    print('WORKING ON FILE:', countfile)
    tmpfile = open(tmpfile,'r')
    Rclusters = {}
    for line in tmpfile:
      line = line.rstrip().rsplit("\t")
      tag  = line[0]
      bar  = line[1]
      r1   = line[2]
      r2   = line[3]
      ID   = tag+'_'+bar
      read = r1+'_'+r2
      if 'N' in read or 'N' in ID: continue
      #filter out small read
      if len(r1) <200 or len(r2) <200: continue
      if ID in Rclusters: Rclusters[ID].append(read)
      else: Rclusters[ID] = [read]
    tmpfile.close()
    for ID in Rclusters.keys():     
      efreeread = EFread(Rclusters[ID],Ecutoff,Scutoff)
      if efreeread != 'bad':
        r1 = efreeread.rsplit('_')[0]
        r2 = efreeread.rsplit('_')[1]
        ofile.write('>'+ID+'_F'+"\n"+r1+"\n")
        ofile.write('>'+ID+'_R'+"\n"+r2+"\n")
  ofile.close()

#PARSING ARGUMENTS
def ParseArg():
  usage  = 'python script/Fastq2ErrorFreeFasta.py -i FOLDER -o DESTINATION FILE -b POSITION RANGE -p POSITION RANGE -F STRING -R STRING -d INTEGER -e FLOAT -s INTEGER'
  desc   = 'Note: This script is for paired-end reads and Barcode and tag should be at the 5\' end before the reads. At this current version, all parameters should be given. A larger the decomposer parameter will decrease memory usage. But too high will cause problem in python. Recommend using decomposer = 3 or 4'
  parser = OptionParser(description=desc, usage=usage)
  parser.add_option("-i",help="folder containing all fastq (.fq.gz or fastq.gz) files (e.g. fastq/)", metavar="FOLDER")
  parser.add_option("-o",help="output error-free read fasta files", metavar="FILE")
  parser.add_option("-b",help="position of the barcode/population ID (e.g. 1-10), 0-0 if no barcode", metavar="POSITION RANGE")
  parser.add_option("-p",help="position of the tag (e.g. 1-10)", metavar="POSITION RANGE")
  parser.add_option("-d",help="decomposer during tag demultiplexing, minimize memory usage", metavar="INTEGER")
  parser.add_option("-F",help="identifier (replaceable with argument in -F) for read 1 files (e.g. _R1_)", metavar="STRING")
  parser.add_option("-R",help="identifier (replaceable with argument in -R) for read 2 files (e.g. _R2_)", metavar="STRING")
  parser.add_option("-e",help="minimum fraction of # of reads within a read cluster without an error at a bp", metavar="FLOAT")
  parser.add_option("-s",help="minimum size of a read cluster", metavar="INTEGER")
  (options, args) = parser.parse_args()

  folder = options.i
  fafile = options.o
  barpos = options.b
  tagpos = options.p
  read1i = options.F
  read2i = options.R
  Decomp = options.d
  Ecutoff = options.e
  Scutoff = options.s
  errorstate = '0'
  if not folder: print('INPUT ERROR: input folder missing'); errorstate = '1'
  if not fafile: print('INPUT ERROR: output fasta file missing'); errorstate = '1'
  if not barpos: print('INPUT ERROR: barcode position missing'); errorstate = '1'
  if not tagpos: print('INPUT ERROR: tag position missing'); errorstate = '1'
  if not read1i: print('INPUT ERROR: read 1 identifier missing'); errorstate = '1'
  if not read2i: print('INPUT ERROR: read 2 identifier missing'); errorstate = '1'
  if not Decomp: print('INPUT ERROR: decomposer missing'); errorstate = '1'
  if not Ecutoff: print('INPUT ERROR: decomposer missing'); errorstate = '1'
  if not Scutoff: print('INPUT ERROR: decomposer missing'); errorstate = '1'
  if errorstate == '1': sys.exit()

  #CHECKING THE FORMAT OF ALL INPUT VARIABLES
  errorstate = '0'
  #CHECKING TAG/BARCODE POSITION 1
  tag = tagpos.rsplit('-')
  bar = barpos.rsplit('-')
  if len(tag) != 2: print('FORMAT ERROR: position of the tag (e.g. 1-10)'); errorstate = '1'
  try: tagstart = int(tag[0])-1
  except ValueError: print('VALUE ERROR: start position of the tag should be an integer'); errorstate = '1'
  try: tagend = int(tag[1])
  except ValueError: print('VALUE ERROR: end position of the tag should be an integer'); errorstate = '1'

  if len(bar) != 2: print('FORMAT ERROR: position of the tag (e.g. 1-10)'); errorstate = '1'
  try: barstart = int(bar[0])-1
  except ValueError: print('VALUE ERROR: start position of the barcode should be an integer'); errorstate = '1'
  try: barend = int(bar[1])
  except ValueError: print('VALUE ERROR: end position of the barcode should be an integer'); errorstate = '1'
  
  try: Decomp = int(Decomp)
  except ValueError: print('VALUE ERROR: decomposer should be an integer'); errorstate = '1'
  
  try: Ecutoff = float(Ecutoff)
  except ValueError: print('VALUE ERROR: Error cutoff should be a float'); errorstate = '1'
  if Ecutoff < 0 or Ecutoff > 1: print('Error cutoff should be within 0 to 1'); errorstate = '1'

  try: Scutoff = int(Scutoff)
  except ValueError: print('VALUE ERROR: Cluster size cutoff should be an integer'); errorstate = '1'
  if Scutoff < 2: print('VALUE ERROR: Cluster size cutoff should be >= 2'); errorstate = '1'

  if errorstate == '1': sys.exit()
  
  #CHECKING TAG/BARCODE POSITION 2
  if barend < barstart: print('VALUE ERROR: for the barcode, start position should be smaller than end position'); errorstate = '1'
  if tagend < tagstart: print('VALUE ERROR: for the tag, start position should be smaller or equal (if no barcode) than end position'); errorstate = '1'
  if errorstate == '1': sys.exit()

  #CHECKING EXISTENCE OF FASTQ FILES
  filenames = glob.glob(folder+'/*.fq.gz')
  filenames.extend(glob.glob(folder+'/*.fastq.gz'))
  if len(filenames) == 0: print('INPUT ERROR: no fastq files are found (make sure they have fq.gz or fastq.gz as suffix)'); sys.exit()
  countr1 = []
  countr2 = []
  for filename in filenames:
    if read1i in filename: countr1.append(filename)
    if read2i in filename: countr2.append(filename)
  if len(countr1)   != len(countr2): print('INPUT ERROR: number of read 1 files is different from number of read 2 files'); sys.exit()
  if len(list(set(countr1).intersection(set(countr2)))) != 0: print('INPUT ERROR: bad read identifier'); sys.exit()
  return Ecutoff,Scutoff,countr1,Decomp,tagstart,tagend,barstart,barend,read1i,read2i,fafile,folder

#CREATE TMP FOLDER
def creattmpfolder(folder):  
  i = 1
  while i: 
    i += 1
    tmp = folder + '/tmp_'+str(i)+'/'
    tmps = glob.glob(tmp)
    if len(tmps) == 0: 
      os.system('mkdir '+tmp)
      break
  print('Created a tmp folder: %s' % tmp)
  return tmp

#CLEAN UP TMP FILES AFTER ANAYLSIS
def cleaning(tmp,fafile):
  os.system('rm -r '+tmp)
  print('tmp folder %s removed' % tmp)
  print('ERROR FREE READS ARE PRINTED IN FILE: %s AS FASTA FORMAT' % fafile)
  print('READ ID IS FORMATED AS >FORWARDTAG-REVERSETAG_FORWARDBARCODE-REVERSEBARCODE_READDIRECTION')
  print('READDIRECTION - F: FORWARD; R: REVERSE')
  
#MAIN#
def main():
  Ecutoff,Scutoff,countr1,Decomp,tagstart,tagend,barstart,barend,read1i,read2i,fafile,folder = ParseArg()
  tmp = creattmpfolder(folder)
  ########################################################################
  #           Step 1: Decompose the reads according to the tag sequence  #
  ########################################################################
  formatread(countr1,Decomp,tagstart,tagend,barstart,barend,tmp,read1i,read2i)
  print('DONE DECOMPOSING')
  print('TOTAL OF:', len(glob.glob(tmp+'*')), 'TMP READ FILES GENERATED')

  ########################################################################
  #           Step 2: Error correction                                   #
  ########################################################################
  print('START ERROR CORRECTION')
  errorcorrect(tmp,fafile,Ecutoff,Scutoff)
  print('DONE ERROR CORRECTION')
  ########################################################################
  #           Step 3: Clean up                                           #
  ########################################################################
  cleaning(tmp,fafile) 

if __name__ == '__main__':
  main()
