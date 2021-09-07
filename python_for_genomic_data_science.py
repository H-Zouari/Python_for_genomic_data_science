
!pip install Biopython
!pip install pandas

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pandas as pd
import time
from scipy.stats import mode
import re
from collections import Counter


# First we parse the records in the fasta file using the SeqIO.parse method in the Biopython library
# It detects identifiers and sequences and generates and seq record object containing them
 
ORFs = []
ids = []
for dna_record in SeqIO.parse("/content/drive/MyDrive/Colab_Notebooks/Coursera/dna2.fasta", 'fasta'):
# the seq object containing the DNA sequence as a string is extracted and added into a new list
  ids.append(str(dna_record.id))
  for frame in range(3):
      trans= str(dna_record.seq[frame:].translate(to_stop=False))
      RFs = list((m.span() for m in re.finditer("M\D*?\*", trans)))
      ORFs.extend([RFs])


def Longest_ORFs(ORFs):
  longest = []
  seqinfo = defaultdict(list)
  for i in range(0,len(ORFs)):
      if len(ORFs[i]) !=0: 
        longest.append(max(ORFs[i], key = lambda sub: abs(sub[1] - sub[0])))
      else:
        longest.append((0,0))

  for x in range(3,len(longest),3):
    seqinfo[ids[int(x/3)-1]] = longest[x-3:x]
  return seqinfo

ORF_Table = pd.DataFrame(Longest_ORFs(ORFs))
ORF_lengths = ORF_Table.applymap(lambda sub: abs(sub[1] - sub[0])*3)
ORF_start = ORF_Table.applymap(lambda start: start[0]*3)

sequences = {}
lengths = []
for record in SeqIO.parse("/content/drive/MyDrive/Colab_Notebooks/Coursera/dna2.fasta", "fasta"):
  sequences[record.id] = len(record.seq)
  lengths.append(len(record.seq))
shortest = min(sequences.values())
longest = max(sequences.values())

Shortest_seqs = [key  for (key, value) in sequences.items() if value == shortest]
Longest_seqs = [key  for (key, value) in sequences.items() if value == longest]

print(len(sequences))
print(f"shortest sequence length: {shortest}")
print(f"longest sequence length: {longest}")
print(f"list of shortest sequences: {Shortest_seqs} with length {shortest}")
print(f"list of longest sequences : {Longest_seqs} with length {longest}")

seqs = ''
for dna_record in SeqIO.parse("/content/drive/MyDrive/Colab_Notebooks/Coursera/dna2.fasta", 'fasta'):
  seqs += str(dna_record.seq)

""" Here i created a wrapper function which calculates execution time to then decorate the two 
functions for finding repeats and compare their perforemance """

def timing(func):
    def wrapper(*arg, **kw):
        start = time.time()
        result = func(*arg, **kw)
        end = time.time()
        return (end  - start), result, func.__name__
    return wrapper
 
@timing
def find_repeats(seqs):
  hits = re.findall(r'(?=(([ATCG]){12}))', seqs)
  hits = [hit[0] for hit in hits]
  freq_motif = mode(hits)
  return freq_motif

print(find_repeats(seqs))
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""" this solution comes from https://github.com/burun/Python-for-Genomic-Data-Science-2015-Coursera/blob/master/Final%20exam/question_8_10.py slightly modified, i used it
to compare perforemance """
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sequences = []
for dna_record in SeqIO.parse("/content/drive/MyDrive/Colab_Notebooks/Coursera/dna2.fasta", 'fasta'):
  sequences.append(dna_record.seq)


def all_repeats(sequence):
    length = len(sequence)
    repeats = []
    for i in range(length):
        repeats.append(sequence[i:i + 12 ])
    return repeats

@timing
def most_frequent(List):
  l_x_repeats = []
  for i in sequences:
      repeats_list = all_repeats(i)
      for j in repeats_list:
          l_x_repeats.append(j)
  occurence_count = Counter(List)
  repeat = occurence_count.most_common(1)[0][0]
  freq = l_x_repeats.count(repeat)
  return repeat , freq

print(most_frequent(l_x_repeats))
