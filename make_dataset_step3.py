import json
import gzip
import re
import sys 

#load calls into a dictionary of calls and sequences into a dictionary (ID is the gene name)

sequences = {}
calls = {}

with open('Arabidopsis_thaliana.TAIR10.52.json') as f:
   for line in f:
      doc = json.loads(line)
      calls[doc['id']] = doc['calls']

#create a list of genes with only 2 isoforms
genes = []
with open("test2.txt") as f:
   for line in f:
      genes.append(line.replace('\n', ''))

def rc(seq):
   replace = { "A":"T", "C":"G", "G":"C", "T":"A", "N":"N",
               "a":"t", "c":"g", "g":"c", "t":"a", "n":"n" }
   return "".join(replace.get(x, "N") for x in seq[::-1])

def dict_from_fasta_file(f):
   ''' 
   Read the whole file in memory. This can crash if the file is big.
   '''
   genome = dict()
   txt = f.read() # Read it all.
   # Split on fasta header
   for item in txt.decode("ascii").split('>'):
      # Separate header from newline-delimited
      # sequence (proseq) and join the sequence.
      if item == "": 
         continue
      header, proseq = item.split('\n', 1)
      seq = proseq.replace('\n', '') 
      # Keep only the first word of the header.
      genome[re.sub("\s.*", "", header)] = seq
   genome["1"] = genome["NC_003070.9"]
   genome["5"] = genome["NC_003076.8"]
   genome["3"] = genome["NC_003074.8"]
   genome["2"] = genome["NC_003071.7"]    
   genome["4"] = genome["NC_003075.7"]
   genome["Mt"] = genome["NC_037304.1"]
   genome["Pt"] = genome["NC_000932.1"]
   del genome["NC_003070.9"]
   del genome["NC_003076.8"]
   del genome["NC_003074.8"]
   del genome["NC_003071.7"]
   del genome["NC_003075.7"]
   del genome["NC_037304.1"]
   del genome["NC_000932.1"]

   return genome

def getExons(genome, strand, chrom, start, end, location, length):
   exon_p = [int(i) for i in location.split(",")]             
   exon_l = [int(i) for i in length.split(",")]
   sequences = []
   for i, (p_, l_) in enumerate(zip(exon_p,exon_l)):
         if strand == "-":
            sequence = rc(genome[chrom][start-1:end][p_:p_+l_])
            sequences.append(sequence)
         elif strand == "+":
            sequence = genome[chrom][start-1:end][p_:p_+ l_]
            sequences.append(sequence)
         else:
            print("strand is not positive or negative")

   return sequences

def getIntrons(genome, strand, chrom, start, end, location, length):
   exon_p = [int(i) for i in location.split(",")]    
   exon_l = [int(i) for i in length.split(",")]
   sequences = []
   for i, (p_, l_) in enumerate(zip(exon_p,exon_l)):
      if i != (len(exon_p)-1):
         if strand == "-":
            sequence = rc(genome[chrom][start-1:end][p_+l_:exon_p[i+1]])
            sequences.append(sequence)
         elif strand == "+":
            sequence = genome[chrom][start-1:end][p_+l_:exon_p[i+1]] 
            sequences.append(sequence)
         else:
            print("strand is not positive or negative")
      else:
         break
   return sequences

def getDifference(sequencesA, SequencesB):
   difference = list(set(sequencesA)-set(sequencesB))
   difference2 = list(set(sequencesB)-set(sequencesA))
   [difference.append(x) for x in difference2]
   return difference

def getJunctions(exons, introns):
   junctionA = [ "".join([a, b]) for a, b in zip(exons, introns)]
   junctionB = ["".join([introns[i], exons[i+1]]) for i in range(len(exons)-2)] 
   return junctionA, junctionB


if __name__ == '__main__':
   output_dictionary = {"id": 0}
   fastagz  = sys.argv[1]
   coordgz = sys.argv[2]
   with gzip.open(fastagz) as f:
      genome = dict_from_fasta_file(f)
   with open(coordgz) as f:
      for line in f:
         items = line.split()
         name = items[0]
         chrom = items[1]
         gene = name.split(".")[0]
         start = int(items[2])
         end = int(items[3])
         strand = items[4]
         exon_p = items[5]
         exon_l = items[6]
         exons_list = []
         
         if gene in genes and genome[chrom]!= None:
            exons = getExons(genome, strand, chrom, start, end, exon_p, exon_l)
            introns = getIntrons(genome, strand, chrom, start, end, exon_p, exon_l)
            junctions = getJunctions(exons, introns)
            if output_dictionary["id"] == gene:
                  sequencesB = output_dictionary["exon_intron"]
                  difference_exons = getDifference(junctions[0], sequencesB)
                  for num, j in enumerate(difference_exons):
                     print(json.dumps({"id": gene + "-exon_intron." + str(num), "seq": j,"    labels": 1}))
                  sequencesC = output_dictionary["intron_exon"]
                  difference_introns = getDifference(junctions[1], sequencesC)
                  for num, j in enumerate(difference_introns):
                     print(json.dumps({"id": gene + "-intron_exon." + str(num), "seq": j,"labels": 1}))

               #sequencesB = output_dictionary["intron"]
               #sequencesC = output_dictionary["exon"]
               #difference_introns = getDifference(introns, sequencesB)
               #difference_exons = getDifference(exons, sequencesC)
               #print(json.dumps({"id": gene, "intron": difference_introns,"exon":difference_exons, "labels": 1}))
            else:
                  output_dictionary = {"id": gene, "exon_intron": junctions[0], "intron_exon": junctions[1]}
         else: continue








