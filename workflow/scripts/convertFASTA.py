"""
converts HMMER hit output file from Stockholm format to FASTA
"""
import os.path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

output_msa = sys.argv[1]
FASTA = sys.argv[2]

basename = os.path.basename(output_msa)
genome_gene = os.path.splitext(basename)[0]

records = []
msa = SeqIO.parse(output_msa, "stockholm")
for hit in msa:
    ungapped_seq = str(hit.seq.ungap("-"))
    ungapped_seq_upper = ungapped_seq.upper()
    record = SeqRecord(
        Seq(ungapped_seq_upper),
        id=hit.id,
        description=genome_gene
        )
    records.append(record)
SeqIO.write(records, FASTA, "fasta")
