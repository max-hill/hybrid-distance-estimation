from typing import List
import math
import re
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def num_gaps(record: SeqRecord) -> int:
    return len(re.findall(r'[\-nX]', str(record.seq)))  # Counting gaps and 'n' or 'X' characters in the seq.


def fraction_max_gap(aln: List[SeqRecord]) -> float:
    """
    Returns the maximum fraction of gaps among the sequwnces in the alignment (used to filter chunks).
    """
    max_gaps = max((num_gaps(record) for record in aln))
    seq_len = len(aln[0].seq)
    return max_gaps / seq_len


def split_an_alignment(aln_path: str, schema='fasta', chunk_size=1000, max_gap=0.3) -> List[List[SeqRecord]]:
    aln: List[SeqRecord] = list(SeqIO.parse(aln_path, schema))
    seq_len = len(aln[0].seq)
    for record in aln:
        assert len(record.seq) == seq_len  # quick check that all sequences are of same length.

    chunks: List[List[SeqRecord]] = []
    chunk_num = math.ceil(seq_len / chunk_size)
    for chunk_i in range(chunk_num):
        left = chunk_i * chunk_size
        right = min((chunk_i + 1) * chunk_size, seq_len)  # right boundary is exclusive
        chunk = []
        for record in aln:
            subseq = SeqRecord(record.seq[left:right], id=record.id, description=record.description, name=record.name)
            chunk.append(subseq)
        if fraction_max_gap(chunk) <= max_gap:  # Make sure that the chunk doesn't have too many gaps.
            chunks.append(chunk)
    return chunks


def compute_distance(record1: SeqRecord, record2: SeqRecord) -> int:
    """
    Finds the number of different characters between the two aligned sequences. Ignores gaps.
    """
    dist = 0
    for i, char1 in enumerate(record1.seq):
        char2 = record2.seq[i]
        if char1 != '-' and char2 != '-':
            if char1 != char2:
                dist += 1
    return dist


def compute_jc_distances(chunks: List[List[SeqRecord]], seq1: str, seq2: str, bar_width=0.001):
    jc_distances = []
    distances = []
    for chunk in chunks:
        rec1 = [record for record in chunk if record.id == seq1][0]
        rec2 = [record for record in chunk if record.id == seq2][0]
        dist = compute_distance(rec1, rec2) / len(rec1.seq)
        distances.append(dist)
        jc_dist = -0.75 * (math.log(1 - (4 / 3) * dist))
        jc_distances.append(jc_dist)

    bars = int(max(distances) // bar_width)
    plt.hist(jc_distances, bins=bars, alpha=0.7)


if __name__ == '__main__':
    aln_path = '../data/set1c.fasta'
    chunks = split_an_alignment(aln_path, 'fasta', 1000, 0.3)
    print(f'Chunks created: {len(chunks)}')
    compute_jc_distances(chunks, 'B589', 'K22')
    compute_jc_distances(chunks, 'C33', 'C46')
    plt.show()
