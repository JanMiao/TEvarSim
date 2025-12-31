import numpy as np
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------- sampling TE insertions with min distance ----------------
def sample_TEins(start: int, end: int, n: int, TEdistance: int):
    """
    sample_TEins: Sample n TE insertion positions between start and end with a minimum distance of TEdistance.
    return: a list of positions
    """
    region_len = end - start
    # strict feasibility check
    min_required = n * TEdistance
    if min_required >= region_len:
        raise ValueError(
            f"TEdistance too large. Please reduce TEdistance or the number of TEs to be simulated.")

    # compressed sampling
    max_start = end - min_required
    base = np.sort(
        np.random.choice(
            np.arange(start, max_start),
            size=n,
            replace=False
        )
    )

    # expand by offsets
    offsets = np.arange(n) * TEdistance
    positions = base + offsets

    return positions

# ---------------- adding background SVs ----------------
def bgSV(bedin:str, bedout:str, nSV:int, ins_ratio:float, fasta_in:str, fasta_out:str):
    """
    Add background SVs (INS/DEL) to a BED file.
    Args:
        bedin: input BED file (only TE insertions/deletions)
        bedout: output BED file with added SVs
        nSV: number of SVs to add
        ins_ratio: ratio of insertions among the background SVs
        fasta_file: FASTA file of sequences with newly added insertions (INS)
    """
    # parse the input bed file
    TEs = []
    with open(bedin, "r") as fin:
        for line in fin:
            fields = line.strip().split("\t")
            chrom, start, end, teID = fields[0], int(fields[1]), int(fields[2]), fields[3]
            TEs.append((chrom, start, end, teID))

    # number of background SVs
    nINS = int(nSV * ins_ratio)
    nDEL = nSV - nINS
    
    # generating DEL and INS randomly
    SVmin, SVmax = 30, 300
    # DELlens = np.random.uniform(low=SVmin, high=SVmax, size=nDEL)
    INSlens = np.random.randint(low=SVmin, high=SVmax, size=nINS)
    # INS sequences
    bgINS_seqs = []
    for idx, ilen in enumerate(INSlens):
        seq = ''.join(random.choices('ATGC', k=ilen))
        bgINS_seqs.append((f"bgINS_{idx}_{ilen}", seq))
    
    # output all INS sequences to new fasta file
    records = list(SeqIO.parse(fasta_in, "fasta"))
    new_records = [SeqRecord(Seq(i), id=j) for j, i in bgINS_seqs]
    records.extend(new_records)
    SeqIO.write(records, fasta_out, "fasta")
    
    # intervals for background SV
    existing_intervals = [(i[1], i[2]) for i in TEs]
    empty_intervals = []
    prev_end = existing_intervals[0][1]
    for start, end in existing_intervals:
        if start > prev_end:
            empty_intervals.append((start - prev_end, prev_end, start))
        prev_end = max(prev_end, end)
    
    # sample DEL and INS positions
    # not allow for too much background SVs
    max_DEL = len(empty_intervals) // 2
    if nDEL > max_DEL:
        raise ValueError(f"Too much background deletions. We recommend to reduce the number of background deletions < {max_DEL}.")
    # select longer empty intervals for DEL
    empty_intervals.sort(reverse=True)
    candidate_DEL_intervals = empty_intervals[:nDEL]
    candidate_INS_intervals = empty_intervals[nDEL:]
    # background SV positions
    # bgSV_positions = []
    # sample DEL positions
    SVmin = min(empty_intervals[nDEL][0], SVmin)
    for idx, i in enumerate(candidate_DEL_intervals):
        length, e_start, e_end = i
        del_len = random.randint(SVmin, min(SVmax, length))
        del_start = random.randint(e_start, e_end - del_len)
        TEs.append((chrom, del_start, del_start + del_len, f"bgDEL_{idx}_{del_len}"))
    # sample INS positions
    sampled_points = set()
    length_eachINS_intervals = [i[0] for i in candidate_INS_intervals]
    while len(sampled_points) < nINS:
        interval = random.choices(candidate_INS_intervals, weights=length_eachINS_intervals, k=1)[0]
        _, start, end = interval
        point = random.randint(start, end)
        sampled_points.add(point)
    for i, j in zip(sampled_points, bgINS_seqs):
        TEs.append((chrom, i, i+1, j[0]))
    
    # output the new bed file
    with open(bedout, "w") as fout:
        TEs.sort()
        for te in TEs:
            fout.write(f"{te[0]}\t{te[1]}\t{te[2]}\t{te[3]}\n")

# ---------------- select TEs with the restrict of minimum number ----------------
def make_min_TE(TE_list: list, nMIN: int, nTE: int, TEtype: set):
    # idct of the TEs
    te = {}
    for i in TE_list:
        tefamily = i[3]
        if tefamily not in TEtype:
            continue
        if tefamily not in te:
            te[tefamily] = [i]
        else:
            te[tefamily].append(i)
    # feasibility check
    if len(te) * nMIN > nTE:
        raise ValueError(f"Cannot satisfy the minimum number of TEs per family with the total number of TEs to be simulated. Please decrease nMIN.")
    # select minimum TEs
    selected_te = []
    for key in te.keys():
        col = te[key]
        if len(col) < nMIN:
            raise ValueError(f"The TE number of family {key} is less than nMIN ({nMIN}). Please decrease nMIN.")
        selected_te.extend(random.sample(col, nMIN))
    if len(selected_te) < nTE:
        remaining = nTE - len(selected_te)
        remaining_pool = [i for i in TE_list if i not in selected_te]
        selected_te.extend(random.sample(remaining_pool, remaining))
    return selected_te

