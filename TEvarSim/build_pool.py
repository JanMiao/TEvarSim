import numpy as np
import random
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import sample_TEins, bgSV, make_min_TE

def is_dna(seq: Seq) -> bool:
    return set(str(seq)) <= set("ATGCN")

def check_TEid(seq: str) -> bool:
    return "#" in seq and "/" in seq

def check_output_file(output_path):
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file '{output_path}' already exists. Please choose a different name.")

def CHRnorm(chrA,chrB):
    ### normalize chrA according the style if chrB
    if chrA.isdigit() and chrB.isdigit():
         return chrA
    elif chrA.isdigit() and not chrB.isdigit():
        if not chrB[3:].isdigit():
            raise ValueError("Chromosome names must be either plain numbers ('1') or prefixed with 'chr' ('chr1')")
        else:
            chrA = chrB[:3] + chrA
    elif not chrA.isdigit() and chrB.isdigit():
            chrA = chrA[3:]
    else:
        if chrA[:3].upper() != chrB[:3].upper():
            raise ValueError("Chromosome names must be either plain numbers ('1') or prefixed with 'chr' ('chr1')")
        else:
            chrA = chrB[:3] + chrA[3:]
    return chrA

class RandomTE:
    def __init__(self, args):
        self.pool_fasta = args.outprefix + ".fa"
        self.out_fasta = args.outprefix + ".bgSV.fa"
        self.DELfile = args.knownDEL
        self.CHR = args.CHR
        self.nTE = args.nTE
        self.ins_ratio = args.ins_ratio
        self.TEtype = args.TEtype
        self.prefix = args.outprefix
        self.DELlen = args.DELlen
        self.TEdistance = args.TEdistance
        self.nSV = args.nSV # number of background SVs
        self.nMIN = args.nMIN # minimum number of TEs for each family
        self.random_seed = args.seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
            random.seed(self.random_seed)

    def _run(self):
        # 1. parse DEL file
        _, ext = os.path.splitext(self.DELfile)
        if ext.lower() == ".txt":
            self.parse_DEL_ucsc()
        elif ext.lower() == ".out":
            self.parse_DEL_repeatmasker()
        else:
            raise ValueError("DEL file name must end with .txt (UCSC) or .out (RepeatMasker).")
        # 2. parse TEpool
        nINS = int(self.nTE * self.ins_ratio)
        self.parse_TEpool(nINS)
        # 3. Remove overlapping regions between INS and DEL
        self.remove_dup_DEL()
        # 4. Generate BED file
        self.build_bed()
        logging.info(f"Generated TE BED file: {self.prefix}.bed")
        # 5. Add background SVs if specified
        if self.nSV > 0:
            bedin = self.prefix + ".bed"
            bedout = self.prefix + ".bgSV.bed"
            bgSV(bedin, bedout, self.nSV, self.ins_ratio, self.pool_fasta, self.out_fasta)
            logging.info(f"Generated BED file and fasta file with background SVs are: {bedout} and {self.out_fasta}")

    def build_bed(self):
        # merge 
        nINS = int(self.nTE * self.ins_ratio)
        nDEL = self.nTE - nINS
        logging.info(f"Generating {nINS} INS and {nDEL} DEL for chromosome {self.CHR}")
        # print(len(self.INS), len(self.DEL))
        if nINS > len(self.INS) or nDEL > len(self.DEL): 
            raise ValueError(f"generated TE count exceeds total; require INS < {len(self.INS)} and DEL < {len(self.DEL)}")
        # check nMIN
        if self.nMIN >= self.nTE:
            raise ValueError(f"minumum number of a TE family ({self.nMIN}) should be less than nTE ({self.nTE})")
        if self.nMIN > 0:
            DELlist = make_min_TE(self.DEL, self.nMIN, nDEL, self.TEtype)
        else:
            DEL_indices = np.random.choice(np.arange(len(self.DEL)), size=nDEL, replace=False)
            DELlist = [self.DEL[i] for i in DEL_indices]
        # merge INS and DEL
        merged = self.INS + DELlist
        merged.sort()
        # bedfile output
        bed_name = self.prefix + ".bed"
        check_output_file(bed_name)
        with open(bed_name, "w") as f:
            for start, end, teID, *_ in merged:
                f.write("\t".join([self.CHR, str(start), str(end), teID]) + "\n")
        logging.info(f"Generated BED file: {bed_name}")

    def parse_DEL_ucsc(self):
        """
        Parse UCSC repeat annotation file (.txt) and return a list of parsed records.
        self.DEL
        """
        self.DEL = []
        if not self.TEtype:
            self.TEtype = {'Alu', 'L1', 'SVA'}
        ## nomalize the chr ID
        #with open(self.DELfile, "r") as f:
        #    first_line = f.readline().strip()
        #    chrom = first_line.split('\t')[5]
        #    print(f"chromosome from DEL file: {chrom}")
        #    self.CHR = CHRnorm(self.CHR, chrom)
        # process file
        with open(self.DELfile) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                chrom = fields[5]
                # nomalize the chr ID
                self.CHR = CHRnorm(self.CHR, chrom)
                # repClass = fields[12][:3]
                repClass = fields[12]
                if chrom == self.CHR and repClass in self.TEtype:
                    start = int(fields[6])
                    end = int(fields[7])
                    if end - start > self.DELlen:
                        teID = f"DEL-{chrom}-{start}-{end}-{fields[12]}-{fields[10]}"
                        TEfamily = fields[12]
                        self.DEL.append((start, end, teID, TEfamily, "DEL"))
                    # self.TEpool[teID] = fasta.fetch(self.CHR, start, end)

    def parse_DEL_repeatmasker(self):
        """
        Parse repeatmasker output file (.out) and return a list of parsed records.
        self.DEL
        """
        self.DEL = []
        if not self.TEtype:
            self.TEtype = {'Alu', 'L1', 'SVA'}
        # process file
        with open(self.DELfile) as f:
            for line in f:
                if not line.strip() or not line.strip()[0].isdigit(): 
                    continue
                fields = line.strip().split()
                chrom = fields[4]
                # nomalize the chr ID
                self.CHR = CHRnorm(self.CHR, chrom)
                if "/" not in fields[10]:
                    continue
                repClass = fields[10].split("/")[1]
                if chrom == self.CHR and repClass in self.TEtype:
                    start = int(fields[5])
                    end = int(fields[6])
                    if end - start > self.DELlen:
                        teID = f"DEL-{chrom}-{start}-{end}-{fields[10]}-{fields[9]}"
                        self.DEL.append((start, end, teID, repClass,"DEL"))
    
    def parse_TEpool(self, nINS):
        self.INS = []
        records = list(SeqIO.parse(self.pool_fasta, "fasta"))
        # The max insert position is the end of the last DEL
        # INSpos = np.random.choice(range(self.DEL[-1][1]), size=nINS, replace=False)
        INSpos = sample_TEins(start=1, end=self.DEL[-1][1], n=nINS, TEdistance=self.TEdistance)
        for record,pos in zip(records,INSpos):
            self.INS.append((pos, pos + 1, record.id, record.id.split("/")[1], "INS"))

    def remove_dup_DEL(self):
        """
        Remove DELs that overlap with INSs.
        """
        i, j = 0, 0
        rm_idx = []

        while i < len(self.INS) and j < len(self.DEL):
            a_start, a_end, *_ = self.INS[i]
            b_start, b_end, *_ = self.DEL[j]

            # If B[j] overlaps with A[i]
            if a_start < b_end and a_end > b_start:
                rm_idx.append(j)

            if a_end < b_end:
                i += 1
            else:
                j += 1
        # remove
        delete_set = set(rm_idx)
        self.DEL = [item for i, item in enumerate(self.DEL) if i not in delete_set]


class TEPoolBuilder:
    def __init__(self, args):
        # base
        self.consensus = args.consensus
        self.nTE = args.nTE
        self.ins_ratio = args.ins_ratio
        self.nINS = round(self.nTE * self.ins_ratio)
        self.prefix = args.outprefix
        # SNP and INDEL
        self.snp_rate = args.snp_rate
        self.indel_rate = args.indel_rate
        self.indel_ins = args.indel_ins
        self.indel_geom_p = args.indel_geom_p
        # truncate
        self.truncated_ratio = args.truncated_ratio
        self.truncated_max_length = args.truncated_max_length
        # polyA
        self.polyA_ratio = args.polyA_ratio
        self.polyA_min = args.polyA_min
        self.polyA_max = args.polyA_max
        # other
        self.random_seed = args.seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
        logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
        
        # data stored
        self.BASES = ['A', 'C', 'G', 'T']
        self.CHOICES_DICT = {b: [c for c in self.BASES if c != b] for b in self.BASES}
        self.current_seq = []
        self.seqID_suffix = ""

    def _run(self):
        records = list(SeqIO.parse(self.consensus, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in consensus FASTA: {self.consensus}")
        logging.info(f"Found {len(records)} sequences in consensus FASTA.")
        # check if all sequence are legal
        for record in records:
            record.seq = record.seq.upper()
            if not is_dna(record.seq):
                raise ValueError(f"Invalid characters found in {record.id}")
            if not check_TEid(record.id):
                raise ValueError(f"Invalid format in TE ID: {record.id}. Please follow the format '>TEname#class/superfamily', e.g., '>AluY#SINE/Alu'")

        # generate random masks for truncation and polyA
        truncated_num = np.random.random(self.nINS)
        polyA_num = np.random.random(self.nINS)
        trunc_mask = truncated_num < self.truncated_ratio
        polyA_mask = polyA_num < self.polyA_ratio
        mutate = np.select(
            [trunc_mask & polyA_mask, trunc_mask, polyA_mask],
            [3, 2, 1],
            default=0
        )

        # map numbers to corresponding functions
        func_map = {
            0: lambda: self.INDEL_mutate().SNP_mutate(),
            1: lambda: self.INDEL_mutate().SNP_mutate().apply_polyA(),
            2: lambda: self.INDEL_mutate().SNP_mutate().apply_truncate(),
            3: lambda: self.INDEL_mutate().SNP_mutate().apply_truncate().apply_polyA()
        }

        out_records = []
        for j in mutate:
            idx = np.random.randint(0, len(records))
            record = records[idx]
            self.current_seq = list(str(record.seq))
            self.seqID_suffix = ""  # reset suffix each loop
            func_map[j]()
            new_id = f"{record.id}_{self.seqID_suffix}"
            new_record = SeqRecord(Seq("".join(self.current_seq)), id=new_id, description="")
            out_records.append(new_record)

        # output
        output_fasta = self.prefix + ".fa"
        check_output_file(output_fasta)
        SeqIO.write(out_records, output_fasta, "fasta")
        logging.info(f"Generated {len(out_records)} sequences -> {output_fasta}")
    

    def apply_truncate(self):
        max_trunc = int(len(self.current_seq) * self.truncated_max_length)
        if max_trunc < 1:
            return self
        trunc_len = np.random.randint(1, max_trunc + 1)  # inclusive upper bound
        self.current_seq = self.current_seq[trunc_len:]
        self.seqID_suffix += f"{trunc_len}truncate"
        return self

    def apply_polyA(self):
        polyA_len = np.random.randint(self.polyA_min, self.polyA_max + 1)
        self.current_seq.extend(["A"] * polyA_len)
        self.seqID_suffix += f"{polyA_len}polyA"
        return self
    
    def SNP_mutate(self):
        L = len(self.current_seq)
        n_snp = np.random.poisson(self.snp_rate * L)
        if n_snp < 1:
            return self
        self.seqID_suffix += f"{n_snp}SNP"
        # SNP positions, ensure not exceeding sequence length
        snp_positions = np.random.choice(L, size=min(n_snp, L), replace=False)
        for pos in snp_positions:
            current_base = self.current_seq[pos]
            if current_base in self.CHOICES_DICT:
                self.current_seq[pos] = np.random.choice(self.CHOICES_DICT[current_base])
        return self
    
    def INDEL_mutate(self):
        # total INDEL number
        L = len(self.current_seq)
        n_indel = np.random.poisson(self.indel_rate * L)
        if n_indel < 1:
            return self
        self.seqID_suffix += f"{n_indel}INDEL"
        if n_indel == 0:
            return self
        # generate INDEL lengths
        indel_len_list = np.random.geometric(self.indel_geom_p, n_indel)
        indel_len_list = np.clip(indel_len_list, a_min=None, a_max=30)
        # select positions for INDELs
        indel_positions = np.random.choice(L, size=min(n_indel, L), replace=False)
        for pos, indel_len in zip(sorted(indel_positions, reverse=True), indel_len_list):
            if np.random.random() < self.indel_ins:
                ins_seq = [np.random.choice(self.BASES) for _ in range(indel_len)]
                self.current_seq[pos:pos] = ins_seq
            else:
                del_end = min(pos + indel_len, len(self.current_seq))
                del self.current_seq[pos:del_end]
        return self

def run(args):
    # generate TE pool
    TEPoolBuilder(args)._run()
    # generate bed file with random TE insertions
    RandomTE(args)._run()


