
import subprocess
import os
import random
from Bio import SeqIO

def coverTE(rptmak_record):
    if rptmak_record[8] == "C":
        non = int(rptmak_record[11][1:-1])
        START = int(rptmak_record[13])
    else:
        non = int(rptmak_record[13][1:-1])
        START = int(rptmak_record[11])
    
    END = int(rptmak_record[12])
    cover = END - START + 1
    full_len = END + non
    return float(cover/full_len)


class panTE:
    def __init__(self, args):
        # args
        self.GFAfile = args.gfa
        self.TEtype = args.TEtype
        self.minLen = args.minLen
        self.lib = args.lib
        self.cov = args.cov
        self.CHR = args.CHR
        self.nTE = args.nTE
        self.ins_ratio = args.ins_ratio
        self.tmpDir = args.tmpDir
        self.outprefix = args.outprefix
        self.random_seed = args.seed
        if self.random_seed is not None:
            random.seed(self.random_seed)

    def _run(self):
        # set default TE types
        if not self.TEtype:
            self.TEtype = {'LINE', 'SINE', 'LTR', 'RC'}
        os.makedirs(self.tmpDir, exist_ok=True)
        # 1. get SV from gfa
        self.extractSV()
        # 2. get INDEL sequences
        self.getINDEL()
        # 3. repeatmasker annotations
        self.anno()
        # 4. generate FAST and BED file
        self.generate()
        
    
    def extractSV(self):
        strings = f"gfatools bubble {self.GFAfile} > {self.tmpDir}/{self.outprefix}_SV.bed"
        subprocess.run(strings, shell=True, check=True)

    def getINDEL(self):
        with open(f"{self.tmpDir}/{self.outprefix}_SV.bed") as f, open(f"{self.tmpDir}/{self.outprefix}_INDEL.fa", "w") as fo:
            for line in f:
                info = line.strip().split("\t")
                CHR = info[0]
                if CHR != self.CHR:
                    continue
                if info[3] != "3":
                    continue
                if int(info[7]) < self.minLen:
                    continue
                START = info[1]
                END = info[2]
                SEQ= info[13]
                if int(END) > int(START):
                    TYPE= "DEL"
                else:
                    TYPE= "INS"
                fo.write(f">{CHR}:{START}-{END}-{TYPE}\n")
                fo.write(SEQ + "\n")
    
    def anno(self):
        strings = f"RepeatMasker -nolow -lib {self.lib} -s -dir {self.tmpDir} {self.tmpDir}/{self.outprefix}_INDEL.fa"
        subprocess.run(strings, shell=True, check=True)
        
    def generate(self):
        RPT = f"{self.tmpDir}/{self.outprefix}_INDEL.fa.out"
        indelFA = f"{self.tmpDir}/{self.outprefix}_INDEL.fa"
        BED = f"{self.outprefix}.bed"
        FASTA = f"{self.outprefix}.fa"

        seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(indelFA, "fasta")}
        outInfo = []
        tmpID = ""
        tmpSeqLen = 0
        with open(RPT) as f:
            for line in f:
                # skip header and empty line
                strp = line.strip()
                if not strp or not strp[0].isdigit():
                    continue
                info = strp.split()
                # TEtype filter
                tetype = info[10].split("/")[0]
                if tetype in self.TEtype:
                    # TE coverage filter
                    if coverTE(info) > self.cov:
                        seqID = info[4]
                        INDELtype = seqID.split("-")[-1]
                        left = int(info[5])
                        right = int(info[6])
                        SeqLen = right - left + 1
                        if INDELtype == "INS": # INS, select the longest TE for insertion
                            if seqID == tmpID:
                                if SeqLen > tmpSeqLen:
                                    outInfo.pop()
                                else:
                                    continue
                            new_seq = seq_dict[seqID][left-1:right]
                            pos, _ = seqID.rsplit("-", 1)
                            seqID = f"{pos}-{info[9]}"
                        else: # del renew seqID
                            new_seq = seq_dict[seqID][left-1:right]
                            CHR, other = seqID.split(":")
                            start, end, _ = other.split("-")
                            start = int(start)
                            end = int(end)
                            seqID = f"{CHR}:{left + start - 1}-{left + end}-{info[9]}"             

                        tmpID = seqID
                        tmpSeqLen = SeqLen
                        outInfo.append((seqID, info[9], new_seq, INDELtype))
        
        # output
        fa = open(FASTA, "w")
        bed = open(BED, "w")
        if self.nTE:
            random.shuffle(outInfo)
            nINS = int(self.nTE * self.ins_ratio)
            nDEL = self.nTE - nINS
            cDEL = 0
            cINS = 0
            for i in outInfo:
                SVtype = i[3]
                if SVtype == "INS":
                    cINS += 1
                else:
                    cDEL += 1
                if cDEL > nDEL and cINS >nINS:
                    break
                else:
                    if cDEL > nDEL and SVtype=="DEL":
                        continue
                    if cINS > nINS and SVtype=="INS":
                        continue
                fa.write(f">{i[0]}\n")
                fa.write(i[2] + "\n")
                chrom, pos = i[0].split(":")
                start, end, _ = pos.split("-", 2)
                bed.write(f"{chrom}\t{int(start) - 1}\t{end}\t{i[0]}\n")
        else:
            for i in outInfo:
                fa.write(f">{i[0]}\n")
                fa.write(i[2] + "\n")
                chrom, pos = i[0].split(":")
                start, end, _ = pos.split("-", 2)
                bed.write(f"{chrom}\t{int(start) - 1}\t{end}\t{i[0]}\n")


def run(args):
    panTE(args)._run()