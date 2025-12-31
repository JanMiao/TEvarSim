from Bio import SeqIO
import subprocess
import sys

class Read:
    def __init__(self, args):
        # self.reference = args.reference
        self.type = args.type
        self.genomeF = args.genome
        self.depth = args.depth
        # self.ngenome = args.ngenome
        # self.outprefix = args.outprefix
        # long reads
        self.Lerror = args.Lerror
        self.Lmean = args.Lmean
        self.Lstd = args.Lstd
        # short reads
        self.length = args.length
        self.Fmean = args.Fmean
        self.Fstd = args.Fstd
        # seed
        self.seed = args.seed

    def _run(self):
        # get genomes
        records = list(SeqIO.parse(self.genomeF, "fasta"))
        # single genome
        if len(records) == 1:
            record = records[0]
            if self.type == "long":
                self.simLong(self.genomeF, record.id)
            else:
                print(f"Simulating short reads for genome: {record.id}")
                self.simShort(self.genomeF, record.id, record.seq)
        else:
            # multiple genomes
            for record in records:
                singleSeq = record.id + ".fa"
                SeqIO.write(record, singleSeq, "fasta")
                if self.type == "long":
                    self.simLong(singleSeq, record.id)
                else:
                    self.simShort(singleSeq, record.id, record.seq)
    
    def simShort(self, genomeFile, seqID, seq):
        genome_length = len(seq)
        depth = int(genome_length * self.depth / (self.length * 2))
        print(f"Simulating short reads for genome: {seqID}")
        if self.seed:
            strings = (
                f"mason_simulator -ir {genomeFile} -o {seqID}_R1.fq.gz -or {seqID}_R2.fq.gz "
                f"-n {depth} --num-threads 10 --illumina-read-length {self.length} --seed {self.seed} "
                f"--fragment-mean-size {self.Fmean} --fragment-size-std-dev {self.Fstd}"
            )
        else:
            strings = (
                f"mason_simulator -ir {genomeFile} -o {seqID}_R1.fq.gz -or {seqID}_R2.fq.gz "
                f"-n {depth} --num-threads 10 --illumina-read-length {self.length} "
                f"--fragment-mean-size {self.Fmean} --fragment-size-std-dev {self.Fstd}"
            )
        subprocess.run(strings, shell=True, check=True)

    def simLong(self, genomeFile, seqID):
        print(f"Simulating long reads for genome: {seqID}")
        if self.seed:
            strings = (
                f"pbsim --strategy wgs --method qshmm --qshmm {sys.prefix}/data/QSHMM-RSII.model --depth {self.depth} "
                f"--accuracy-mean {1-self.Lerror} --length-mean {self.Lmean} --length-sd {self.Lstd} "
                f"--genome {genomeFile} --prefix {seqID} --seed {self.seed}"
                )
        else:
            strings = (
                f"pbsim --strategy wgs --method qshmm --qshmm {sys.prefix}/data/QSHMM-RSII.model --depth {self.depth} "
                f"--accuracy-mean {1-self.Lerror} --length-mean {self.Lmean} --length-sd {self.Lstd} "
                f"--genome {genomeFile} --prefix {seqID}"
                )
        subprocess.run(strings, shell=True, check=True)

def run(args):
    Read(args)._run()