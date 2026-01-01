import pysam

def load_truth_vcf(vcf_file, sampleID, INSonly, TEtype):
    variants = {}
    vcf = pysam.VariantFile(vcf_file)

    if sampleID not in vcf.header.samples:
        raise ValueError(f"Sample ID '{sampleID}' not found in VCF header")
    
    for record in vcf:
        varID = record.id
        if INSonly:
            if varID.startswith("DEL"):
                continue
            TEfamily = varID.split("/")[1].split("_")[0]
            if TEfamily == TEtype:
                sample = record.samples[sampleID]
                gt = sample.get('GT')
                # filter positions withou any variants
                if any(allele is None for allele in gt) or all(allele == 0 for allele in gt):
                    continue
                if record.chrom in variants:
                    variants[record.chrom].append((record.pos, gt))
                else:
                    variants[record.chrom] = [(record.pos, gt)]
        else:
            if varID.startswith("DEL"):
                TEfamily = varID.split("/")[1].split("-")[0]
            else:
                TEfamily = varID.split("/")[1].split("_")[0]
            if TEfamily == TEtype:
                sample = record.samples[sampleID]
                gt = sample.get('GT')  
                # filt positions withou any variants
                if any(allele is None for allele in gt) or all(allele == 0 for allele in gt):
                    continue
                if record.chrom in variants:
                    variants[record.chrom].append((record.pos, gt))
                else:
                    variants[record.chrom] = [(record.pos, gt)]
    return variants

def load_vcf(vcf_file, sampleID):
    variants = {}
    vcf = pysam.VariantFile(vcf_file)

    if sampleID not in vcf.header.samples:
        raise ValueError(f"Sample ID '{sampleID}' not found in VCF header")
    
    for record in vcf.fetch():    
        sample = record.samples[sampleID]
        gt = sample.get('GT')  
        # filt positions withou any variants
        if any(allele is None for allele in gt) or all(allele == 0 for allele in gt):
            continue
        if record.chrom in variants:
            variants[record.chrom].append((record.pos, gt))
        else:
            variants[record.chrom] = [(record.pos, gt)]
    return variants

def load_bed(bed_file):
    variants = {}
    with open(bed_file, 'r') as fin:
        for line in fin:
            chrom, start, _, gt, *_ = line.strip().split()
            # tuple of genotypes
            if "/" in gt and "." not in gt:
                gt = tuple(map(int, gt.split("/")))
            else:
                gt = (None,)
            start = int(start)
            if chrom in variants:
                variants[chrom].append((start, gt))
            else:
                variants[chrom] = [(start, gt)]
    return variants

def calculate_metrics(tp, fp, fn):
    '''
    Calculate F1, precision, recall from confusion counts.
    '''
    # Parameter check
    assert tp >= 0 and fp >= 0 and fn >= 0, "tp, fp, fn must be non-negative integers"
    # Calculate Recall and Precision
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    # Calculate F1
    if (recall + precision) > 0:
        F1 = 2 * recall * precision / (recall + precision)
    else:
        F1 = 0
    return F1, precision, recall


class CompareVCF:
    def __init__(self, args):
        self.nHap = args.nHap
        self.truth_file = args.truth
        self.compare_file = args.pred
        self.predType = args.predType
        self.TEtype = args.TEtype
        self.INSonly = args.INSonly
        self.max_dist = args.max_dist
        self.truthID = args.truthID
        self.predID = args.predID
        self.MatchFile = args.outprefix
        # inner variables
        self.nMatch = 0
    
    def _run(self):
        # ground truth data
        if self.nHap > 1:
            self.convert_to_ploidy()
            bench_var = load_truth_vcf("polyhap.vcf", self.truthID, self.INSonly, self.TEtype)
        else:
            bench_var = load_truth_vcf(self.truth_file, self.truthID, self.INSonly, self.TEtype)
        # prediction data
        if self.predType == "VCF":
            pred_var = load_vcf(self.compare_file, self.predID)
        else:
            self.predID = "sample" if self.predID is None else self.predID
            pred_var = load_bed(self.compare_file)
        # metrics
        tp, fp, fn, gtDiff = 0, 0, 0, 0
        # different chromosomes
        bench_chr = set(bench_var.keys())
        pred_chr = set(pred_var.keys())
        fn_chr = bench_chr - pred_chr
        if len(fn_chr) > 0:
            fn += sum(len(bench_var[chrom]) for chrom in fn_chr)
        fp_chr = pred_chr - bench_chr
        if len(fp_chr) > 0:
            fp += sum(len(pred_var[chrom]) for chrom in fp_chr)
        # common chromosomes
        common_chr = pred_chr & bench_chr
        print(f"Common chromosomes: {common_chr}")
        print(f"Nmber of {self.TEtype} in truth VCF: {sum(len(bench_var[chrom]) for chrom in bench_chr)}")
        print(f"Number of {self.TEtype} in prediction VCF: {sum(len(pred_var[chrom]) for chrom in pred_chr)}")
        for chrom in common_chr:
            bench_pos, bench_gt = zip(*bench_var[chrom])
            compare_pos, compare_gt = zip(*pred_var[chrom])
            tp_tmp, fp_tmp, fn_tmp, gtDiff_tmp = self.calculate_confusion_counts(bench_pos, bench_gt, compare_pos, compare_gt, chrom)
            tp += tp_tmp
            fp += fp_tmp
            fn += fn_tmp
            gtDiff += gtDiff_tmp

        self.F1, self.precision, self.recall = calculate_metrics(tp, fp, fn)
        print(f"For {self.TEtype} occurrence sites:")
        print(f"True Positive: {tp}, False Positive: {fp}, False negative: {fn}")
        print(f"Recall: {self.recall:.4f}, Precision: {self.precision:.4f}, F1: {self.F1:.4f}")
        print(f"The genotype accuracy for matched {self.TEtype}: {1 - gtDiff/self.nMatch}")
    
    # 应该把这个函数放到VCF文件生成里
    def convert_to_ploidy(self):
        # self.nHap
        # self.truth_file
        with open(self.truth_file, 'r') as fin, open("polyhap.vcf", 'w') as fout:
            for line in fin:
                if line.startswith('##'):
                    # 直接输出注释行
                    fout.write(line)
                elif line.startswith('#CHROM'):
                    # 处理表头
                    fields = line.strip().split('\t')
                    fixed_cols = fields[:9]
                    samples = fields[9:]
                    # 检查样本数能否被ploidy整除
                    if len(samples) % self.nHap != 0:
                        raise ValueError(f"Error: number of haplotypes ({len(samples)}) is not divisible by ploidy ({self.nHap})")
                    # 合并样本名
                    merged_samples = []
                    for i in range(0, len(samples), self.nHap):
                        group = samples[i:i+self.nHap]
                        merged_name = "_".join(group)
                        merged_samples.append(merged_name)
                    fout.write('\t'.join(fixed_cols + merged_samples) + '\n')
                else:
                    # 处理数据行
                    fields = line.strip().split('\t')
                    fixed_cols = fields[:9]
                    genotypes = fields[9:]
                    merged_gts = []
                    for i in range(0, len(genotypes), self.nHap):
                        group = genotypes[i:i+self.nHap]
                        # 只保留GT字段里的数字，假设只有GT且格式简单
                        # 如果是复杂格式，需要再扩展
                        merged_gt = '|'.join(group)
                        merged_gts.append(merged_gt)
                    fout.write('\t'.join(fixed_cols + merged_gts) + '\n')

    def calculate_confusion_counts(self, bench_pos, bench_gt, compare_pos, compare_gt, chrom):
        # print("Calculating confusion counts...")
        # get matched data
        i, j = 0, 0
        b_len = len(bench_pos)
        c_len = len(compare_pos)
        match_idxb = []
        match_idxc = []
        fo = open(self.MatchFile + ".csv", "a")
        while i < b_len and j < c_len:
            pos_a = bench_pos[i]
            pos_b = compare_pos[j]
            if pos_b > pos_a + self.max_dist:
                i += 1
            elif pos_a - self.max_dist <= pos_b <= pos_a + self.max_dist:
                outLine = [str(i) for i in [bench_pos[i], bench_gt[i], compare_pos[j], compare_gt[j]]]
                fo.write(chrom + "," + ",".join(outLine) + "\n")
                match_idxb.append(i)
                match_idxc.append(j)
                i += 1
                j += 1
            else:
                j += 1
        fo.close()
        
        # calculate confusion counts
        tp = len(match_idxb)
        self.nMatch += tp
        #print(f"Number of matched TEs: {self.nMatch}")
        # FP: only in comparison VCF
        fp = c_len - tp
        # FN: only in bench vcf
        fn = b_len - tp
        bgt = [bench_gt[i] for i in match_idxb]
        cgt = [compare_gt[i] for i in match_idxc]
        # genotype accuracy
        if len(bgt[0]) != len(cgt[0]):
            # print("Warning: the number of allele is different between truth and prediction files")
            # print("Warning: filling the missing alleles in truth VCF, e.g., 1 -> 1/1")
            for i in range(len(bgt)):
                bgt[i] = bgt[i] + bgt[i]
        genotype_same = sum(1 for a, b in zip(bgt, cgt) if sorted(a) == sorted(b))
        gt_diff = tp - genotype_same
        
        return tp, fp, fn, gt_diff

def run(args):
    CompareVCF(args)._run()

