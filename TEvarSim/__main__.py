# __main__.py
import argparse
import sys
from . import build_pool, TE_real, simulate, compare_vcf, reads, TEpan
from . import __version__

def main():
    parser = argparse.ArgumentParser(prog="tevarsim", 
                                     description=f"TEvarSim: A genome simulation tool for transposable element (TE) variants\nVersion: {__version__}")
    parser.add_argument("--version", "-V", action="version", version=f"%(prog)s {__version__}")
    
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. TErandom
    p1 = subparsers.add_parser("TErandom", 
                               help="Generate pTE position from known deletion sites and random TE insertion")
    # Base
    p1.add_argument("--consensus", "-C", type=str,  required=True,
                    help="Path to the TE consensus FASTA file. The sequenceIDs in the FASTA header should be >TEname#class/superfamily, e.g., >AluY#SINE/Alu")
    p1.add_argument("--knownDEL", "-L", type=str, required=True, 
                    help="Input known TE deletion file (RepeatMasker .out or UCSC .txt)")
    p1.add_argument("--CHR", "-H", type=str, required=True,
                    help="Chromosome to simulate TE insertions on (e.g., chr21 or 21)")
    p1.add_argument("--TEtype", "-e", type=str, action="append",
                    help="Which TE super families to be extracted from the TE deletion file, with the default set as Alu, L1, ERV, and SVA")
    p1.add_argument("--nTE", "-N", type=int, default=100, 
                    help="Number of polymorphic TE (pTE) insertions to simulate (default: 100)")
    p1.add_argument("--ins-ratio", "-R", type=float, default=0.6, 
                    help="Proportion of insertion events among all simulated pTE (0-1, default: 0.6)")
    p1.add_argument("--outprefix", "-O", type=str, default="TErandom", 
                    help="Output prefix for the generated TE pool FASTA file and the bed file (default: TErandom)")
    p1.add_argument("--DELlen", type=int, default=100,
                    help="A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)")
    p1.add_argument("--nMIN", type=int, default=0,
                    help="A minimum number of TE deletions for each TE super family to be simulated (default: 0)")
    p1.add_argument("--TEdistance", type=int, default=500,
                    help="A minimum length of distance between two TE insertions (default: 500 bp)")
    p1.add_argument("--nSV", type=int, default=0, help="Number of background structural variants to simulate (default: 0)")
    # SNP and INDEL
    p1.add_argument("--snp-rate", "-S", type=float, default=0.02, 
                    help="SNP mutation rate per base (default: 0.02)")
    p1.add_argument("--indel-rate", "-I", type=float, default=0.005, 
                    help="Indel mutation rate per base (default: 0.005)")
    p1.add_argument("--indel-ins", "-r", type=float, default=0.4, 
                    help="Proportion of insertion events among INDELs (0-1, default: 0.4)")
    p1.add_argument("--indel-geom-p", "-G", type=float, default=0.7, 
                    help="Parameter 'p' of geometric distribution for indel lengths (default: 0.7)")
    # Truncation
    p1.add_argument("--truncated-ratio", "-T", type=float, default=0.3, 
                    help="Proportion of TE sequences to truncate (0-1, default: 0.3)")
    p1.add_argument("--truncated-max-length", "-K", type=float, default=0.5, 
                    help="Maximum proportion of sequence length to truncate (0-1, default: 0.5)")
    # PolyA
    p1.add_argument("--polyA-ratio", "-A", type=float, default=0.8, 
                    help="Proportion of TE sequences to add polyA tail (0-1, default: 0.8)")
    p1.add_argument("--polyA-min", "-M", type=int, default=5, 
                    help="Minimum polyA tail length (default: 5)")
    p1.add_argument("--polyA-max", "-X", type=int, default=20, 
                    help="Maximum polyA tail length (default: 20)")
    # Other
    p1.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    p1.set_defaults(func=build_pool.run)

    # 2. TEreal
    p2 = subparsers.add_parser("TEreal", 
                               help="Generate pTE position from Known TE insertion and deletion")
    # Input
    p2.add_argument("--knownINS", "-K", type=str, required=True, 
                    help="Input known TE insertion file (e.g., MEI_Callset)")
    p2.add_argument("--knownDEL", "-L", type=str, required=True, 
                    help="Input known TE deletion file (RepeatMasker .out or UCSC .txt)")
    p2.add_argument("--TEtype", "-e", type=str, action="append",
                    help="TEs to be extracted from the TE deletion file, with the default set as LINE, SINE, LTR, and Helitron.")
    p2.add_argument("--CHR", "-C", type=str, required=True,
                    help="Chromosome to simulate TE insertions on (e.g., chr21 or 21)")
    p2.add_argument("--DELlen", type=int, default=100,
                    help="A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)")
    p2.add_argument("--nMIN", type=int, default=0,
                    help="A minimum number of TE deletions for each TE super family to be simulated (default: 0)")
    p2.add_argument("--nSV", type=int, default=0, help="Number of background structural variants to simulate (default: 0)")
    # Output
    p2.add_argument("--outprefix", "-O", type=str, default="TEreal", 
                    help="Output prefix for generated BED file (default: 'TEreal')")
    p2.add_argument("--nTE", "-N", type=int,
                    help="Number of polymorphic TE (pTE) insertions to simulate (default: all TEs)")
    p2.add_argument("--ins-ratio", "-R", type=float, default=0.4, 
                    help="Proportion of insertion events among all simulated pTE (0-1, default: 0.4)")
    # Other
    p2.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    p2.set_defaults(func=TE_real.run)

    # 3. TE from pangenome graph
    p3 = subparsers.add_parser("TEpan", 
                               help="Generate pTE position from Pangenome graph")
    # Input
    p3.add_argument("--gfa", "-G", type=str, required=True, 
                    help="GFA file of the pangenome graph")
    p3.add_argument("--lib", "-L", type=str, required=True, 
                    help="RepeatMasker library file")
    p3.add_argument("--CHR", "-H", type=str, required=True, 
                    help="Chromosome to simulate TE insertions on (e.g., chr21 or 21)")
    p3.add_argument("--minLen", "-I", type=int, default= 250,
                    help="Minimum length of structural variants to consider (default: 250)")
    p3.add_argument("--TEtype", "-e", type=str, action="append",
                    help="TEs to be extracted from the RepeatMasker annotation file, with the default set as LINE, SINE, LTR, and Helitron.")
    p3.add_argument("--cov", "-C", type=float, default= 0.5,
                    help="Minimum TE coverage to consider a structural variant as TE (0-1, default: 0.5)")
    p3.add_argument("--tmpDir", "-T", type=str, default="tmp_TEpan", 
                    help="Temporary directory for intermediate files (default: tmp_TEpan)")
    p3.add_argument("--nTE", "-N", type=int, 
                    help="Number of polymorphic TE (pTE) insertions to simulate (default: all TEs)")
    p3.add_argument("--ins-ratio", "-R", type=float, default=0.4, 
                    help="Proportion of insertion events among all simulated pTE (0-1, default: 0.4)")
    # Output
    p3.add_argument("--outprefix", "-O", type=str, default="TEpan", 
                    help="Prefix for output files (VCF + modified genome FASTA)")
    p3.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    p3.set_defaults(func=TEpan.run)


    # 4. simulate
    p4 = subparsers.add_parser("Simulate", 
                               help="Simulate TE insertions/deletions and generate VCF and modified genome FASTA")
    # Input
    p4.add_argument("--ref", "-F", type=str, required=True, 
                    help="Reference genome FASTA file")
    p4.add_argument("--pool", "-P", type=str, required=True, 
                    help="FASTA file of TE sequences generated from 'TEpool'")
    p4.add_argument("--bed", "-B", type=str, required=True, 
                    help="BED file containing TE positions (can be generated by 'TEreal')")
    # Output
    p4.add_argument("--outprefix", "-O", type=str, default="Sim", 
                    help="Prefix for output files (VCF + modified genome FASTA)")
    # Options
    p4.add_argument("--num", "-N", type=int, required=True, 
                    help="Number of genomes to simulate")
    p4.add_argument("--diverse", "-I", action="store_true",
                    help="Introduce sequence diversity among individuals for the same TE event")
    p4.add_argument("--diverse_config", "-c", type=str,
                    help="Path to a configuration file for introducing sequence diversity among individuals for the same TE event")
    p4.add_argument("--af-min", "-A", type=float, default=0.1, 
                    help="Minimum allele frequency for simulated TE variants (default: 0.1)")
    p4.add_argument("--af-max", "-X", type=float, default=0.9, 
                    help="Maximum allele frequency for simulated TE variants (default: 0.9)")
    p4.add_argument("--tsd-min", "-M", type=int, default=5, 
                    help="Minimum TSD length (default: 5)")
    p4.add_argument("--tsd-max", "-Y", type=int, default=20, 
                    help="Maximum TSD length (default: 20)")
    p4.add_argument("--sense-strand-ratio", "-S", type=float, default=0.5, 
                    help="Proportion of TE variants in the sense strand (default: 0.5)")
    # Other
    p4.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    p4.set_defaults(func=simulate.run)

    # 5. compare
    p5 = subparsers.add_parser("Compare", help="Compare predicted VCF to ground truth VCF")
    # Input
    p5.add_argument("--truth", "-T", type=str, required=True, help="Ground truth VCF file")
    p5.add_argument("--pred", "-P", type=str, required=True, help="Predicted file to compare")
    p5.add_argument("--predType", "-p", type=str,  choices=["VCF", "BED"], default="VCF", 
                    help="Type of the predicted file (VCF or BED, default: VCF)")
    # Output
    p5.add_argument("--outprefix", "-O", type=str, required=True, help="Output matched TEs")
    # Options
    p5.add_argument("--truthID", "-I", type=str, required=True, help="Sample ID in the truth VCF")
    p5.add_argument("--predID", "-J", type=str, required=True, help="Sample ID in the predicted VCF")
    p5.add_argument("--TEtype", "-e", type=str, default=None, help="TE type in truth VCF to consider in the comparison")
    p5.add_argument("--INSonly", action="store_true",  help="Only compare insertions in truth VCF")
    p5.add_argument("--nHap", "-N", type=int, default=2,  help="Number of haplotypes in the genome (default: 2)")
    p5.add_argument("--max_dist", "-M", type=int, default=100, 
                    help="Maximum allowed distance (bp) to consider two variants as matching")
    p5.set_defaults(func=compare_vcf.run)

    # 6. Read simulation
    p6 = subparsers.add_parser("Readsim", 
                               help="generate short or long reads from the simulated genome")
    # general
    p6.add_argument("--type", "-T", choices=["short", "long"],
                    type=str, required=True, help="Simulate short reads or long reads")
    p6.add_argument("--genome", "-G", type=str, required=True, 
                    help="The file contains genomes where reads simulated from")
    p6.add_argument("--depth", "-P", type=int, required=True, 
                    help="Depth of simulated reads")
    #p6.add_argument("--ngenome", "-N", type=int, default= 1, 
    #                help="Number of genomes for simultaneous reads simulation")
    #p6.add_argument("--outprefix", "-O", type=str, required=True, 
    #                help="prefix of output files")    
    # long reads settings
    p6.add_argument("--Lerror", "-E", type=float, default= 0.15, 
                    help="sequencing error rate for long reads (default: 0.15)")
    p6.add_argument("--Lmean", "-M", type=int,  default= 9000,
                    help="average size of read length (only for long reads, default: 9000 bp)")
    p6.add_argument("--Lstd", "-S", type=int, default= 7000,  
                    help="read length standard deviation (only for long reads, default: 7000 bp)")
    # short reads settings
    p6.add_argument("--length", "-i", type=int, default= 150, 
                    help="read length (only for short reads, default: 150 bp)")
    p6.add_argument("--Fmean", "-m", type=int, default= 500, 
                    help="average size of fragment length (only for short reads, default: 500 bp)")
    p6.add_argument("--Fstd", "-s", type=int, default= 30,  
                    help="Fragment size standard deviation (only for short reads, default: 30 bp)")
    # random seed
    p6.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    p6.set_defaults(func=reads.run)



    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    if args.command == "simulate":
        if args.diverse_config and not args.diverse:
            parser.error("--diverse_config requires --diverse to be set")
    args.func(args)


if __name__ == "__main__":
    main()
