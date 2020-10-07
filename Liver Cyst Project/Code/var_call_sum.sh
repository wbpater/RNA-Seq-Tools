#!/bin/bash

##### Variables

genome=
region=
bamfolder=
outdir=$(pwd)/
annotation=
verbose=0

##### Functions
usage()
{
  echo '
Usage: var_call_sum.sh [options]... [FILE]

This script performs variant calling on a number of BAM files and makes a summarised output for all input files
Required files can be dowloaded from the following sources:

Genome FASTA: gencode whole genome
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.p12.genome.fa.gz

Region FASTA (for calling exonic and splice variants): navigate to a gene in the
UCSC genome browser and click on the major splicing variant,
click on "Genomic sequence" select 5" UTR Exons, CDS Exons, 3" UTR Exons.
Select One FASTA record per region, with 3 extra bases upstream and downstream.
Save the FASTA  as <genename>.fasta (will be used for naming files) or copy to a flat file.

ClinVar vcf: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/ look for the vcf_<genome_build> folder and the clinvar.vcf.gz file therin (build CRCh38 at time of writing)
(https://www.ncbi.nlm.nih.gov/variation/docs/ClinVar_vcf_files/)

Required arguments:
-b|--bamfolder folder with BAM files
			Bam files to call variants for sorted by coordinate (default: None)

-g|--genome FASTA file  The genome used for mapping that created the BAM files.
                        (default: None)

-a|--annotation VCF file
			The annotation file for calling known variants
			(default: None)
Optional arguments:
-h|--help               Show this help message and exit.

-v|--verbose		Set to see processing messages
-r|--region FASTA file	FASTA file that contains the region you want to call,
			designed to accept UCSC "Genomic Sequence Near Gene" FASTA files
			which allow for filtering on exon/intron. (default: None)

-o|--outdir		Directory where the output of the variant calling is stored
			(default: current working directory)

Exit codes:
0			Operation complete
1			Syntax error
128			Invalid argument

'
}

##### Main

# Command line processor

while [ "$1" != "" ]; do
    case $1 in
	-b | --bamfolder )	shift
				bamfolder=$1
				;;
	-g | --genome )		shift
				genome=$1
				;;
	-a | --annotation )	shift
				annotation=$1
				;;
	-r | --region )		shift
				region=$1
				;;
	-o | --outdir )		shift
				outdir=$1
				;;
	-v | --verbose )	verbose=1
				;;
	-h | --help )		usage
				exit 0
				;;
	* )			usage
				exit 1
    esac
    shift
done

# Check required arguments are given
if [ "$bamfolder" == "" ]; then
  echo "Specify a folder with bam files, sorted by coordinate"
  exit 128
fi

if [ "$genome" == "" ]; then
  echo "Specify the genome that was used for mapping the bam files"
  exit 128
fi

if [ "$annotation" == "" ]; then
  echo "Specify an annotation VCF file that holds the variants you want to annotate"
  exit 128
fi

if [ $verbose = 1 ]; then
  echo "Processing bam files from $bamfolder
that were mapped to $genome
selecting region based on $region
using $annotation for annotation
writing output to $outdir
";
fi

# Wait for user input if the output directory already exists
if [ -f "$outdir/*.vcf" ]; then
  read -n 1 -s -r -p 'The output directory you indicated already exists and has vcf files in it.
Please note that files may be overwritten and results may be messed up.
IT IS HIGHLY RECOMMENDED THAT YOU ABORT
Press any key to continue or Ctrl-c to abort...' ;
fi

# make the output directory
mkdir -p $outdir

# grab coordinates from the region FASTA and write to BED formatted file
egrep "chr" $region | cut -f 2 -d " "  \
               | cut -f 2 -d "=" \
               | tr ":" " " \
               | tr "-" " " \
               | sort  \
	       > $outdir/$(basename $region .fasta).bed

## Create and filter VCF files

for file in $bamfolder/*.bam ; do
  # If index not present, index bam file
  if [ ! -f "$file.bai" ]; then
    if [ $verbose = 1 ]; then echo "Indexing $file"; fi
    samtools index $file ;
  fi
  if [ $verbose = 1 ]; then echo "Writing vcf files for $file"; fi
  # Variant calling for the input BAM files on the target region, written to vcf format
  freebayes -f $genome -t $outdir/$(basename $region .fasta).bed \
  -v $outdir/$(basename $file Aligned.sortedByCoord.out.bam).vcf  $file &
done
wait

# Filtering on quality score, read depth and number of reads supporting the variant
for file in $outdir/*.vcf ; do
  if [ $verbose = 1 ]; then echo "Filtering $file on quality score and depth of variant bases "; fi
  bcftools view -i 'QUAL > 50 && INFO/DP > 1 && (SAF+SAR) > (INFO/DP*0.3)' \
  $file > $outdir/$(basename $file .vcf).filtered.vcf &
done
wait

# Change chromosome naming (chr1 -> 1) and zip the filtered vcf files
for file in $outdir/*.filtered.vcf ; do
  if [ $verbose = 1 ]; then echo "Formatting, zipping and indexing $file"; fi
  sed 's/chr//' $file > $outdir/$(basename $file .filtered.vcf)_no_chr.vcf
  cat $outdir/$(basename $file .filtered.vcf)_no_chr.vcf > $file
  bgzip -f $file
  bcftools index $file.gz
  rm $outdir/$(basename $file .filtered.vcf)_no_chr.vcf &
done

# Annotate the vcf files with ClinVar data on the variant
for file in $outdir/*filtered.vcf.gz ; do
  if [ $verbose = 1 ]; then echo "Annotating $file"; fi
  bcftools annotate --threads 32 -c RS,CLNDN,CLNHGVS,CLNSIG,MC \
  -a $annotation $file \
  > $outdir/$(basename $file .filtered.vcf.gz).output.vcf &
done
wait

## Generate summary

# write header
echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tREFRD\tALTRD\tRS\tCLNDN\tGENPOS\tCLINSIG\tTYPE\tSAMPLE" > $outdir/vcf_summary_$(basename $region .fasta).txt

# write vcf entries
for file in $outdir/*.output.vcf; do
  if [ $verbose = 1 ]; then echo "Adding records of $file to summary "; fi
  bcftools query  \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%SRF\t%SRR\t(%SRF+%SRR)\t%SAF\t%SAR\t(%SAF+%SAR)\t%RS\t%CLNDN\t%CLNHGVS\t%CLNSIG\t%MC\n' \
  $file | awk '{print $1,$2,$3,$4,$5,$6,$7+$8,$10+$11,$13,$14,$15,$16,$17,"FILENAME"}' OFS='\t' | sed "s/FILENAME/$(basename $file .output.vcf)/"  -  \
  >> $outdir/vcf_summary_$(basename $region .fasta).txt;
done

# Sort vcf entries (not header) on position
if [ $verbose = 1 ]; then echo "Sorting summary"; fi
cat $outdir/vcf_summary_$(basename $region .fasta).txt | awk 'NR<2{print;next}{print | "sort -k 2"}' > $outdir/sorted.temp
cat $outdir/sorted.temp > $outdir/vcf_summary_$(basename $region .fasta).txt
rm $outdir/sorted.temp

if [ $verbose = 1 ]; then echo "Done"; fi

exit 0
