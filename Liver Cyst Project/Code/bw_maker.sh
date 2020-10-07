#!/bin/bash

##### Variables

bamfile=
outdir=$(pwd)/bw_files/
alignmentType=paired
separateFWRV=0
genomeSize=
calcSize=
genomeVersion=
genome=

##### Functions

usage()
{
  echo '
Usage: UCSC_tracklayer.sh [options]

This script generates bigWig files from input bam files.

Required arguments:
-b|--bamfile             Bam file that will be used to make (a) bigwig(s).

--genomeVersion          Version of the genome used to create the bam file.

Optional arguments:
-h|--help                Show this help message and exit.

-v|--verbose             Set to see processing messages. (default: not set)

-o|--outdir              Directory where the bigwig files will be saved.
                         (default: current_directory/bw_files/)

-a|--alignmentType       By which type of alignment (single- or paired-end
                         the bam files were created. (default: paired)

-s|--separateRFRV        Forward and reverse strands are separated into two
                         .bw files if set. (default: not set)

-g|--genomeSize          Effective mappable size of the genome used to create
                         the bam file, needed for correct normalisation.
                         (default: None)

-c | --calcSize          Calculates the effective genome size if set.
                         Requires a file be provided with genome. (default: not set)

--genome                 Filepath to the genome used to create the bam file.
                         Required if calcSize is set. (default: None)

Exit codes:
0                       Operation complete
1                       Syntax error
128                     Invalid argument

'
}

##### Main

# Command line processor
while [ "$1" != "" ]; do
    case $1 in
        -b | --bamfile )      shift
                                bamfile=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
        -a | --alignmentType )  shift
                                alignmentType=$1
                                ;;
        -s | --separateRFRV )   separateFWRV=1
                                ;;
        -g | --genomeSize )     shift
                                genomeSize=$1
                                ;;
	--genomeVersion )	shift
				genomeVersion=$1
				;;
        -c | --calcSize )       calcSize=1
                                ;;
        --genome )              shift
                                genome=$1
                                ;;
        -v | --verbose )        verbose=1
                                ;;
        -h | --help )           usage
                                exit 0
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

## Checks
if [ "$bamfile" == "" ] || [ "$genomeVersion" == "" ]; then
  echo "
You must provide a bam file and the version of the genome used
to create the file (e.g. hg19, hg38).
"
  exit 128
fi

# If outdir not exist, make outdir
if [ ! -d "$outdir" ]; then
  echo "Creating outdir"
  mkdir -p $outdir
fi

# Check if calcSize is on when providing genome
if [ "$calcSize" == "" ] && [ "$genome" != "" ]; then
  echo "
You have provided a genome file, but have not set calcSize.
Assuming you want to calculate the effective genome size, calcSize has been set to true.
"
  calcSize=1
fi

# Check genome is calcSize is on
if [ "$calcSize" == 1 ] && [ "$genome" == "" ]; then
  echo "
You have set calcSize, but have not provided a genome file.
Please provide the filepath to the genome fasta file used for creating the input bamfile.
"
  exit 128
fi

if [ "$alignmentType" != "single" ] && [ "$alignmentType" != "paired" ]; then
  echo "Please specify alignmentType 'single' or 'paired'"
  exit 128 ;
fi

## Bigwig generation
# Index input file
samtools index $bamfile &

# Get chromosome sizes
if [ ! -f "$outdir/$genomeVersion.chrom.sizes" ]; then
  fetchChromSizes $genomeVersion > $outdir/$genomeVersion.chrom.sizes &
fi


# Calculate the effective genome size from the readlength of the bamfile
# using the genome the bamfiles were mapped to
if [ "$calcSize" == 1 ]; then
  if [ "$verbose" == 1 ]; then echo "Calculating effective genome size"; fi
  readlength=$(samtools view $bamfile | head  -n 1| cut -f 10 | wc -c)
  genomeSize=$(unique-kmers.py -q -k $readlength $genome 2>&1 > /dev/null | tail -n 1 | cut -f 7 -d " ")
  echo "Effective genome size: $genomeSize"
fi
wait

# Make bedgraphs for further processing, adding options where required
name=$(basename $bamfile Aligned.sortedByCoord.out.bam).bedgraph

if [ "$verbose" == 1 ]; then echo "Creating bedgraph(s)"; fi

bamCoverage -b $bamfile -o $outdir/$(if [ "$separateFWRV" == 1 ]; then echo "FW_"; fi)$name \
-of bedgraph -p max --effectiveGenomeSize $genomeSize --normalizeUsing BPM --skipNonCoveredRegions \
$(if [ "$alignmentType" == "paired" ];then echo "--samFlagInclude 64"; fi) \
$(if [ "$separateFWRV" == 1 ]; then echo "--filterRNAstrand forward"; fi) &

# If separateFWRV is set, run the same command as above again, but filter the reverse strand
  if [ "$separateFWRV" == 1 ]; then
    bamCoverage -b $bamfile -o $outdir/RV_$name \
    -of bedgraph -p max --effectiveGenomeSize $genomeSize --normalizeUsing BPM --skipNonCoveredRegions \
    $(if [ "$alignmentType" == "paired" ];then echo "--samFlagInclude 64"; fi) --filterRNAstrand reverse &
  fi ;
wait

# Process the bedgraph files into the final bigwig
if [ "$separateFWRV" == 1 ]; then
  # make the REV strands into the inverse
  if [ "$verbose" == 1 ]; then echo "Inverting reverse strand" ; fi
  awk '{print $1, $2, $3, $4 * -1}' OFS='\t' $outdir/RV_$name > \
  $outdir/$(basename RV_$name .bedgraph).tmp
  cat $outdir/$(basename RV_$name .bedgraph).tmp > $outdir/$(basename RV_$name)
  rm $outdir/$(basename RV_$name .bedgraph).tmp

  if [ "$verbose" == 1 ]; then echo "Clipping and sorting"; fi
  for file in {$outdir/FW_$name,$outdir/RV_$name}; do
    # Clip records not mapping to chromosomes
    awk '$1 ~ /chr/ {print}' $file > $outdir/$(basename $file .bedgraph).clipped.bedgraph
    # Sort clipped bedgraphs
    bedSort $outdir/$(basename $file .bedgraph).clipped.bedgraph $outdir/$(basename $file)
    rm $outdir/$(basename $file .bedgraph).clipped.bedgraph &
  done
  wait

  if [ "$verbose" == 1 ]; then echo "Converting to bw" ; fi
  for file in {$outdir/FW_$name,$outdir/RV_$name} ;do
    # Convert bedgraph to bigwig
    bedGraphToBigWig $file $outdir/$genomeVersion.chrom.sizes \
    $outdir/$(basename $file .bedgraph).bw &
  done
  wait

  rm $outdir/FW_$name
  rm $outdir/RV_$name
fi

# Process the bedgraph file into the final bigwig
if [ "$separateFWRV" == 0 ]; then
  if [ "$verbose" == 1 ]; then echo "Clipping and sorting"; fi
  # Clip records not mapping to chromosomes
  awk '$1 ~ /chr/ {print}' $outdir/$name > $outdir/$(basename $name .bedgraph).clipped.bedgraph
  # Sort clipped bedgraphs
  bedSort $outdir/$(basename $name .bedgraph).clipped.bedgraph $outdir/$(basename $name)
  rm $outdir/$(basename $name .bedgraph).clipped.bedgraph

  if [ "$verbose" == 1 ]; then echo "Converting to bw" ; fi
  # Convert bedgraph to bigwig
  bedGraphToBigWig $outdir/$(basename $name) $outdir/$genomeVersion.chrom.sizes \
  $outdir/$name.bw

  rm $outdir/$(basename $name)
fi
exit 0
