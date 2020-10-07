#!/bin/bash

##### Variables
indexdir=
read1=
read2=
mode="star"
outdir=
genome=
annotation=
kmer=31
libtype="A"
verbose=0
mapmode=
cores=$(nproc)

##### Functions
usage()
{
  echo '
Usage: mapper.sh [options]... [FILE]

This is a wrapper script for using STAR and salmon mapping algorithms. It is able to generate the needed indices and
can map fastq files both single end and paired end. If requested the scripts can also provide FASTQC output along with
the mapping results.

This scripts is made with the intention of using gencode Fasta and annotation files. You can obtain these files from:
https://www.gencodegenes.org/human/

Fastq files need to provided in gzipped (.gz) format, Fasta and GTF annotation files need to be unzipped.

Arguments:
-h|--help               Show this help message and exit.

-v|--verbose            Set to see processing messages. (default: not set)

-i|--indexdir           Directory where the index is stored or generated if mode is set to index.
                        (default: current_working_directory/star_index or current_working_directory/salmon_index)

-1|--read1              Gzipped fastq file with reads. Must contain the forwards reads if
                        mapping mate-pair sequenced reads. Mapping will be single-end
                        if -2 is not set. Required for mapping modes (default: None)

-2|--read2		Gzipped fastq file with reads from the second mate. Mapping will
			be paired-end if this is set. (default: None)

-m|--mode		Mode the script will operate in. Options are "star" "star_index"
			"salmon" and "salmon_index". (default: star)

-o|--outdir             Directory where the output of the mapping algorithm is stored.
                        (default: current_working_directory/star_output or current_working_directory/salmon_output)

-g|--genome		Genome or transcriptome fasta file. Genome file required for star_index mode.
			Transcriptome file required for salmon_index mode. (default: None)

-a|--annotation		GTF annotation file. Required for star_index mode. (default: None)

-k|--kmer		Minimum acceptable length for matches. Used when building the salmon index.
			Must always be an odd number. Around 0.5 * read length is recommended (default: 31)

-l|--libtype		The orientation that the mate pair is in. Only for paired-end mapping with salmon.
			See salmon documentation for more information. (default: A)

-c|--cores		The amount of cores the operation should run on. (default: all cores)

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
        -i | --indexdir )	shift
                                indexdir=$1
                                ;;
        -1 | --read1 )          shift
                                read1=$1
                                ;;
        -2 | --read2 )     	shift
                                read2=$1
                                ;;
        -m | --mode )  	   	shift
                                mode=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
	-g | --genome )		shift
				genome=$1
				;;
	-a | --annotation )	shift
				annotation=$1
				;;
	-k | --kmer )		shift
				kmer=$1
				;;
	-l | --libtype )	shift
				libtype=$1
				;;
	-c | --cores )		shift
				cores=$1
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

## Check if arguments are valid

# Mode set within parameters?
if [ "$mode" != "star" ] && [ "$mode" != "star_index" ] && [ "$mode" != "salmon" ] && [ "$mode" != "salmon_index" ]; then
  echo "Please set mode to one of the following options:
star
star_index
salmon
salmon_index
"
  exit 128
fi

# Arguments for star index present?
if [ "$mode" == "star_index" ]; then
  if [ "$genome" == "" ] || [ "$genome" == "*.gz" ] ; then
    echo "Please provide an unzipped genome fasta for creating the star index"
    exit 128
  fi
  if [ "$annotation" == "" ] || [ "$annotation" == "*.gz" ]; then
    echo "Please provide an unzipped GTF annotation file for creating the star index"
    exit 128
  fi
fi

# Arguments for salmon index present?
if [ "$mode" == "salmon_index" ]; then
  if [ "$genome" == "" ] || [ "$genome" == "*.gz" ] ; then
    echo "Please provide an unzipped transcriptome fasta for creating the salmon index"
    exit 128
  fi
fi

# Mapping type?
if [ "$read1" != "" ] && [ "$read2" == "" ]; then
  mapmode=1
elif [ "$read1" != "" ] && [ "$read2" != "" ]; then
  mapmode=2
fi

# Arguments for star mapping present?
if [ "$read1" == "" ] && [ "$mode" == "star" ]; then
  echo "Please provide a gzipped fastq file for star to map."
  exit 128
elif [ "$read1" == "*.fastq" ] && [ "$mode" == "star" ]; then
  echo "Please provide a gzipped fastq file for star to map."
  exit 128
elif [ "$indexdir" == "" ] && [ "$mode" == "star" ]; then
  echo "Please provide and index directory. You can generate an index using --mode star_index"
  exit 128
fi

# Arguments for salmon mapping present?
if [ "$read1" == "" ] && [ "$mode" == "salmon" ]; then
  echo "Please provide a gzipped fastq file for salmon to map."
  exit 128
elif [ "$read1" == "*.fastq" ] && [ "$mode" == "salmon" ]; then
  echo "Please provide a gzipped fastq file for salmon to map."
  exit 128
elif [ "$indexdir" == "" ] && [ "$mode" == "salmon" ]; then
  echo "Please provide and index directory. You can generate an index using --mode salmon_index"
  exit 128
fi


## Set default arguments

# Set output directory if none is provided
if [ "$mode" == "star" ] && [ "$outdir" == "" ]; then
  outdir=$(pwd)/star_output
elif [ "$mode" == "salmon" ] && [ "$outdir" == "" ]; then
  outdir=$(pwd)/salmon_output
fi
outdir=$outdir/

# Set index directory if none is provided
if [ "$mode" == "star_index" ] && [ "$indexdir" == "" ]; then
  indexdir=$(pwd)/star_index
elif [ "$mode" == "salmon_index" ] && [ "$indexdir" == "" ]; then
  indexdir=$(pwd)/salmon_index
fi
indexdir=$indexdir/

## Indexing STAR
if [ "$mode" == "star_index" ]; then
  mkdir -p $indexdir
  STAR --runMode genomeGenerate \
  --runThreadN $cores --genomeDir $indexdir --genomeFastaFiles $genome  --sjdbGTFfile $annotation
fi

## Indexing Salmon
if [ "$mode" == "salmon_index" ]; then
  mkdir -p $indexdir
  salmon index --type quasi --perfectHash --gencode \
  -p $cores -i $indexdir -t $genome -k $kmer
fi

## Mapping STAR
if [ "$mode" == "star" ] ; then
  mkdir -p $outdir
  echo "Starting mapping for $read1 $read2"
  STAR --quantMode GeneCounts --twopassMode Basic --alignSoftClipAtReferenceEnds No \
  --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimScoreMin 1 --outFilterIntronMotifs RemoveNoncanonical \
  --outFilterType BySJout --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat\
  --runThreadN $cores --genomeDir $indexdir --readFilesIn $read1 $read2 --outFileNamePrefix $outdir/$(basename $read1 .fastq.gz)
fi

## Mapping Salmon
if [ "$mode" == "salmon" ] && [ "$mapmode" == 1 ]; then
  mkdir -p $outdir
  salmon quant --validateMappings --rangeFactorizationBins 4 \
  -p $cores -i $indexdir -l $libtype -r $read1 -o $outdir/$(basename $read1 .fastq.gz)
elif [ "$mode" == "salmon" ] && [ "$mapmode" == 2 ]; then
  mkdir -p $outdir
  salmon quant --validateMappings --rangeFactorizationBins 4 \
  -p $cores -i $indexdir -l $libtype -1 $read1 -2 $read2 -o $outdir/$(basename $read1 _R1.fastq.gz)
fi

exit 0
