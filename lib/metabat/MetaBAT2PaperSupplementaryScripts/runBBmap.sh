#!/bin/bash

bbmap=$(which bbmap.sh)

set -x
set -e

if [ $# -lt 2 ]
then
  echo "USAGE: $0 reference reads.fastq [read2.fastq]" 1>&2
  exit 1
fi

ref=$(realpath $1)
refdir=$ref.d
in1=$(realpath $2)
if [ ! -f "$in1" ]
then
  echo "Could not find $in1" 1>&2
  exit 1
fi

interleaved=auto
in2=
if [ -f "$3" ]
then
  in2=$(realpath $3)
  if [ ! -f "$in2" ]
  then
     echo "Could not find $in2" 1>&2
     exit 1
  fi
else
  # detect interleaved
  read1=
  read2=
  if [ "$in1" == "${in1%gz}" ]
  then 
    read1=$(head -1 $in1 | awk '{print $1}' | sed 's/\/[12]//;')
    read2=$(head -5 $in1 | tail -1 | awk '{print $1}' | sed 's/\/[12]//;')
  else
    read1=$(gunzip -c $in1 | head -1 | awk '{print $1}' | sed 's/\/[12]//;')
    read2=$(gunzip -c $in1 | head -5 | tail -1 | awk '{print $1}' | sed 's/\/[12]//;')
  fi
  if [ "$read1" == "$read2" ]
  then
    interleaved=true
    echo "Detected interleaved fastq" 
  fi

  in2=
fi

bbopts=${bbopts:="samversion=1.4 local=t kbp=f minhits=2 minratio=0.8 maxindel=50 mdtag=true requirecorrectstrand=false trd=t interleaved=$interleaved usemodulo=t "}

if [ -n "${idfilter}" ]
then
  bbopts="${bbopts} idfilter=${idfilter}"
fi


mkdir -p $refdir

inputs="in=$in1"
if [ -f "$in2" ]
then
  # split reads
  inputs="$inputs in2=$in2"
fi


out=${ref##*/}-${in1##*/}.sam
outu="outu=$out.unmap.fastq.gz"
if [ "${outputunmapped}" == 'f' ]
then
  outu="outputunmapped=${outputunmapped}"
fi

  (cd $refdir && [ -d ref/genome ] ||  $bbmap $bbopts ref=$ref build=1)
  tmpout=$TMPDIR/${out##*/}
  ( cd $refdir && \
     $bbmap $bbopts $inputs build=1 outm=$tmpout.tmp.sam $outu && \
     samtools view -Sbu $tmpout.tmp.sam | \
       samtools sort -l 0 -m 1G -@ 12 -T $tmpout.dummy | \
         samtools calmd -u - $ref | \
           samtools view -b -@ 12 - > $out.tmp.bam && \
     mv $out.tmp.bam $out.bam  && \
     samtools index $out.bam && \
     rm $tmpout.tmp.sam || echo "sorting and indexing $sam failed" 1>&2
 )


