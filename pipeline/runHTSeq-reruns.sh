#!/bin/bash -l

#set -ex

mail="nathaniel.street@umu.se"
in=/mnt/picea/projects/aspseq/24_Bellini_cambium/BelliniCambium/star
out=/mnt/picea/projects/aspseq/24_Bellini_cambium/BelliniCambium/htseq
HTSeq=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3

module load bioinfo-tools htseq

for f in `find $in -name "*merged.bam"`;
do
    bash $UPSCb/pipeline/runHTSeq.sh $out $f $HTSeq
done
