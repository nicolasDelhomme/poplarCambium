#!/bin/bash -l

#set -ex

mail="nathaniel.street@umu.se"
in=/mnt/picea/projects/aspseq/24_Bellini_cambium/BelliniCambium/raw
out=/mnt/picea/projects/aspseq/24_Bellini_cambium/BelliniCambium
genome=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/indices/STAR/Potri03
gff3=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_gene_exons.gff3
HTSeq=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3
start=1
end=8

module load bioinfo-tools FastQC samtools star htseq sortmerna

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fq.gz"`; do echo "${f//_[1,2].fq.gz/}" ; done | sort | uniq | while read line;
do
    bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome -G $gff3 -H $HTSeq BelliniCambium $mail ${line}_1.fq.gz ${line}_2.fq.gz $out
done
