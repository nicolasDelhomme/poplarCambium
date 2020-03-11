#!/bin/bash -l

#set -ex
proj=u2015056
mail="nicolas.delhomme@umu.se"
in=/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium/raw
out=/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium
genome=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/indices/STAR/Potri03
gtf=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gtf/Ptrichocarpa_v3.0_210_gene_exons.gtf
HTSeq=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3
start=7
end=8

module load bioinfo-tools FastQC Trimmomatic samtools star/2.4.0f1 htseq sortmerna

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

line=$in/"FCC4AHUACXX-TRDPOPwhfTASRAAPEI-220_L1"
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome -G $gtf -H $HTSeq $proj $mail ${line}_1.fq.gz ${line}_2.fq.gz $out
