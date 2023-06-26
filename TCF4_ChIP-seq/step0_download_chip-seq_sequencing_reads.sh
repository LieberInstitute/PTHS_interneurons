PROJDIR=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq
TMP=/scratch/bnphan/fastq
mkdir -p $TMP

cd $PROJDIR

############################################################
## 1) grab Forrest et al 2018 TCF4 ChIP-seq in SH-SY5Y cells

geofetch -i GSE96915 \
--name "Forrest2018_ChIPseq" \
--metadata-root $PROJDIR \
--add-convert-modifier \
--just-metadata \
--refresh-metadata \
--discard-soft

SRARAW=/scratch/bnphan/ncbi/sra
SRAFQ=$PROJDIR/Forrest2018_ChIPseq/fastq
mkdir -p $SRAFQ $TMP

for i in {2..5}; 
do 
LINE=$(awk -v var=$i 'NR==var{print $0}' ${PROJDIR}/Forrest2018_ChIPseq/GSE96915_PEP/GSE96915_PEP_raw.csv)
## slot 28 is where the SRR numbers are
SRR=$(echo $LINE|cut -d',' -f28) 
echo $SRR
## slopt 40 where the human readable names are
SAMPLE=$(echo $LINE|cut -d',' -f40)
## go through and make fastq
if [[ ! -s $SRAFQ/$SRR.fastq.gz ]]; then
    fasterq-dump $SRARAW/$SRR.sra -O $TMP --temp $TMP -p
    ## gzip the files manually to conserve space
    gzip $TMP/$SRR*
    ## move to final destination
    rsync -Pav --remove-source-files $TMP/$SRR* $SRAFQ
fi
done

############################################################
## 2) grab Xia et al 2018 TCF4 ChIP-seq in SH-SY5Y cells

PROJ=Xia2018_ChIPseq
geofetch -i GSE112704 \
--name $PROJ \
--metadata-root $PROJDIR \
--add-convert-modifier \
--refresh-metadata \
--discard-soft

SRAFQ=$PROJDIR/$PROJ/fastq
TMP=/scratch/bnphan/fastq
mkdir -p $SRAFQ $TMP

for i in {2..14}; 
do 
LINE=$(awk -v var=$i 'NR==var{print $0}' ${PROJDIR}/$PROJ/GSE112704_PEP/GSE112704_PEP_raw_subtable.csv)
## slot 3 is where the SRR numbers are
SRR=$(echo $LINE|cut -d',' -f3 | tr -d "[:blank:]" ) 
echo $SRR
## go through and make fastq
if [[ ! -s $SRAFQ/$SRR.fastq.gz ]]; then
    fasterq-dump $SRARAW/$SRR.sra -O $TMP --temp $TMP -p
    ## gzip the files manually to conserve space
    gzip $TMP/$SRR*
    ## move to final destination
    rsync -Pav --remove-source-files $TMP/$SRR* $SRAFQ
fi
done


############################################################
## 3) grab ENCODE v3 TCF4-GFP ChIP-seq in SK-N-SH cells

PROJ=ENCODE2021_ChIPseq
geofetch -i GSE177840 \
--name $PROJ \
--metadata-root $PROJDIR \
--add-convert-modifier \
--discard-soft

geofetch -i GSE176975 \
--name $PROJ \
--metadata-root $PROJDIR \
--add-convert-modifier \
--discard-soft

SRABAM=$PROJDIR/$PROJ/bam
mkdir -p $SRABAM 

## download aligned bam files for GFP-tagged TCF4 ChIP-seq
wget -P $TMP --tries=0 https://www.encodeproject.org/files/ENCFF752LTA/@@download/ENCFF752LTA.bam
wget -P $TMP --tries=0 https://www.encodeproject.org/files/ENCFF146HUB/@@download/ENCFF146HUB.bam

## download aligned bam files for controls of GFP-TCF4 
wget -P $TMP --tries=0 https://www.encodeproject.org/files/ENCFF134IOI/@@download/ENCFF134IOI.bam
wget -P $TMP --tries=0 https://www.encodeproject.org/files/ENCFF981WVH/@@download/ENCFF981WVH.bam
wget -P $TMP --tries=0 https://www.encodeproject.org/files/ENCFF160ORG/@@download/ENCFF160ORG.bam

rsync -Pav --remove-source-files $TMP/ENCFF* $SRABAM

############################################################
## 4) transfer the ChIP-seq results from Rannals et al. 2016
PROJ=Rannals2016_ChIPseq
SRAFQ=$PROJDIR/$PROJ/fastq
SRABAM=$PROJDIR/$PROJ/bam
mkdir -p $SRAFQ $SRABAM

## send over the fastq files
rsync -Pav bphan@jhpce-transfer01.jhsph.edu:/dcs04/legacy-dcs01-ajaffe/Brady/ChIPseq/TCF4_core/FASTQ/* $SRAFQ

## send over the aligned BAM files
rsync -Pav bphan@jhpce-transfer01.jhsph.edu:/dcs04/legacy-dcs01-ajaffe/Brady/ChIPseq/TCF4_core/BAM/* $SRABAM







################################################################
## 1) grab Wang et al TCF4 ChIP-seq in human MGE-like organoids

geofetch -i GSE235071 \
--name "Wang2022_ChIPseq" \
--metadata-root $PROJDIR \
--add-convert-modifier \
--just-metadata \
--refresh-metadata \
--discard-soft

SRARAW=/scratch/bnphan/ncbi/sra
SRAFQ=$PROJDIR/Wang2022_ChIPseq/fastq
mkdir -p $SRAFQ $TMP

for i in {2..3}; 
do 
LINE=$(awk -v var=$i 'NR==var{print $0}' ${PROJDIR}/Wang2022_ChIPseq/GSE235071_PEP/GSE235071_PEP_raw.csv)
## slot 28 is where the SRR numbers are
SRR=$(echo $LINE|cut -d',' -f6) 
echo $SRR
## slopt 40 where the human readable names are
SAMPLE=$(echo $LINE|cut -d',' -f1)
## go through and make fastq
if [[ ! -s $SRAFQ/$SRR.fastq.gz ]]; then
    prefetch $SRR
    fasterq-dump $SRARAW/$SRR.sra -O $TMP --temp $TMP -p
    ## gzip the files manually to conserve space
    gzip $TMP/$SRR*
    ## move to final destination
    rsync -Pav --remove-source-files $TMP/$SRR* $SRAFQ
fi
done

