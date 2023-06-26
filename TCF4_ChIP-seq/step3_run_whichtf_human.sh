PROJDIR=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq
TMP=/scratch/bnphan
WTF_HOME=/home/bnphan/src/whichtf
conda activate whichtf

##############################################
## 1) Forrest et al. 2018 SH-SY5Y TCF4 ChIP-seq
PROJ=Forrest2018_ChIPseq
cd ${PROJDIR}/${PROJ}

mkdir -p ${PROJDIR}/${PROJ}/whichtf

## WhichTF on ENCODE IDR peaks
PEAK=peak/idr_reproducibility/idr.optimal_peak.regionPeak.gz
TMP=whichtf/${PROJ}.tmp.narrowPeak
OUT=whichtf/${PROJ}.idr.optimal_peak.whichtf.tsv
gunzip < $PEAK > $TMP
${WTF_HOME}/bin/WhichTF $TMP hg38 --outFile ${OUT}

## WhichTF on ChIP-R peaks
CHIPR=chipr/${PROJ}.chipr_optimal.bed
OUT2=whichtf/${PROJ}.chipr.optimal_peak.whichtf.tsv
${WTF_HOME}/bin/WhichTF $CHIPR hg38 --outFile ${OUT2}


#######################################
## 3) ENCODE 2021 SK-N-SH TCF4 ChIP-seq
PROJ=ENCODE2021_ChIPseq
cd ${PROJDIR}/${PROJ}

## WhichTF on ENCODE IDR peaks
PEAK=peak/idr_reproducibility/idr.optimal_peak.regionPeak.gz
TMP=whichtf/${PROJ}.tmp.narrowPeak
OUT=whichtf/${PROJ}.idr.optimal_peak.whichtf.tsv
gunzip < $PEAK > $TMP
${WTF_HOME}/bin/WhichTF $TMP hg38 --outFile ${OUT}

## WhichTF on ChIP-R peaks
CHIPR=chipr/${PROJ}.chipr_optimal.bed
OUT2=whichtf/${PROJ}.chipr.optimal_peak.whichtf.tsv
${WTF_HOME}/bin/WhichTF $CHIPR hg38 --outFile ${OUT2}


#######################################
## 3) Xia 2018 SH-SY5Y TCF4 ChIP-seq
PROJ=Xia2018_ChIPseq
cd ${PROJDIR}/${PROJ}

## WhichTF on ENCODE IDR peaks
PEAK=peak/idr_reproducibility/idr.optimal_peak.regionPeak.gz
TMP=whichtf/${PROJ}.tmp.narrowPeak
OUT=whichtf/${PROJ}.idr.optimal_peak.whichtf.tsv
gunzip < $PEAK > $TMP
${WTF_HOME}/bin/WhichTF $TMP hg38 --outFile ${OUT}

## WhichTF on ChIP-R peaks
CHIPR=chipr/${PROJ}.chipr_optimal.bed
OUT2=whichtf/${PROJ}.chipr.optimal_peak.whichtf.tsv
${WTF_HOME}/bin/WhichTF $CHIPR hg38 --outFile ${OUT2}





