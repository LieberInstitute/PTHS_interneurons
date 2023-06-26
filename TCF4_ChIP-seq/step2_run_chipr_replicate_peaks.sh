PROJDIR=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq
TMP=/scratch/bnphan
conda activate chipr

################################
## 1) Forrest et al. 2018 SH-SY5Y TCF4 ChIP-seq
PROJ=Forrest2018_ChIPseq
cd ${PROJDIR}/${PROJ}
mkdir -p ${PROJDIR}/${PROJ}/chipr

## copy over the replicate chipseq peaks
cp ${PROJDIR}/${PROJ}/peak/*/*300K.regionPeak.gz ${PROJDIR}/${PROJ}/chipr
gunzip --force ${PROJDIR}/${PROJ}/chipr/*300K.regionPeak.gz
INPUT2=$(ls ${PROJDIR}/${PROJ}/chipr/*300K.regionPeak | tr '\n' ' ')

# run chipR
chipr -i $INPUT2 -m 2 -o chipr/${PROJ}.chipr

############################################
## 2) Xia et al. 2018 SH-SY5Y TCF4 ChIP-seq
## this dataset is bad. almost no IDR peaks
PROJ=Xia2018_ChIPseq
cd ${PROJDIR}/${PROJ} && mkdir -p ${PROJDIR}/${PROJ}/chipr

## copy over the replicate chipseq peaks
cp ${PROJDIR}/${PROJ}/peak/*/*300K.regionPeak.gz ${PROJDIR}/${PROJ}/chipr
gunzip --force ${PROJDIR}/${PROJ}/chipr/*300K.regionPeak.gz
INPUT2=$(ls ${PROJDIR}/${PROJ}/chipr/*300K.regionPeak | tr '\n' ' ')

# run chipR
chipr -i $INPUT2 -m 2 -o chipr/${PROJ}.chipr


#######################################
## 3) ENCODE 2021 SK-N-SH TCF4 ChIP-seq
PROJ=ENCODE2021_ChIPseq
cd ${PROJDIR}/${PROJ} && mkdir -p ${PROJDIR}/${PROJ}/chipr

## copy over the replicate chipseq peaks
cp ${PROJDIR}/${PROJ}/peak/*/*300K.regionPeak.gz ${PROJDIR}/${PROJ}/chipr
gunzip --force ${PROJDIR}/${PROJ}/chipr/*300K.regionPeak.gz
INPUT2=$(ls ${PROJDIR}/${PROJ}/chipr/*300K.regionPeak | tr '\n' ' ')

# run chipR
chipr -i $INPUT2 -m 2 -o chipr/${PROJ}.chipr




