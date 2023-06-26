PROJDIR=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq
TMP=/scratch/bnphan
conda activate chipseq

################################
## 1) Forrest et al. 2018 SH-SY5Y TCF4 ChIP-seq
PROJ=Forrest2018_ChIPseq
mkdir -p $TMP/$PROJ && cd $TMP/$PROJ
rsync -Pav $PROJDIR/$PROJ/* $TMP/$PROJ

INPUT_JSON=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq/${PROJ}/${PROJ}.json
caper run /home/bnphan/src/chip-seq-pipeline2/chip.wdl \
-i "${INPUT_JSON}" \
--singularity /home/bnphan/src/chip-seq-pipeline2/chip-seq-pipeline_v2.2.1.sif \
--max-concurrent-tasks 1 --db in-memory

META_JSON=/scratch/bnphan/Forrest2018_ChIPseq/chip/d91d42aa-b60a-421b-853a-f2e05780005a/metadata.json

## copy over the pipeline outputs
croo --tmp-dir $TMP --out-dir $PROJDIR/$PROJ \
--method copy $META_JSON
rm -r $TMP/$PROJ


############################################
## 2) Xia et al. 2018 SH-SY5Y TCF4 ChIP-seq
## this dataset is bad. almost no IDR peaks
PROJ=Xia2018_ChIPseq
mkdir -p $TMP/$PROJ && cd $TMP/$PROJ
rsync -Pav $PROJDIR/$PROJ/* $TMP/$PROJ

INPUT_JSON=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq/${PROJ}/${PROJ}.json
caper run /home/bnphan/src/chip-seq-pipeline2/chip.wdl \
-i "${INPUT_JSON}" --singularity /home/bnphan/src/chip-seq-pipeline2/chip-seq-pipeline_v2.2.1.sif \
--max-concurrent-tasks 1 --db in-memory

META_JSON=/scratch/bnphan/Xia2018_ChIPseq/chip/383895ee-8055-44ce-ae8f-024cdb4d5fe4/metadata.json

## copy over the pipeline outputs
croo --tmp-dir $TMP --out-dir $PROJDIR/$PROJ \
--method copy $META_JSON
rm -r $TMP/$PROJ


#######################################
## 3) ENCODE 2021 SK-N-SH TCF4 ChIP-seq
PROJ=ENCODE2021_ChIPseq
mkdir -p $TMP/$PROJ && cd $TMP/$PROJ
rsync -Pav $PROJDIR/$PROJ/* $TMP/$PROJ

INPUT_JSON=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq/${PROJ}/${PROJ}.json
/home/bnphan/src/caper/bin/caper run /home/bnphan/src/chip-seq-pipeline2/chip.wdl \
-i "${INPUT_JSON}" --singularity /home/bnphan/src/chip-seq-pipeline2/chip-seq-pipeline_v2.2.1.sif \
--max-concurrent-tasks 1

META_JSON=/scratch/bnphan/ENCODE2021_ChIPseq/chip/a626ec8c-395d-4be7-a686-39cf3972bc8b/metadata.json


## copy over the pipeline outputs
croo --tmp-dir $TMP --out-dir $PROJDIR/$PROJ \
--method copy $META_JSON
rm -r $TMP/$PROJ


################################
## 4) Wang et al. 2022 HMGOs TCF4 ChIP-seq
PROJ=Wang2022_ChIPseq
mkdir -p $TMP/$PROJ && cd $TMP/$PROJ
rsync -Pav $PROJDIR/$PROJ/* $TMP/$PROJ

INPUT_JSON=/projects/pfenninggroup/mouseCxStr/PTHS_interneurons/TCF4_ChIP-seq/${PROJ}/${PROJ}.json
caper run /home/bnphan/src/chip-seq-pipeline2/chip.wdl \
-i "${INPUT_JSON}" --singularity /home/bnphan/src/chip-seq-pipeline2/chip-seq-pipeline_v2.2.1.sif \
--max-concurrent-tasks 1 --db in-memory

META_JSON=/scratch/bnphan/Wang2022_ChIPseq/chip/d91d42aa-b60a-421b-853a-f2e05780005a/metadata.json

## copy over the pipeline outputs
croo --tmp-dir $TMP --out-dir $PROJDIR/$PROJ \
--method copy $META_JSON
rm -r $TMP/$PROJ

