#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# RNA-seq quantification workflow
#
# fastp
#   -> TopHat2
#   -> Cuffquant
#   -> Cuffnorm for Root and Shoot
#   -> Cuffdiff for Root samples 02, 03 and 04
###############################################################################

###############################################################################
# 1. Variables
###############################################################################

THREADS=16
FASTP_THREADS=3
DRY_RUN=true

RAWFASTQ_DIR="../readyfastq"

REFERENCE="BioResource/Reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="BioResource/GeneModel/Arabidopsis_thaliana.TAIR10.44.gtf"
MASK_GTF="BioResource/GeneModel/Arabidopsis_thaliana.TAIR10.44.mask.gtf"

# 모든 분석 단계에 사용하는 Root 샘플
ROOT_CONTROL_SAMPLES=(
    Root-Control-01
    Root-Control-02
    Root-Control-03
    Root-Control-04
)

ROOT_PAO1_SAMPLES=(
    Root-PA01-01
    Root-PA01-02
    Root-PA01-03
    Root-PA01-04
)

# 모든 분석 단계에 사용하는 Shoot 샘플
SHOOT_CONTROL_SAMPLES=(
    Shoot-Control-01
    Shoot-Control-02
    Shoot-Control-03
    Shoot-Control-04
)

SHOOT_PAO1_SAMPLES=(
    Shoot-PA01-01
    Shoot-PA01-02
    Shoot-PA01-03
    Shoot-PA01-04
)

ALL_SAMPLES=(
    "${ROOT_CONTROL_SAMPLES[@]}"
    "${ROOT_PAO1_SAMPLES[@]}"
    "${SHOOT_CONTROL_SAMPLES[@]}"
    "${SHOOT_PAO1_SAMPLES[@]}"
)

# Root DEG 분석에만 사용하는 샘플
ROOT_DEG_CONTROL_SAMPLES=(
    Root-Control-02
    Root-Control-03
    Root-Control-04
)

ROOT_DEG_PAO1_SAMPLES=(
    Root-PA01-02
    Root-PA01-03
    Root-PA01-04
)

###############################################################################
# 2. Create output directories
###############################################################################

CMD=(
    mkdir -p \
        00.Rawdata \
        01.Clean \
        02.Align \
        03.Cuffquant \
        04.Cuffnorm/Root \
        04.Cuffnorm/Shoot \
        05.Cuffdiff/Root \
        logs
)
printf '[CMD]'
printf '%q ' "${CMD[@]}"
printf '\n'

if [[ "${DRY_RUN}" == false ]]
then
    "${CMD[@]}"
fi

###############################################################################
# 3. Check and link FASTQ files
###############################################################################

for SAMPLE in "${ALL_SAMPLES[@]}"
do
    R1="${RAWFASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${RAWFASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]
    then
        echo "[ERROR] FASTQ pair not found: ${SAMPLE}" >&2
        exit 1
    fi

    CMD=(
        ln -sfn "$(realpath "${R1}")" \
            "00.Rawdata/${SAMPLE}_R1.fastq.gz"
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '\n'

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}"
    fi

    CMD=(
        ln -sfn "$(realpath "${R2}")" \
            "00.Rawdata/${SAMPLE}_R2.fastq.gz"
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '\n'

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}"
    fi
done

###############################################################################
# fastp
###############################################################################

for SAMPLE in "${ALL_SAMPLES[@]}"
do
    CLEAN_R1="01.Clean/${SAMPLE}_R1.fastq.gz"
    CLEAN_R2="01.Clean/${SAMPLE}_R2.fastq.gz"
    FASTP_JSON="01.Clean/${SAMPLE}.fastp.json"
    FASTP_HTML="01.Clean/${SAMPLE}.fastp.html"

    if [[ -s "${CLEAN_R1}" &&
          -s "${CLEAN_R2}" &&
          -s "${FASTP_JSON}" &&
          -s "${FASTP_HTML}" ]]
    then
        echo "[PASS][fastp] ${SAMPLE}"
        continue
    fi

    echo "[RUN ][fastp] ${SAMPLE}"

    CMD=(
        fastp \
            --thread "${FASTP_THREADS}" \
            --in1 "00.Rawdata/${SAMPLE}_R1.fastq.gz" \
            --in2 "00.Rawdata/${SAMPLE}_R2.fastq.gz" \
            --out1 "${CLEAN_R1}" \
            --out2 "${CLEAN_R2}" \
            --json "${FASTP_JSON}" \
            --html "${FASTP_HTML}" \
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '> %q 2>&1\n' "logs/${SAMPLE}.fastp.log"

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}" > "logs/${SAMPLE}.fastp.log" 2>&1
    fi
done

###############################################################################
# TopHat2
###############################################################################

for SAMPLE in "${ALL_SAMPLES[@]}"
do
    ALIGN_DIR="02.Align/${SAMPLE}"
    BAM="${ALIGN_DIR}/accepted_hits.bam"

    if [[ -s "${BAM}" ]]
    then
        echo "[PASS][TopHat2] ${SAMPLE}"
        continue
    fi

    echo "[RUN ][TopHat2] ${SAMPLE}"

    if [[ "${SAMPLE}" == *"-Control-"* ]]
    then
        GROUP="Control"
    else
        GROUP="PA01"
    fi

    if [[ "${DRY_RUN}" == false && -d "${ALIGN_DIR}" ]]
    then
        echo "[INFO][TopHat2] Removing incomplete output: ${ALIGN_DIR}"
        rm -rf "${ALIGN_DIR}"
    fi

    CMD=(
        tophat \
            -o "${ALIGN_DIR}" \
            -p "${THREADS}" \
            -r 250 \
            --mate-std-dev 50 \
            --library-type fr-unstranded \
            -G "${GTF}" \
            --rg-id "${GROUP}" \
            --rg-sample "${SAMPLE}" \
            "${REFERENCE}" \
            "01.Clean/${SAMPLE}_R1.fastq.gz" \
            "01.Clean/${SAMPLE}_R2.fastq.gz" \
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '> %q 2>&1\n' "logs/${SAMPLE}.tophat.log"

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}" > "logs/${SAMPLE}.tophat.log" 2>&1
        echo "[DONE][TopHat2] ${SAMPLE}"
    fi
done


###############################################################################
# Cuffquant
###############################################################################

for SAMPLE in "${ALL_SAMPLES[@]}"
do
    CUFFQUANT_DIR="03.Cuffquant/${SAMPLE}"
    CXB="${CUFFQUANT_DIR}/abundances.cxb"

    if [[ -s "${CXB}" ]]
    then
        echo "[PASS][Cuffquant] ${SAMPLE}"
        continue
    fi

    echo "[RUN ][Cuffquant] ${SAMPLE}"

    if [[ "${DRY_RUN}" == false && -d "${CUFFQUANT_DIR}" ]]
    then
        echo "[INFO][Cuffquant] Removing incomplete output: ${CUFFQUANT_DIR}"
        rm -rf "${CUFFQUANT_DIR}"
    fi

    CMD=(
        cuffquant \
            --output-dir "${CUFFQUANT_DIR}" \
            --frag-bias-correct "${REFERENCE}" \
            --multi-read-correct \
            --num-threads "${THREADS}" \
            --library-type fr-unstranded \
            --mask-file "${MASK_GTF}" \
            "${GTF}" \
            "02.Align/${SAMPLE}/accepted_hits.bam"
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '> %q 2>&1\n' "logs/${SAMPLE}.cuffquant.log"

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}" > "logs/${SAMPLE}.cuffquant.log" 2>&1
        echo "[DONE][Cuffquant] ${SAMPLE}"
    fi
done


###############################################################################
# Function for creating comma-separated CXB input
###############################################################################

make_cxb_input() {
    local SAMPLE
    local CXB_FILES=()

    for SAMPLE in "$@"
    do
        CXB_FILES+=("03.Cuffquant/${SAMPLE}/abundances.cxb")
    done

    local IFS=,
    echo "${CXB_FILES[*]}"
}


###############################################################################
# Cuffnorm: Root expression table
###############################################################################

ROOT_CONTROL_CXB=$(
    make_cxb_input "${ROOT_CONTROL_SAMPLES[@]}"
)

ROOT_PAO1_CXB=$(
    make_cxb_input "${ROOT_PAO1_SAMPLES[@]}"
)

ROOT_CUFFNORM_DIR="04.Cuffnorm/Root"
ROOT_CUFFNORM_RESULT="${ROOT_CUFFNORM_DIR}/genes.fpkm_table"

if [[ -s "${ROOT_CUFFNORM_RESULT}" ]]
then
    echo "[PASS][Cuffnorm] Root"
else
    echo "[RUN ][Cuffnorm] Root"

    if [[ "${DRY_RUN}" == false && -d "${ROOT_CUFFNORM_DIR}" ]]
    then
        rm -rf "${ROOT_CUFFNORM_DIR}"
    fi

    if [[ "${DRY_RUN}" == false ]]
    then
        mkdir -p "${ROOT_CUFFNORM_DIR}"
    fi

    CMD=(
        cuffnorm \
            --output-dir "${ROOT_CUFFNORM_DIR}" \
            --labels Root-Control,Root-PA01 \
            --num-threads "${THREADS}" \
            --library-type fr-unstranded \
            --library-norm-method classic-fpkm \
            --output-format simple-table \
            "${GTF}" \
            "${ROOT_CONTROL_CXB}" \
            "${ROOT_PAO1_CXB}"
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '> %q 2>&1\n' "logs/Root.cuffnorm.log"

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}" > logs/Root.cuffnorm.log 2>&1
        echo "[DONE][Cuffnorm] Root"
    fi
fi


###############################################################################
# Cuffnorm: Shoot expression table
###############################################################################

SHOOT_CONTROL_CXB=$(
    make_cxb_input "${SHOOT_CONTROL_SAMPLES[@]}"
)

SHOOT_PAO1_CXB=$(
    make_cxb_input "${SHOOT_PAO1_SAMPLES[@]}"
)

SHOOT_CUFFNORM_DIR="04.Cuffnorm/Shoot"
SHOOT_CUFFNORM_RESULT="${SHOOT_CUFFNORM_DIR}/genes.fpkm_table"

if [[ -s "${SHOOT_CUFFNORM_RESULT}" ]]
then
    echo "[PASS][Cuffnorm] Shoot"
else
    echo "[RUN ][Cuffnorm] Shoot"

    if [[ "${DRY_RUN}" == false && -d "${SHOOT_CUFFNORM_DIR}" ]]
    then
        rm -rf "${SHOOT_CUFFNORM_DIR}"
    fi

    if [[ "${DRY_RUN}" == false ]]
    then
        mkdir -p "${SHOOT_CUFFNORM_DIR}"
    fi

    CMD=(
        cuffnorm \
            --output-dir "${SHOOT_CUFFNORM_DIR}" \
            --labels Shoot-Control,Shoot-PA01 \
            --num-threads "${THREADS}" \
            --library-type fr-unstranded \
            --library-norm-method classic-fpkm \
            --output-format simple-table \
            "${GTF}" \
            "${SHOOT_CONTROL_CXB}" \
            "${SHOOT_PAO1_CXB}"
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '> %q 2>&1\n' "logs/Shoot.cuffnorm.log"

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}" > logs/Shoot.cuffnorm.log 2>&1
        echo "[DONE][Cuffnorm] Shoot"
    fi
fi

###############################################################################
# Cuffdiff: Root DEG analysis
###############################################################################

ROOT_DEG_CONTROL_CXB=$(
    make_cxb_input "${ROOT_DEG_CONTROL_SAMPLES[@]}"
)

ROOT_DEG_PAO1_CXB=$(
    make_cxb_input "${ROOT_DEG_PAO1_SAMPLES[@]}"
)

ROOT_CUFFDIFF_DIR="05.Cuffdiff/Root"
ROOT_CUFFDIFF_RESULT="${ROOT_CUFFDIFF_DIR}/gene_exp.diff"

if [[ -s "${ROOT_CUFFDIFF_RESULT}" ]]
then
    echo "[PASS][Cuffdiff] Root DEG"
else
    echo "[RUN ][Cuffdiff] Root DEG"

    if [[ "${DRY_RUN}" == false && -d "${ROOT_CUFFDIFF_DIR}" ]]
    then
        rm -rf "${ROOT_CUFFDIFF_DIR}"
    fi

    if [[ "${DRY_RUN}" == false ]]
    then
        mkdir -p "${ROOT_CUFFDIFF_DIR}"
    fi

    CMD=(
        cuffdiff \
            --output-dir "${ROOT_CUFFDIFF_DIR}" \
            --labels Root-Control,Root-PA01 \
            --frag-bias-correct "${REFERENCE}" \
            --multi-read-correct \
            --num-threads "${THREADS}" \
            --library-type fr-unstranded \
            --library-norm-method classic-fpkm \
            --dispersion-method pooled \
            "${GTF}" \
            "${ROOT_DEG_CONTROL_CXB}" \
            "${ROOT_DEG_PAO1_CXB}"
    )
    printf '[CMD]'
    printf '%q ' "${CMD[@]}"
    printf '> %q 2>&1\n' "logs/Root.cuffdiff.log"

    if [[ "${DRY_RUN}" == false ]]
    then
        "${CMD[@]}" > logs/Root.cuffdiff.log 2>&1
        echo "[DONE][Cuffdiff] Root DEG"
    fi
fi
