#!/bin/bash -
#===============================================================================
#
#          FILE: cellranger.expr.sh
#
#         USAGE: ./cellranger.expr.sh
#
#   DESCRIPTION: 
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Shenbing shan (), 1186426740@qq.com
#  ORGANIZATION: PICB
#       CREATED: 2021年06月18日 16时35分51秒
#      REVISION:  ---
#===============================================================================

#set -o nounset                                  # Treat unset variables as an error



source ~/.bashrc


sampleid=$1

#settings
MAIN_PATH=/path_to_dir/projects/LCBMNew
WORK_DIR=${MAIN_PATH}

expr_fq_path=${WORK_DIR}/fq/${sampleid}
#bcr_fq_path=${WORK_DIR}/fq/${sampleid}/bcr
#tcr_fq_path=${WORK_DIR}/fq/${sampleid}/tcr


expr_ref=/path_to_dir/PublicData/reference/cellranger/refdata-gex-GRCh38-2020-A
#vdj_ref=/path_to_dir/PublicData/reference/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

#prepare dir
dir_mk ${WORK_DIR}/cellranger/${sampleid}/expr/run_count
#dir_mk ${WORK_DIR}/cellranger/${sampleid}/bcr/run_count
#dir_mk ${WORK_DIR}/cellranger/${sampleid}/tcr/run_count
#move to out dir
cd ${WORK_DIR}/cellranger/${sampleid}/expr
#expr
cellranger count --id=run_count --fastqs=${expr_fq_path} --sample=${sampleid} --transcriptome=${expr_ref}

#for sample with BCR and TCR info only.
#BCR
#cd ${WORK_DIR}/cellranger/${sampleid}/bcr
#cellranger vdj --chain=IG --id=run_count --fastqs=${bcr_fq_path} --sample=${sampleid}-10XBCR --reference=${vdj_ref}
#TCR
#cd ${WORK_DIR}/cellranger/${sampleid}/tcr
#cellranger vdj --chain=TR --id=run_count --fastqs=${tcr_fq_path} --sample=${sampleid}-10XTCR --reference=${vdj_ref}
