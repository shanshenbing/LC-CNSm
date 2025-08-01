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




#settings
MAIN_PATH=
WORK_DIR=${MAIN_PATH}/celescope

#new data
mymappfile=${MAIN_PATH}/script/celescope/mapSampleided.txt

genomeDir=/dir_to_path/software/CeleScope/data/hs.gencode.v32_GRCh38


conda activate celescope
#prepare dir
dir_mk ${WORK_DIR}/result
cd ${WORK_DIR}/result
#Generate scripts for each sample
multi_rna \
	--mapfile ${mymappfile} \
	--genomeDir ${genomeDir}\
	--thread 8\
	--mod shell