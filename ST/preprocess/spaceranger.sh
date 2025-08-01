#!/bin/bash -
###
 # @Author: shanshenbing 1186426740@qq.com
 # @Date: 2025-06-23 16:38:12
 # @LastEditors: shanshenbing 1186426740@qq.com
 # @LastEditTime: 2025-06-23 18:05:34
 # @FilePath: 
 # @Description: 
 # 
 # Copyright (c) 2025 by shanshenbing, All Rights Reserved. 
### 
##SBATCH -p q1,q2,q3
#SBATCH -p low,high,fat
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=10G
#SBATCH --mem=300G


sample_id=$1

source ~/.bashrc
#settings
MAIN_PATH=/path_2_dir/project/LMST

fq=/path_2_dir/data/raw_upload/20250305/st_batch1/raw_data/${sample_id}
#images
sample_tif=/path_2_dir/images/${sample_id}/${sample_id}.tiff
cell_tif=$(ls /path_2_dir/images/${sample_id}/*_${sample_id}.tif)
json_file=$(ls /path_2_dir/images/${sample_id}/*.json)


ref_file_dir=/path_2_dir/PublicData/softwareData/spaceranger/refdata-gex-GRCh38-2020-A/
probe=/path_2_dir/PublicData/softwareData/spaceranger/probe/Human_Transcriptome_v2/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv

filename=$(basename $json_file)
if [[ $filename =~ ^([^-]+-[^-]+)-([^-]+)- ]]; then
    cellid="${BASH_REMATCH[1]}"  
    zoneid="${BASH_REMATCH[2]}"
else
    echo "No match found"
    break
fi

echo "sample id: $sample_id"
echo "fq:  $fq"
echo "image: $sample_tif"
echo "ref: $ref_file_dir"
echo "json: $json_file"
echo "cell_tif: $cell_tif"
echo "slide id: $cellid"
echo "area id: $zoneid"

dir_mk ${MAIN_PATH}/data/ST/spaceranger/${sample_id}

spaceranger count \
	--id=${sample_id} \
	--transcriptome=${ref_file_dir} \
	--fastqs=${fq} \
	--sample=${sample_id} \
	--image=${sample_tif} \
	--slide=${cellid} \
	--area=${zoneid} \
	--localcores=32 \
	--localmem=200 \
	--create-bam true \
	--probe-set=${probe} \
	--cytaimage=${cell_tif} \
	--loupe-alignment=${json_file} \
    --output-dir=${MAIN_PATH}/data/ST/spaceranger/${sample_id}
