#!/bin/bash
set -xe

if [[ ! -e ./ERR3240275/ERR3240275_1.fastq.gz ]]; then
 mkdir -p ERR3240275
 cd ERR3240275
 echo "null" > null.txt
 ls | grep -v -E 'ERR3240275_1.fastq.gz' | grep -v -E 'ERR3240275_2.fastq.gz' | xargs rm -r
 cd ..
fi

if [[ ! -e ./SRR8315715/SRR8315715_1.fastq.gz ]]; then
 mkdir -p SRR8315715
 cd SRR8315715
 echo "null" > null.txt
 ls | grep -v -E 'SRR8315715_1.fastq.gz' | grep -v -E 'SRR8315715_2.fastq.gz' | xargs rm -r
 cd ..
fi

echo "null" > null.txt
ls | grep -v -E 'test.sh' | grep -v -E 'ERR3240275' | grep -v -E 'SRR8315715' | xargs rm -r

cwltool  --tmp-outdir-prefix=${PWD}/tmp_cwl/ --tmpdir-prefix=${PWD}/tmp_cwl/ --rm-tmpdir ../workflow/createindex.cwl ../workflow/createindex.job.yaml

cd ERR3240275
if [[ ! -e ./ERR3240275_1.fastq.gz ]]; then
#  fasterq-dump -e 8 ERR3240275 -p
#  pigz ERR3240275_*.fastq
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/005/ERR3240275/ERR3240275_1.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/005/ERR3240275/ERR3240275_2.fastq.gz
fi
cwltool --tmp-outdir-prefix=${PWD}/tmp_cwl/ --tmpdir-prefix=${PWD}/tmp_cwl/ --rm-tmpdir ../../workflow/VIRTUS.PE.cwl ../../workflow/VIRTUS.PE.job.yaml
cd ..


cd SRR8315715
if [[ ! -e ./SRR8315715_1.fastq.gz ]]; then
#  fasterq-dump -e 8 SRR8315715 -p
#  pigz SRR8315715_*.fastq
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR831/005/SRR8315715/SRR8315715_1.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR831/005/SRR8315715/SRR8315715_2.fastq.gz
fi
cwltool --tmp-outdir-prefix=${PWD}/tmp_cwl/ --tmpdir-prefix=${PWD}/tmp_cwl/ --rm-tmpdir ../../workflow/VIRTUS.SE.cwl ../../workflow/VIRTUS.SE.job.yaml
cd ..

echo Successfully completed test.sh!
