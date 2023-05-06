 


###### for one sample
WorkDir="/NAS/home/chengc/Prj_SCLC/count"
cd $WorkDir
PatientNumber="P1N"
FastqData="/NAS/home/chengc/Prj_SCLC/fq/P1N/"
GRCh38="/NAS/home/chengc/reference/GRCh38/GRCh38/"

cellranger count --id=${PatientNumber} \
--fastqs=${FastqData} \
--transcriptome=${GRCh38} \
--localcores=24 --localmem=100









