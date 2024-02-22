#Trim Galore 

for fn in ./Ilaria_Palmisano_SOUK008651/*S{1..16};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
trim_galore --quality 20 --phred33 --fastqc \
    --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
    --stringency 3 --length 50 \
    -paired ${fn}/*R1_001.fastq.gz ${fn}/*R2_001.fastq.gz \
    -o ./trimgalore/${samp}_trim 
done