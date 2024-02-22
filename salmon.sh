#Salmon
for fn in ./Bmal/trimgalore/*S{1..16}_trim;
do
samp=`basename ${fn} _trim` 
echo "Processing sample ${samp}"
salmon quant -i ../../leptin/live/salmon_index/m27index_salmon -l A \
    -1 <(gunzip -c ${fn}/${samp}_R1_001_val_1.fq.gz) \
    -2 <(gunzip -c ${fn}/${samp}_R2_001_val_2.fq.gz) \
    -p 8 --validateMappings --gcBias --seqBias -o ./Bmal/salmon/${samp}_quant
done