prefix=SRR5876943
prefix2=SRR5876944
outpre=Pvag

java -jar ~/opt/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
  -threads 60 \
  -phred33 \
  -summary SummaryFile \
  ./${prefix}_1.fastq.gz \
  ./${prefix}_2.fastq.gz \
  ./${prefix}_1.paired.fq.gz  \
  ./${prefix}_1.unpaired.fq.gz  \
  ./${prefix}_2.paired.fq.gz  \
  ./${prefix}_2.unpaired.fq.gz \
  ILLUMINACLIP:/mnt/Data_disk/zhoubiaofeng/opt/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:20 \
  TRAILING:20 \
  SLIDINGWINDOW:10:20 \
  MINLEN:36


if [ -d trinity_${outpre} ];then /bin/rm -rf trinity_${outpre};fi 
Trinity --seqType fq \
  --max_memory 150G \
  --left ../${prefix}/${prefix}_1.paired.fq.gz,../${prefix2}/${prefix2}_1.paired.fq.gz \
  --right ../${prefix}/${prefix}_2.paired.fq.gz,../${prefix2}/${prefix2}_2.paired.fq.gz \
  --CPU 70 \
  --output trinity_${outpre} \
  --min_kmer_cov 2

~/opt/cdhit/cd-hit-v4.6.7-2017-0501/cd-hit-est \
  -i ./trinity_Pvag.Trinity.fasta \
  -o ./Unigene_Pvag.Trinity_CdsHit0.95.fasta \
  -T 70 \
  -M 0 \
  -c 0.95

 
  if [ -d TransDecoder_LongOrfs ];then /bin/rm -rf TransDecoder_LongOrfs/ TransDecoder_LongOrfs.__checkpoints/ TransDecoder_LongOrfs.__checkpoints_longorfs/; fi 
~/opt/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs \
  -t ./Unigene_Pvag.Trinity_CdsHit0.95.fasta \
  -m 70 \
  --output_dir TransDecoder_LongOrfs
  rm *cmds

if [ -f diamond_blastp.outfmt6 ];then rm diamond_blastp.outfmt6;fi
diamond blastp \
  -d /mnt/Data_disk/zhoubiaofeng/database/Uniprot/uniprot_sprot.fasta \
  -q TransDecoder_LongOrfs/longest_orfs.pep \
  --evalue 1e-5 \
  --max-target-seqs 1 > diamond_blastp.outfmt6

~/opt/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict \
    -t  Unigene_Pvag.Trinity_CdsHit0.95.fasta \
    --retain_blastp_hits diamond_blastp.outfmt6 \
    -T 30 \
    -O TransDecoder_LongOrfs
rm *cmds