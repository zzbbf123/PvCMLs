spname=Pvag
genome=/mnt/Data_disk/zhoubiaofeng/database/genome/genome_database/Paspalum_vaginatum/Phytozome/PhytozomeV13/Pvaginatum/v3.1/Pvaginatum_672_v3.0.fa
GFF="/mnt/Data_disk/zhoubiaofeng/database/genome/genome_database/Paspalum_vaginatum/Phytozome/PhytozomeV13/Pvaginatum/v3.1/Pvaginatum_672_v3.1.gene_exons.gff3"
readsdir="/mnt/Data_disk/zhoubiaofeng/database/haibin/2016PlantCell_Cold/"

## genome build
if [ ! -d reference ];then mkdir reference;fi
ln -sf ${genome} reference/${spname}.fa
if [ ! -f reference/${spname}.6.ht2 ];then
hisat2-build reference/${spname}.fa reference/${spname}
fi 

## single_seq_build
if [ ! -d raw_data ];then mkdir raw_data;fi
for i in SRR4280403 SRR4280409 SRR4280414 SRR4280420 SRR4280426 SRR4280427;do 
   ln -sf ${readsdir}/${i}/${i}.fastq.gz raw_data/${i}.fq.gz
done

for i in SRR4280403 SRR4280409 SRR4280414 SRR4280420 SRR4280426 SRR4280427;do
  fastqc raw_data/${i}.fq.gz & 
done

for f in SRR4280403 SRR4280409 SRR4280414 SRR4280420 SRR4280426 SRR4280427;do
java -jar /mnt/Data_disk/zhoubiaofeng/opt/Trimmomatic-0.38/trimmomatic-0.38.jar \
  SE -threads 1 -phred33 \
  -summary raw_data/${i}.log \
  raw_data/${f}.fq.gz \
  raw_data/${f}.clean.fq.gz \
  ILLUMINACLIP:/mnt/Data_disk/zhoubiaofeng/opt/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:5 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 &
done

if [ ! -d hisat2Result ];then mkdir hisat2Result;fi
  f=$1
  #hisat2 -x /mnt/Data_disk/zhoubiaofeng/running/YMZ/pipeline_for_gf_id/Pvag/phytozome_CML/transcription/refRNAseq/2020BMC_salt/reference/${spname} -U raw_data/${f}.clean.fq.gz -S hisat2Result/${f}.sam -p 5 & 
  #samtools view -bS hisat2Result/${f}.sam -o hisat2Result/${f}.bam 
  #samtools sort -@ 8 hisat2Result/${f}.bam -o hisat2Result/${f}.sorted.bam
  export PATH=${PATH}:/mnt/Data_disk/zhoubiaofeng/opt/stringtie-2.1.6.Linux_x86_64

## Transcriptome assembly
export PATH=${PATH}:/mnt/Data_disk/zhoubiaofeng/opt/stringtie-2.1.6.Linux_x86_64
#for f in SRR4280403 SRR4280409 SRR4280414 SRR4280420 SRR4280426 SRR4280427;do  stringtie -p 20 -G ${GFF} -o hisat2Result/${f}.gtf hisat2Result/${f}.sorted.bam ; done
#for f in SRR4280403 SRR4280409 SRR4280414 SRR4280420 SRR4280426 SRR4280427;do echo hisat2Result/${f}.gtf ;done > gtf.txt

# Merges transcripts into a non-redundant set of transcripts
#stringtie --merge -p 20 -G ${GFF} -o hisat2Result/merged.gtf gtf.txt

if [ ! -d ballgown ];then mkdir ballgown ; fi
# Expression level estimation
for f in SRR4280403 SRR4280409 SRR4280414 SRR4280420 SRR4280426 SRR4280427;do 
  stringtie -e -B -p 8 -G hisat2Result/merged.gtf -o ballgown/${f}/${f}_expression.gtf hisat2Result/${f}.sorted.bam &
done

export PATH=${PATH}:/mnt/Data_disk/zhoubiaofeng/opt/stringtie-2.1.6.Linux_x86_64
#prepDE.py3 -i sample.info

for i in `cat CML.id`;do grep ${i} hisat2Result/merged.gtf | awk -F '\t' '$3 == "transcript"{print $9}';done | awk -F' ' '{print $4"\t"$2}' | sed -e 's#\"##g' -e 's#\;##g'| sort -k1,1nr > CML.transcript.id