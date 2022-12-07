#### Blast search
if [ ! -d prepare_d ];then mkdir prepare_d ;fi  # for store the prepare data
if [ ! -d db ];then mkdir db ;fi  # for store the blast db

# prepare GFF3
while read line;do 
  spname=Pvag
  gff=`echo ${line}| cut -d '=' -f 2 `
  uniq=../redundenceData/${spname}.gene.locus.pep.len.sort.uniq

  # get gene region
  awk -F '[\t;]' '$3=="gene"{for(i=1;i<=NF;i++){if($i ~/^ID=/) printf $1"\t"$4"\t"$5"\t"$i"\n"}}' ${gff} | sed 's#ID=##g'| sed 's#\.v3\.1##g' > pre.gff

  # add length info
  awk -F '\t' 'NR == FNR{d[$1]=$3;next} NR!=FNR && d[$4] != "" {printf "%s\t%s\t%s\t%s\n",$1,d[$4],$2,$3}' ${uniq} pre.gff > prepare_d/${spname}.gff
done  < query_annotation

# prepare fasta
# only for research species
 spname=$1
 if [ ! -f db/${spname}.pep.fasta ];then cd db; ln -s ../../redundenceData/${spname}.longest.uniq.fasta ${spname}.pep.fasta; cd ..;fi
 if [ ! -f db/${spname}.pep.fasta.ndb ];then
 makeblastdb -parse_seqids -hash_index -dbtype prot -in db/${spname}.pep.fasta -out db/${spname}.pep.fasta
 fi
 blastp -query ../redundenceData/${spname}.longest.uniq.fasta -db db/${spname}.pep.fasta -outfmt 6 -out db/${spname}.blast.txt
 -num_threads 10 -evalue 1e-5 & # self blast

spname1=$1   # db
spname2=$2   # query
if [ ! -f db/${spname1}.pep.fasta ];then cd db; ln -s ../../redundenceData/${spname1}.longest.uniq.fasta ${spname1}.pep.fasta;cd ..;fi
if [ ! -f db/${spname1}.pep.fasta.ndb ];then
    makeblastdb -parse_seqids -hash_index -dbtype prot -in db/${spname1}.pep.fasta -out db/${spname1}.pep.fasta
fi
blastp -query ../redundenceData/${spname2}.longest.uniq.fasta -db db/${spname1}.pep.fasta -outfmt 6 -out db/${spname2}_${spname1}.blast.txt -num_threads 10 -evalue 1e-5   # among blast

# Note: prepare .rep file
  spname=$1
  gff=./prepare_d/${spname}.gff
  if [ -f prepare_d/${spname}.rep.gff ];then rm prepare_d/${spname}.rep.gff ;fi
  awk -F '\t' '{gsub(/\.[0-9]+/, "", $2); print}' ${gff} | sed 's# #\t#g' > ./prepare_d/${spname}.rep.gff

for f in `ls db/*blast.txt `;do echo 
  if [ -f prepare_d/${f}.rep.blast ];then rm prepare_d/${f}.rep.blast;fi
  blast=db/${f}.blast.txt
  awk -F '\t' '{gsub(/\.[0-9]+/, "", $1); gsub(/\.[0-9]+/, "", $2); print}' ${blast} | sed 's# #\t#g' > ./prepare_d/${f}.blast
  python only_best.py `echo ${i##*/} | sed 's#.blast.txt##g'` &  
done

## run MCScanX
if [ ! -d prepare_d ];then echo "do not detected dir prepare_d, please check your file";fi
MS="/mnt/Data_disk/zhoubiaofeng/opt/MCScanX"
MS_downstream="/mnt/Data_disk/zhoubiaofeng/opt/MCScanX/downstream_analyses"

# generate comb / can be 1 or 2 or 3 ..
n1=Osat
n2=Atha
n3=Pvag
f=${n1}_${n2}_${n3}
# manually change the .rep.blast file name
cat prepare_d/${n2}_${n1}.rep.blast prepare_d/${n3}_${n2}.rep.blast prepare_d/${n3}_${n1}.rep.blast > prepare_d/${f}.rep.blast
cat prepare_d/${n1}.rep.gff prepare_d/${n2}.rep.gff prepare_d/${n3}.rep.gff > prepare_d/${f}.rep.gff

${MS}/MCScanX prepare_d/${f}.rep
${MS_downstream}/detect_collinear_tandem_arrays \
  -g prepare_d/${f}.rep.gff \
  -b prepare_d/${f}.rep.blast \
  -c prepare_d/${f}.rep.collinearity \
  -o prepare_d/${f}.rep.detect_collinear_tandem_arrays

## run Circos
## get subdata from genome-wide collinear : only for CaCA pair; 
rename=Osat_Atha_Pvag
gf_name=CML
array=()   # this step is for save the combs
length=`echo ${rename} | awk -F '_' '{print NF}'`
for((i=1;i<=$length;i++));do array[${i}]=`echo ${rename} | cut -d '_' -f ${i}`;done # add in
echo "species used: ${array[@]}"

## run circos only plot CML genes
# circos config file store at : /pkgs/circos-0.69.8-hdfd78af_1/etc/
cp circos_my.conf circos_${gf_name}.conf; sed -i 's#karyotype_file_name#karyotype.'${rename}'.txt#g' circos_${gf_name}.conf
sed -i 's#link_file_name#link_'${gf_name}'#g' circos_${gf_name}.conf
cp link.conf link_${gf_name}.conf; sed -i 's#collinearity_block_file#'${gf_name}'.collinearity.blocks#g' link_${gf_name}.conf
circos -conf circos_${gf_name}.conf