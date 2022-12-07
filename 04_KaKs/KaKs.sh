### prepare data for kaks calculation
if [ ! -d prepare_d ];then mkdir prepare_d;fi 
n1=Osat
n2=Atha
n3=Pvag
rename=${n1}_${n2}_${n3}
gf_name=CML

#### prepare homo list 

  collinear_f=./prepare_d/${rename}.rep.collinearity
  for gene in `cat ./prepare_d/${gf_name}.geneid`;do grep ${gene} ${collinear_f} | grep -v ^"#" | awk -F '\t' '{print $2,$3}' ;done    | grep -v Atha | sed 's#Osat_##g' > ./prepare_d/${gf_name}.${n1}_${n3}.homo #only one is fine
  
  for gene in `cat ./prepare_d/${gf_name}.geneid`;do grep ${gene} ${collinear_f} | grep -v ^"#" | awk -F '\t' '{print $2,$3}' ;done    | grep -v Osat | sed 's#Atha_##g'> ./prepare_d/${gf_name}.${n2}_${n3}.homo #only one is fine
  
  for gene in `cat ./prepare_d/${gf_name}.geneid`;do grep ${gene} ${collinear_f} | grep -v ^"#" | awk -F '\t' '{print $2,$3}' ;done    | grep -v Pavag | sed 's#Atha_##g' | sed 's#Osat_##g' > ./prepare_d/${gf_name}.${n1}_${n2}.homo #only one is fine


  if [ -f ./prepare_d/${gf_name}.homo.checked.list ] ;then rm ./prepare_d/${gf_name}.homo.checked.list;fi 
  while read line;do 
    homo1=`echo ${line} | cut -d ' ' -f 1 `
    homo2=`echo ${line} | cut -d ' ' -f 2 `
    ishomo=0
    if [[ `grep $homo1 ./prepare_d/${gf_name}.geneid` -lt 1 ]] || [[ `grep ${homo2} ./prepare_d/${gf_name}.geneid` -lt 1 ]];then 
       echo ${line} >> ./prepare_d/${gf_name}.checked.homo
    else 
       echo -e "${line} is not CML homology ! one of them donot belong to CML family ."
    fi 
  done  < ./prepare_d/${gf_name}.homo


#### prepare cds fasta   ; need manually check
  cds_f3=/mnt/Data_disk/zhoubiaofeng/database/genome/download_genome_seq_forCML/Paspalum_vaginatum/Phytozome/PhytozomeV13/Pvaginatum/v3.1/Pvaginatum_672_v3.1.cds.fa
  cds_f1=/mnt/Data_disk/zhoubiaofeng/database/genome/Orazy_sative/release7/all.cds
  cds_f2=/mnt/Data_disk/zhoubiaofeng/database/genome/Orazy_sative/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.cds.reduc.fa
  gff_f2=/mnt/Data_disk/zhoubiaofeng/database/genome/Orazy_sative/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff
  for g in `grep ${n1} ./prepare_d/CML.geneid  |sed 's#'${n1}'_##g'`;do grep ${g} redundenceData/${n1}.gene.locus.pep.len.sort.uniq | awk -F '\t' '{print $3}';done > ./prepare_d/${n1}.geneid     # same as protein id
  grep ${n2} ./prepare_d/CML.geneid  |sed 's#'${n2}'_##g' > ./prepare_d/${n2}.geneid # have done some changes
  for g in `grep Pavag ./prepare_d/CML.geneid`;do grep ${g} redundenceData/${n3}.gene.locus.pep.len.sort.uniq | awk -F '\t' '{print $3}';done  > ./prepare_d/${n3}.geneid  # same as protein id
  python get_geneCDS.py ${cds_f1} ./prepare_d/${n1}.geneid
  python get_geneCDS.py ${cds_f2} ./prepare_d/${n2}.geneid
  python get_geneCDS.py ${cds_f3} ./prepare_d/${n3}.geneid


#### prepare pep fasta  manually check 
  pep_f1=../../phytozome_TIFY/redundenceData/${n1}.longest.uniq.fasta
  pep_f2=../../phytozome_TIFY/redundenceData/${n2}.longest.uniq.fasta
  pep_f3=../../phytozome_TIFY/redundenceData/${n3}.longest.uniq.fasta
  for g in `grep ${n1} ./prepare_d/CML.geneid | sed 's#'${n1}'_##g'`;do grep ${g} redundenceData/${n1}.gene.locus.pep.len.sort.uniq | awk -F '\t' '{print $2,$3}';done > ./prepare_d/${n1}.pep.geneid
  for g in `grep ${n2} ./prepare_d/CML.geneid | sed 's#'${n2}'_##g'`;do grep ${g} redundenceData/${n2}.gene.locus.pep.len.sort.uniq | awk -F '\t' '{print $2,$3}';done > ./prepare_d/${n2}.pep.geneid
  for g in `grep Pavag ./prepare_d/CML.geneid `;do grep ${g} redundenceData/${n3}.gene.locus.pep.len.sort.uniq | awk -F '\t' '{print $2,$3}';done > ./prepare_d/${n3}.pep.geneid
  python get_genePEP.py ${pep_f1} ./prepare_d/${n1}.pep.geneid
  python get_genePEP.py ${pep_f2} ./prepare_d/${n2}.pep.geneid
  python get_genePEP.py ${pep_f3} ./prepare_d/${n3}.pep.geneid



###### run kaks in paml
n1=Osat
n2=Atha
n3=Pvag
rename=${n1}_${n2}_${n3}
gf_name=CML


sp1=$1
sp2=$2
cat ./prepare_d/${n1}.geneid.fa ./prepare_d/${n2}.geneid.fa ./prepare_d/${n3}.geneid.fa > ./prepare_d/CML.cds.fa
cat ./prepare_d/${n1}.pep.geneid.fa ./prepare_d/${n2}.pep.geneid.fa ./prepare_d/${n3}.pep.geneid.fa > ./prepare_d/CML.pep.fa

#### Alignment
if [ -d CML_${sp1}_${sp2} ];then rm -rf CML_${sp1}_${sp2};fi
~/opt/ParaAT2.0/ParaAT.pl \
  -h ./prepare_d/CML.${sp1}_${sp2}.homo \
  -n ./prepare_d/CML.cds.fa \
  -a ./prepare_d/CML.pep.fa \
  -p proc \
  -o CML_${sp1}_${sp2} \
  -f paml \
  -m mafft

#### run paml with yn00 model
  sp1=$1
  sp2=$2
  cp paml-yn00-run-pipeline-master/yn00_pre.ctl CML_${sp1}_${sp2}/yn00.ctl2
  cp paml-yn00-run-pipeline-master/run.yn00_v2.sh CML_${sp1}_${sp2}/
  cd CML_${sp1}_${sp2}
  sh run.yn00_v2.sh
  cd ..
