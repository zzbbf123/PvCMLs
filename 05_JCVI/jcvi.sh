type1=mRNA
key1=transcript_id
key2=ID

####### jcvi
while read line;do 
  spname=`echo ${line} | cut -d '=' -f 1`
  gffpath=`echo ${line} | cut -d '=' -f 2`
  if [ ${spname} == "!Osat" ];then
     $python -m jcvi.formats.gff bed --type=${type1} --key=${key2} ${gffpath} > ${spname}.bed
     $python -m jcvi.formats.bed uniq ${spname}.bed
  fi 

  if [ ${spname} == '!Atha' ];then 
     $python -m jcvi.formats.gff bed --type=${type1} --key=${key2} ${gffpath} > ${spname}.bed
     $python -m jcvi.formats.bed uniq ${spname}.bed
  fi 
  if [ ${spname} == 'Pvag' ];then
     $python -m jcvi.formats.gff bed --type=${type1} --key=${key2} ${gffpath} > ${spname}.bed
     $python -m jcvi.formats.bed uniq ${spname}.bed
  fi
done < gff.list

if [ -f command_line ];then rm command_line;fi 

while read line;do
  spname=`echo ${line} | cut -d '=' -f 1`
  peppath=`echo ${line} | cut -d '=' -f 2`
  if [ -f ${spname}.pep ];then rm ${spname}.pep ;fi 
  echo "seqkit grep -f <(cut -f 4 ${spname}.uniq.bed) ${peppath} | seqkit seq -i > ${spname}.pep" 
done < pep.list >> command_line

while read line;do
  spname=`echo ${line} | cut -d '=' -f 1`
  cdspath=`echo ${line} | cut -d '=' -f 2`
  if [ -f ${spname}.cds ];then rm ${spname}.cds ;fi 
  echo "seqkit grep -f <(cut -f 4 ${spname}.uniq.bed) ${cdspath} | seqkit seq -i > ${spname}.cds" 
done < cds.list >> command_line

bash command_line 

mkdir -p cds && cd cds
while read line;do
  spname=`echo ${line} | cut -d '=' -f 1`
  ln -sf ../${spname}.cds ${spname}.cds
  ln -sf ../${spname}.bed ${spname}.bed
done < ../gff.list

cd cds && $python -m jcvi.compara.catalog ortholog --no_strip_names --notex --cpus=40  Osat Atha