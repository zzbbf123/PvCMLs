if [ ! -d db ];then mkdir db fi
if [ -f searchall.sh ];then rm searchall.sh fi

while read line; do
  e_line=`echo ${line} | tr "=" "\n"`
  array=($e_line)
  query=${array[1]}
  out=${array[0]}
  for db in `cat db.list`
  do
    if [ -f ${db}_${out}.out ];then
       rm ${db}_${out}.out
    fi
    if [ ! -f db/${db}.hmm ];then
      hmmbuild --amino db/${db}.hmm db/${db}.txt
    fi
    echo "hmmsearch -E 1e-5 db/${db}.hmm ${query} > ${db}_${out}.out " >> searchall.sh
  done
done < query_list

sh searchall.sh

rm *pos *phase1
while read line 
do
e_line=`echo ${line} | tr "=" "\n"`
array=($e_line)
query=${array[1]}
out=${array[0]}
rm use_list

## integration of different results of Pfam id 
for i in *${out}.out
do
        if [ `grep "No targets" ${i} |wc -l` -eq 1 ];then
                echo "pass"
        else
                echo ${i} >> use_list
        fi
done

sed -i 's#.out##g' use_list

### 
for i in `cat use_list`
do
start=$((`grep -n ^'Score' ${i}.out | cut -d ':' -f1` + 4))
end=$((`grep -n ^'Domain annotation' ${i}.out | cut -d ':' -f1` - 1))
echo $start $end
sed -n "${start},${end}p" ${i}.out > ${i}.phase1
sed -i "s# #\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
sed -i "s#\t\t#\t#g" ${i}.phase1
cut -d $'\t' -f10 ${i}.phase1 |sort -k2 |uniq | grep -v "^$" > ${i}.pos
done

done <  query_list

while read line 
do
e_line=`echo ${line} | tr "=" "\n"`
array=($e_line)
query=${array[1]}
out=${array[0]}

if [ ! -f ${query%.*}"_oneline.fasta" ];then
   python ~/.script/oneline.py ${query}
else
   echo "Oneline file exist"
fi

protein_data=`echo ${query%.*}"_oneline.fasta"`
if [ -f Tify_candidate_${out}.fasta ];then rm Tify_candidate_${out}.fasta ;fi
for n in `grep -nf PF06200_seed_${out}.pos ${protein_data}| cut -d ':' -f 1`
do
sed -n "$n,$(($n + 1))p" ${protein_data} >> Tify_candidate_${out}.fasta
done

done < query_list