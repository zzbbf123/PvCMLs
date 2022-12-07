input=$1
output=${input%.*}
echo ${output}
muscle -in ${input} -out ${output}