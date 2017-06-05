
for analysis in $(find $PWD -name "commandline.sh"); do

    cd $(dirname $analysis)
    sh $analysis

done
