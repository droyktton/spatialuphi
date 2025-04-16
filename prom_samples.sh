samples=$1

echo "samples="$samples

for((i=1;i<$samples;i++))
do 
    echo $i; 
    ./a.out 1024 1000000 $i; 
    cp inst_sofq.dat "inst_sofq_"$i"_.dat"; 
    cp cm.dat "cm_"$i"_.dat";
done

paste inst_sofq_*_.dat | awk '{acum=0; for(i=0;i<NF;i++){acum+=$i}; if(NF>0) print acum*1.0/NF; else print;}' \
> "sofq_"$samples"samples.dat"

paste cm_*_.dat | awk '{acum=0; for(i=0;i<NF;i++){acum+=$i}; if(NF>0) print acum*1.0/NF; else print;}' \
> "cm_"$samples"samples.dat"

