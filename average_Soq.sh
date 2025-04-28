samples=$(ls inst_sofq????.dat | wc -l)

paste inst_sofq????.dat | \
awk '
{
    acum=0; acum2=0; 
    for(i=0;i<NF;i+=2){
        acum+=$i; acum2+=$(i+1);
    }; 
    if(NF>0) print acum*0.5/NF,acum2*0.5/NF; 
    else print;
}' \
> "sofq_"$samples"samples.dat"

#paste cm_*_.dat | awk '{acum=0; for(i=0;i<NF;i++){acum[i%6]+=$i}; if(NF>0){for(i=0;i<6;i++) print acum[i]*6/NF}; else print;}' \
#> "cm_"$samples"samples.dat"
