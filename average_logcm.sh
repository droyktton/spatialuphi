samples=$(ls logcm????.dat | wc -l)

cols=$(awk '{ print NF; exit }' logcm1234.dat)

paste logcm????.dat | \
awk -v cols=$cols '
{
    for(i=0;i<cols;i++){
        acum[i]=0;
        counter[i]=0;
    };
    for(i=1;i<=NF;i+=cols){
        acum[(i-1)%cols]+=$i;
        counter[(i-1)%cols]+=1;
    }; 
    if(NF>0){
        for(i=0;i<cols;i++){
            print acum[i]/counter[i];
        }
    }  
    else print;
}' \
> "logcm_"$samples"samples.dat"

#paste cm_*_.dat | awk '{acum=0; for(i=0;i<NF;i++){acum[i%6]+=$i}; if(NF>0){for(i=0;i<6;i++) print acum[i]*6/NF}; else print;}' \
#> "cm_"$samples"samples.dat"
