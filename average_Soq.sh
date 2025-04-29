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
> "sofq_samples.dat"

#paste cm_*_.dat | awk '{acum=0; for(i=0;i<NF;i++){acum[i%6]+=$i}; if(NF>0){for(i=0;i<6;i++) print acum[i]*6/NF}; else print;}' \
#> "cm_"$samples"samples.dat"

#set multi lay 2,2; do for[z in "1.2 1.3 1.4 1.5"]{set tit "z=".z; plot [:][0.1:] for[i=7:12] 'sofq_21samples.dat' index i u ($0*(t=2**i)**(2*0.5/z)):($2/t**((1+2*0.5)/z)) w l t sprintf('t=%1.2e',2**i);}; unset multi
