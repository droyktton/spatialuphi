samples=$(ls logcm????.dat | wc -l)

cols=$(awk '{ print NF; exit }' logcm1234.dat)

echo "#samples="$samples
echo "#cols="$cols

paste logcm????.dat | \
awk -v cols=$cols \
'
{
    if(NR>1 && NF>=9)
    {
        for(i=0;i<9;i++)
        {
            counter[i]=0;
            acum[i]=0;
            acum2[i]=0;
        } 
        for(j=1;j<=NF;j++)
        {
            counter[(j-1)%9]++; 
            acum[(j-1)%9]+=$j;
            acum2[(j-1)%9]+=$j*$j;
        }; 
        for(i=0;i<9;i++){
            acum[i]=acum[i]/counter[i];
            acum2[i]=acum2[i]/counter[i];
            acum2[i]=sqrt(acum2[i] - acum[i]*acum[i]);
            printf("%f %f ",acum[i], acum2[i]);     
        } 
        printf("%d\n", counter[8]); 
    }    
}'

#plot '< bash average_logcm.sh' u 1:6 w lp, "" u 1:7 w lp, for[z in "1.25 1.5 1.75"] 0.1*x**(2*0.5/z) t 'z='.z