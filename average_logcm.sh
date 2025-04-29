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
            counter[i]=0;acum[i]=0;
        } 
        for(j=1;j<=NF;j++)
        {
            counter[(j-1)%9]++; 
            acum[(j-1)%9]+=$j;
        }; 
        for(i=0;i<9;i++) 
        printf("%f ",acum[i]/counter[i]); 
        printf("%d\n", counter[8]); 
    }    
}'

#plot '< bash average_logcm.sh' u 1:6 w lp, "" u 1:7 w lp, for[z in "1.25 1.5 1.75"] 0.1*x**(2*0.5/z) t 'z='.z