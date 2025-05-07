samples=$(ls logcm????.dat | wc -l)

cols=$(awk '{ print NF; exit }' logcm1234.dat)

echo "#samples="$samples
echo "#cols="$cols

paste logcm????.dat | \
gawk -v cols=$cols \
'
function abs(x) { return x < 0 ? -x : x }
BEGIN{
	printf("t velu evelu velphi evelphi cmu ecmu cmphi ecmphi cmu2 ecmu2 cmphi2 ecmphi2 maxu emaxu maxphi emaxphi");
}
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
            acum2[i]=sqrt(abs(acum2[i] - acum[i]*acum[i]))/sqrt(counter[i]);
            printf("%f %f ",acum[i], acum2[i]);         
        } 
        printf("%d\n", counter[8]); 
    }    
}' > "logcm_samples.dat"


#plot 'logcm_samples.dat' u 1:11:12 w lp, "" u 1:13:14 w lp, for[z in "1.25 1.5 1.75 2.0"] 0.1*x**(2*0.5/z) t 'z='.z
