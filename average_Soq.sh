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


#set multi lay 2,2; do for[z in "1.3 1.4 1.5 1.6"]{set tit "z=".z; plot [:][:] for[i=14:19] 'sofq_samples.dat' index i u ($0*(t=2**i)**(2*0.5/z)):($2/t**((1+2*0.5)/z)) w l t sprintf('t=%1.2e',2**i), 1e10/x**2}; unset multi


# Grafico de KPZ scaling
#set multi lay 1,2; set logs; set xla "q t^{1/z}"; set yla "S(q,t)/t^{(1+2 zeta)/z}"; do for[z in "1.5"]{set tit "zeta=1/2, z=".z; plot [:1e7][1e-2:] for[i=14:19] 'sofq_samples.dat' index i u ($0*(t=2**i)**(2*0.5/z)):($2/t**((1+2*0.5)/z)) w lp t sprintf('t=%1.2e',2**i), 3e10/x**2}; set title "raw data"; plot [:1000][1e6:] for[i=14:19] 'sofq_samples.dat' index i u ($0):($2) w lp t sprintf('t=%1.2e',2**i), 3e10/x**2; unset multi
