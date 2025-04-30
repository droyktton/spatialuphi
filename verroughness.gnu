lista="256 512 1024 2048 4096 8192 16384 32768 65536"; 


set multi lay 1,2; 

set xtics ("0.1" 0.1, "1" 1, "10" 10, "100" 100, "1000" 1000, "10000" 10000)
set ytics ( "0.1" 0.1, "10" 10, "100" 100, "1000" 1000, "10000" 10000)
set tit 'raw data'; 
set xla 't'; set yla 'sigma^2'; 
plot [][:1000] for[L in lista] sprintf('L%s/logcm_samples.dat',L) u 1:11:12 w lp t "sigma^2{u}, L=".L, \
for[z in "1.25 1.5 1.75"] 0.13*x**(1.0/z) t "zeta=1/2, z=".z lw (z=="1.5")?(4):(1); 


set xtics ("0.00001" 0.00001, "0.0001" 0.0001, "0.001" 0.001, "0.01" 0.01, "0.1" 0.1, "0.5" 0.5, "1" 1, "5" 5)
set ytics ("0.00001" 0.00001, "0.0001" 0.0001, "0.001" 0.001, "0.01" 0.01, "0.1" 0.1, "0.5" 0.5, "1" 1, "5" 5)
set xla 't/L^z'; 
set yla 'sigma^2/L^{zeta}'; 
set tit "zeta=1/2, z=3/2"; 
plot [][0.000001:0.1] for[L in lista] sprintf('L%s/logcm_samples.dat',L) u ($1/L**1.5):(($1>2)?($11/L):(1./0)):($12/L) w error t "sigma^2{u}, L=".L, \
0.15*x**(2*0.5/1.5) t 't^{2zeta_{KPZ}/z_{KPZ}}' lw 2, 0.01*x**(2*0.5/2.0) t 't^{2zeta_{EW}/z_{EW}}' lw 2

unset multi
