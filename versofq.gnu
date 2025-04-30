lista="256 512 1024 2048 4096 8192 16384 32768 65536"; 

set multi lay 1,2; 

set logs; set xla "q t^{1/z}"; 
set yla "S(q,t)/t^{(1+2 zeta)/z}"; 
do for[z in "1.5"]{
    set tit "zeta=1/2, z=".z; 
    plot [:1e7][1e-2:] for[i=14:19] 'L65536/sofq_samples.dat' index i u ($0*(t=2**i)**(2*0.5/z)):($2/t**((1+2*0.5)/z)) w lp t sprintf('t=%1.2e',2**i), 3e10/x**2
}; 
set title "raw data"; plot [:1000][1e6:] for[i=14:19] 'L65536/sofq_samples.dat' index i u ($0):($2) w lp t sprintf('t=%1.2e',2**i), 3e10/x**2; 

unset multi
