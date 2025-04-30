## Runs many seed for one field

for((seed=1234;seed<=1234+20;seed++))
do 
	echo "seed="$seed; 
	./spatialuphi 131072 1.0005 1.0005 1.01 1000000 $seed; 
	mv inst_sofq.dat "inst_sofq"$seed".dat"; 
	mv logcm.dat "logcm"$seed".dat"; 
done
