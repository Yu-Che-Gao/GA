cat $1.path.out > $1.constructiongraph.gnu ; 
./constructiongraph.sh $1;
gnuplot $1.constructiongraph.gnu 
