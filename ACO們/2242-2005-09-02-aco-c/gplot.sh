./objconv.sh $1;
./varconv.sh $1 ;
./eval.sh $1 ;
./stddevplot.sh $1 ;
./iter.sh $1 ;
./runstats.sh $1 ;
./time.sh $1 ;
 gnuplot $1.obj.gnu ;
 gnuplot $1.var.gnu;
 gnuplot $1.dev.gnu;
 gnuplot $1.run.gnu;
 gnuplot $1.eval.gnu;
 gnuplot $1.iter.gnu;
 gnuplot $1.time.gnu
./path.sh $1 ;
 gnuplot $1.constructiongraph.gnu;

