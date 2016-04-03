usage="$0 <nants> <num-iterations> <number-runs> <evaporation-constant> <input-file>"
case $# in
[01234]) echo  "wrong number of arguments  
		Usage is :$usage "; exit 2  >> /dev/null;;
esac


make cx
inpfile=$5

./cln.sh $inpfile
echo "#define ANTS  $1 
#define ncmax $2 
#define runs $3
#define gamma $4 
#define rho    $4 
" > par.h 
make 
comp=$?
echo -e "\n comp is $comp" 
test 0 -eq $comp || ( echo -e "\n the file was not compiled properly \n program is exiting \n\n "; exit 2 )
test  0 -eq $comp  && ( echo "the file compiled properly " ; ./aco  $inpfile  ;  ./gplot.sh $inpfile ;  latex   $inpfile.table  ; dvipdf $inpfile.table.dvi  )# ; d $inpfile) )#
comp=$?
test  0 -eq $comp  &&  echo -e "\n\n
ACO code by Saurabh Samdani \n 
Input  file was $inpfile \n
The code was file compiled properly and executed \n
Please check the output in $inpfile.table.pdf" 

