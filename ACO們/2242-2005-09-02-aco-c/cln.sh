usage="$0 <filename whose output files are to be removed> "
case $# in
[0234]) echo  "wrong number of arguments  
		Usage is :$usage "; exit 2  >> /dev/null;;
esac

rm -v $1.*
