#
# boost or unboost code
#

# arg1: b - boost, u - unboost

declare -a fnames=("trainhmm.cpp" "HMMProblem.cpp" "HMMProblem.h")


if [ "$1" = "b" ]
then
    echo "Boosting..."
    # comment un-boost
	for i in {0..2}
	do
	    perl -i -p -e 's/^(.+)\/\/UNBOOST$/\/\/\1\/\/UNBOOST/g' ${fnames[i]} 
	done
	perl -i -p -e 's/^(.+)#UNBOOST$/#\1#UNBOOST/g' makefile
    # un-comment boost
	for i in {0..2}
	do
	    perl -i -p -e 's/^\/\/(.+)\/\/BOOST$/\1\/\/BOOST/g' ${fnames[i]} 
	done
	perl -i -p -e 's/^#(.+)#BOOST$/\1#BOOST/g' makefile
elif [ "$1" = "u" ]
then
    echo "Un-boosting"
    # un-comment un-boost
	for i in {0..2}
	do
	    perl -i -p -e 's/^\/\/(.+)\/\/UNBOOST$/\1\/\/UNBOOST/g' ${fnames[i]} 
	done
	perl -i -p -e 's/^#(.+)#UNBOOST$/\1#UNBOOST/g' makefile
    # comment boost
	for i in {0..2}
	do
	    perl -i -p -e 's/^(.+)\/\/BOOST$/\/\/\1\/\/BOOST/g' ${fnames[i]} 
	done
	perl -i -p -e 's/^(.+)#BOOST$/#\1#BOOST/g' makefile
else
    echo -e "Wrong parameter.\nUsage:\n\tboost [u|b]"
fi