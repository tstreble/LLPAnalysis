paths=$(echo $PATH | tr ":" "\n") # split on :
pathunset=true
for pp in $paths
do
    if [[ $pp == *"LLPAnalysis/bin"* ]] ; then
        echo "Warning: path to LLP binaries already set to: $pp"
        echo "         ... Ignoring this setup command"
        pathunset=false
    fi
done

if [ "$pathunset" = true ] ; then
    export THISDIR=`pwd`
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THISDIR}/lib

    if [ -n "${DYLD_LIBRARY_PATH}" ] ; then
    export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THISDIR}/lib
    fi

    export PATH=${PATH}:${THISDIR}/bin
fi