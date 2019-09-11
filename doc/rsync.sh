if [ -z $1 ]
then
    echo "Need to specify a host..."
else
    rsync -PaW --inplace --del doxygen-doc/* $1:compass-doc/
fi
