echo "searching HDF5 files (will be removed)"
find . -name "*.h5"
while true; do
    read -p "Do you ve to remove them? (y/n)" yn
    case $yn in	
        [Yy]* ) find . -name "*.h5" -exec rm {} \;; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

