
#!/bin/bash

for x in level3/water-adsorption-in-silica/*/plot*.ipynb;
do 
    jupyter-nbconvert --to script ${x}
done

for x in level3/water-adsorption-in-silica/*/plot*.py;
do 
    DIR="$(dirname "${x}")" ; FILE="$(basename "${x}")"
    cd ${DIR}
    python3 ${FILE}
    rm ${FILE}
    cd ../../../
done
