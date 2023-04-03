make clean
make html
cp source/index_replace.html build/html/index.html
cp -r source/redirecting/2Dmaterials/ build/html/tutorials/
cp -r source/redirecting/bulkfluids/ build/html/tutorials
cp -r source/redirecting/freeenergy/ build/html/tutorials/
cp -r source/redirecting/montecarlo/ build/html/tutorials/
cp -r source/redirecting/reaxff/ build/html/tutorials/

