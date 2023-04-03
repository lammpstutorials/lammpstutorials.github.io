make clean
make html
cp ./source/index_replace.html ./build/html/index.html
cp -rf ./source/redirecting/2Dmaterials/ ./build/html/tutorials/
cp -rf ./source/redirecting/bulkfluids/ ./build/html/tutorials
cp -rf ./source/redirecting/freeenergy/ ./build/html/tutorials/
cp -rf ./source/redirecting/montecarlo/ ./build/html/tutorials/
cp -rf ./source/redirecting/reaxff/ ./build/html/tutorials/

