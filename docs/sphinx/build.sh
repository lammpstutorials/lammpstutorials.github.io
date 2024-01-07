git submodule update  --remote

make clean
make html
cp source/index_replace.html build/html/index.html
