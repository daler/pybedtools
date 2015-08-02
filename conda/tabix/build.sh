set -x
echo $(pwd)
make
BIN=$PREFIX/bin
mkdir -p $BIN
cp tabix bgzip tabix.py $BIN
