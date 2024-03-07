set -e

cd  semd
rm -rf out
mkdir out
cd out 
cmake ..
make

echo "Congratulations on the successful build of SEMD! The output path is out/tests/test_semd!"