rm -rf build;
mkdir build;
cd build;
cmake ..;
make;
mv src/huygens ..;
cd ..;
rm -rf build;