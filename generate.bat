@echo off

powershell rm -recurse -force ./build

mkdir build
cd build
	
cmake .. 
cmake --build . --config Release

cd ..

pause