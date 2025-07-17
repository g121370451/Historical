cmake.exe -DCMAKE_BUILD_TYPE=Debug --preset debug -S .\ -B .\cmake-build-debug-preset
cmake.exe --build .\cmake-build-debug-preset --target experiment_program -j 4