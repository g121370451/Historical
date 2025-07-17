cmake.exe -DCMAKE_BUILD_TYPE=Debug --preset release -S .\ -B .\cmake-build-release-preset
cmake.exe --build .\cmake-build-release-preset --target experiment_program -j 4