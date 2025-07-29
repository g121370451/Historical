cmake -DCMAKE_BUILD_TYPE=Debug --preset debug -S .\ -B ./cmake-build-debug-preset
cmake --build ./cmake-build-debug-preset --target experiment_program -j 4