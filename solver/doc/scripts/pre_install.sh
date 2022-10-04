#!/bin/sh

cmake . -B ./build -G "Ninja Multi-Config"
cmake --build ./build --config Release
pushd solver/doc
doxygen
ls
ls xml
popd
