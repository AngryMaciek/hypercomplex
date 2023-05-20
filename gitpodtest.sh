# unit test setup for gitpod.io
cd hypercomplex
mkdir ../.test/unit/hypercomplex
cp Hypercomplex.hpp ../.test/unit/hypercomplex/Hypercomplex.hpp
cp Polynomial.hpp ../.test/unit/hypercomplex/Polynomial.hpp
cd ../.test/unit
clang++ -O0 -Wall --std=c++17 -o test test.cpp -lmpfr -lgmp -fconstexpr-steps=4294967295
./test -d yes -w NoAssertions --use-colour yes --benchmark-samples 100 --benchmark-resamples 100000
rm -rf hypercomplex test
