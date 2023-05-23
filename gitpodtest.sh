# unit test setup for gitpod.io
cd hypercomplex
mkdir ../.test/unit/hypercomplex
cp Hypercomplex.hpp ../.test/unit/hypercomplex/Hypercomplex.hpp
cp Polynomial.hpp ../.test/unit/hypercomplex/Polynomial.hpp
cd ../.test/unit
g++ -O0 -Wall --std=c++17 -o test test.cpp -lmpfr -lgmp
./test -d yes -w NoAssertions --use-colour yes --benchmark-samples 100 --benchmark-resamples 100000
rm -rf hypercomplex test
