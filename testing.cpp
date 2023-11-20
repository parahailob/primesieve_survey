#include "sieve.cpp"
#include <chrono>
#include <iostream>
#include <functional>
#include <chrono>

using namespace std;

double test_sieve(function<std::vector<uint32_t>(uint32_t)> sieve, uint32_t n){
    double cycles = 10;
    auto start = chrono::high_resolution_clock::now();
    // cout << "started" << endl;
    for (int i = 0; i < cycles; i++) sieve(n);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    return (double)(duration.count() / cycles / 1000);
}


int main(){
    // vector<uint32_t> res = Gries_Misra(100);
    // for(uint32_t v: res) cout << v << endl;
    for(uint32_t n = 1000000; n <= 1000000000; n *= 10){
        cout << "n = " << n << '\n';
        cout << "Naive method: " << test_sieve(naive, n) << " ms\n";
        cout << "Luo method: " << test_sieve(Luo, n) << " ms\n";
        cout << "Eratosthenes method: " << test_sieve(Eratosthenes, n) << " ms\n";
        cout << "Sundaram method: " << test_sieve(Sundaram, n) <<  " ms\n";
        cout << "Gries-Misra method: " << test_sieve(Gries_Misra, n) << " ms\n";
        cout << "Pritchard method: " << test_sieve(Pritchard, n) << " ms\n";
        cout << "Mairson method: " << test_sieve(Mairson, n) << " ms\n";
    }
}