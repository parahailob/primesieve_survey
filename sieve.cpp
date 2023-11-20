#include <vector>
#include <cmath>
#include <iostream>
#include <numeric>


std::vector<uint32_t> naive(uint32_t n){
    n++;
    std::vector<bool> is_prime(n, true);
    std::vector<uint32_t> result;
    is_prime[0] = false;
    is_prime[1] = false;
    
    for(uint32_t i = 2; i < n; i += 1){
        if(is_prime[i]){
            result.push_back(i);
            for(uint32_t j = 2*i; j < n; j += i){
                is_prime[j] = false;
            }
        }
    }

    return result;
}


std::vector<uint32_t> Eratosthenes(uint32_t n){
    n++;
    std::vector<bool> is_prime(n, true);
    is_prime[0] = false;
    is_prime[1] = false;
    
    // eliminating all even numbers
    for(uint32_t i = 4; i < n; i += 2) is_prime[i] = false;

    // for each prime, eliminate all of its multiples starting from its square
    for(uint32_t i = 3; i < std::sqrt(n); i += 2){
        if(is_prime[i]){
            for(uint32_t j = i*i; j < n; j += i){
                is_prime[j] = false;
            }
        }
    }
    std::vector<uint32_t> result;
    for(uint32_t i = 2; i < n; i++){
        if(is_prime[i]) result.push_back(i);
    }

    return result;
}


std::vector<uint32_t> Sundaram(uint32_t n){
    uint32_t k = (n-3) / 2 + 1;

    std::vector<bool> is_prime(k, true);
    for(uint32_t i = 0; i < (((int)std::sqrt(n) - 3) / 2 + 1); i++){
        uint32_t p = 2 * i + 3;
        uint32_t s = (p * p - 3) / 2;

        for(uint32_t j = s; j < k; j += p){
            is_prime[j] = false;
        }
    }
    std::vector<uint32_t> result = {2};
    for(uint32_t i = 0; i < k; i++){
        if(is_prime[i]) result.push_back(2*(i+1) + 1);
    }

    return result;
}


std::vector<uint32_t> Mairson(uint32_t n){
    // creating doubly linked list
    std::vector<uint32_t> RLINK(n+1);
    std::vector<uint32_t> LLINK(n+1);
    for(uint32_t i = 1; i < n; i++) RLINK[i] = i+1;
    RLINK[n] = 0;
    for(uint32_t i = 2; i <= n; i++) LLINK[i] = i-1;
    LLINK[1] = 0;
    std::vector<uint32_t> del(n+1, 0); // array of numbers to be eliminated each cycle
    uint32_t prime = 2; int factor = 2;
    while(prime <= std::sqrt(n)){
        uint32_t pointer = 0;
        while(prime * factor <= n){
            pointer++;
            del[pointer] = prime * factor;
            factor = RLINK[factor];
        }
        for(uint32_t i = 1; i <= pointer; i++){
            RLINK[LLINK[del[i]]] = RLINK[del[i]];
            LLINK[RLINK[del[i]]] = LLINK[del[i]];
        }
        prime = RLINK[prime];
        factor = prime;
    }

    // output the primes
    std::vector<uint32_t> result;
    uint32_t p = RLINK[1];
    while(p != 0){
        result.push_back(p);
        p = RLINK[p];
    }

    return result;
}


std::vector<uint32_t> Luo(uint32_t n){
    std::vector<bool> S(n - (uint32_t(n/2)+uint32_t(n/3)-uint32_t(n/6)) - 1, true);
    // does not contain multiples of 2 and 3
    uint32_t c = 0, k = 1, t = 2, q = uint32_t(std::sqrt(n) / 3), M = n / 3;
    for(uint32_t i = 1; i <= q; i++){
        k = 3-k;
        c = c + 4*k*i;
        uint32_t j = c;
        uint32_t ij = 2*i*(3 - k) + 1;
        t = t + 4*k;
        while(j <= M){
            S[j-1] = false;
            j = j + ij;
            ij = t - ij; 
        }
    }

    std::vector<uint32_t> result;
    if(n > 1) result.push_back(2);
    if(n > 2) result.push_back(3);
    for(uint32_t i = 0; i < S.size(); i++){
       if(S[i]){
        if(i % 2 == 0){
            result.push_back(5 + 6*uint32_t(i / 2));
        } else {
            result.push_back(5 + 6*uint32_t(i / 2) + 2);
        }
       }
    }

    return result;
}


std::vector<uint32_t> Gries_Misra(uint32_t n){
    std::vector<uint32_t> s(n+2);
    std::iota(s.begin(), s.end(), 1);
    uint32_t p = 2;
    uint32_t q; long long x;
    while(p <= std::sqrt(n)){
        q = p;
        while(p*q <= n){
            x = p*q;
            while(x <= n){
                uint32_t pred = std::min(x-1, (long long)s[x-1]);
                s[pred] = s[x];
                s[s[x]-1] = pred;
                x = p*x;
            }
            q = s[q];
        }
        p = s[p];
    }

    std::vector<uint32_t> result;
    uint32_t i = 2;
    while(i <= n){
        result.push_back(i);
        i = s[i];
    }

    return result;
}



void Extend (uint32_t w[], uint32_t &w_end, uint32_t &length, uint32_t n, bool d[], uint32_t &w_end_max) {
    /* Rolls full wheel W up to n, and sets length=n */
    uint32_t i, j, x;
    i = 0; j = w_end;
    x = length + 1; /* length+w[0] */
    while (x <= n) {
        w[++j] = x; /* Append x to the ordered set W */
        d[x] = false;
        x = length + w[++i];
    }
    length = n; w_end = j;
    if (w_end > w_end_max) w_end_max = w_end;
}

void Delete (uint32_t w[], uint32_t length, uint32_t p, bool d[], uint32_t &imaxf) {
    /* Deletes multiples p*w[i] of p from W, and sets imaxf to last i for deletion */
    uint32_t i, x;
    i = 0;
    x = p; /* p*w[0]=p*1 */
    while (x <= length) {
        d[x] = true; /* Remove x from W; */
        x = p*w[++i];
    }
    imaxf = i-1;
}

void Compress(uint32_t w[], bool d[], uint32_t to, uint32_t &w_end) {
    /* Removes deleted values in w[0..to], and if to=w_end, updates w_end, otherwise pads with zeros on right */
    uint32_t i, j;
    j = 0;
    for (i=1; i <= to; i++) {
        if (!d[w[i]]) {
            w[++j] = w[i];
        }
    }
    if (to == w_end) {
        w_end = j;
    } else {
        for (uint32_t k=j+1; k <= to; k++) w[k] = 0;
    }
}

std::vector<uint32_t> Pritchard(uint32_t N) {
    /* finds the nrPrimes primes up to N, printing them if printPrimes */
    std::vector<uint32_t> result;
    uint32_t *w = new uint32_t[N/4+5];
    bool *d = new bool[N+1];
    uint32_t w_end, length;
    /* representation invariant (for the main loop): */
    /* if length < N (so W is a complete wheel), w[0..w_end] is the ordered set W; */
    /* otherwise, w[0..w_end], omitting zeros and values w with d[w] true, is the ordered set W, */
    /* and no values <= N/p are omitted */
    uint32_t w_end_max, p, imaxf;
    /* W,k,length = {1},1,2: */
    w_end = 0; w[0] = 1;
    w_end_max = 0;
    length = 2;
    /* Pr = {2}: */
    result.push_back(2);
    p = 3;
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr = the first k primes */
    /* (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p*p <= N) {
        /* Append p to Pr: */
        result.push_back(p);
        if (length < N) {
            /* Extend W with length to minimum of p*length and N: */
            Extend (w, w_end, length, std::min(p*length,N), d, w_end_max);
        }
        Delete(w, length, p, d, imaxf);
        Compress(w, d, (length < N ? w_end : imaxf), w_end);
        /* p = next(W, 1): */
        p = w[1];
        if (p == 0) break; /* next p is after zeroed section so is too big */
        /* k++ */
    }
    if (length < N) {
        /* Extend full wheel W,length to N: */
        Extend (w, w_end, length, N, d, w_end_max);
    }
    /* gather remaining primes: */
    for (uint32_t i=1; i <= w_end; i++) {
        if (w[i] == 0 || d[w[i]]) continue;
        result.push_back(w[i]);
    }

    return result;
}
