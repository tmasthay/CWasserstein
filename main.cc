#include "wass_split.hh"
#include "helper.hh"
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;
using namespace helper;
using number=double;
using Vec=vector<number>;

Vec gauss(number mu, number sig, number a, number b, int N){
    Vec v(N,0);
    number denominator_inv = 1.0 / (2 * sig * sig);
    for(int i = 0; i < N; i++){
        number x = a + i * (b-a) / N;
        number numerator = (x-mu)*(x-mu);
        v.at(i) = exp(-numerator * denominator_inv);
    }
    return v;
}

int main(void){
    //time information
    number a = -10;
    number b = 10;
    int nt = 1000;
    Vec t(linspace(a,b,nt));

    //probability sampling info
    int np=1000;
    Vec p(linspace(0.0, 1.0, np));

    //data info
    int num_rec = 10;
    number mu = 0.0;
    number sig = 1.0;
    number shift = 3.0;

    //vector initialization
    Vec data(num_rec * nt, 0);
    Vec test(num_rec * nt, 0);

    //initialization for constant trace
    Vec v(gauss(mu, sig, a, b, nt));
    Vec v_shift(gauss(mu+shift, sig, a,b,nt));

    //move data around properly
    for(int i_r = 0; i_r < num_rec; i_r++){
        for(int i_t = 0; i_t < nt; i_t++){
            data.at(i_t + i_r * nt) = v.at(i_t);
            test.at(i_t + i_r * nt) = v_shift.at(i_t);
        } 
    }

    cout << data.size() << endl;
    cout << num_rec * nt << endl;
    cout << num_rec * t.size() << endl;

    WassSplit<number> w(data, t, p, num_rec);

    number value = w.eval(test);

    cout << "W_2^2 = " << value << endl;

    return 0;
}
