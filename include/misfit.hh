#ifndef MISFIT_H
#define MISFIT_H

#include<valarray>
#include <vector>

using namespace std;

//throw this in a namespace;
//template<class T> using Ctn=valarray<T>;

template< class T >
class Misfit {
public:
        virtual T eval(const valarray<T> &m) = 0; //purely virtual --> mandatory
        virtual T grad(valarray<T> m, 
            valarray<T> m_hat) {}; //optional
        virtual T hessian(valarray<T> m, 
           valarray<T> m_tilde, 
           valarray<T> m_hat) {}; //optional
};

#endif
