#ifndef MISFIT_H
#define MISFIT_H

#include <valarray>
#include <vector>

using namespace std;
template<class T> using Ctn=vector<T>;

template< class T >
class Misfit {
public:
        virtual T eval(Ctn<T> m) = 0; //purely virtual --> mandatory
        virtual T grad(Ctn<T> m, 
            Ctn<T> m_hat) {}; //optional
        virtual T hessian(Ctn<T> m, 
           Ctn<T> m_tilde, 
           Ctn<T> m_hat) {}; //optional
};

#endif
