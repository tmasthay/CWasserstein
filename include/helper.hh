#include<iostream>
#include<vector>

using namespace std;
template<class T> using Ctn=vector<T>;

namespace helper{


//behaves just like Python np.linspace    
template<class T>
Ctn<T> linspace(T a, T b, int N){
    Ctn<T> v(N,0);
    v.at(0) = a;
    for(int i = 1; i < N; i++){
        v.at(i) = a + i * (b-a) / (N-1);
    }
    return v;
}


//take slices from indices
template<class T>
Ctn<T> slice(Ctn<T> v, Ctn<int> indices){
    Ctn<T> w(indices.size(), 0);
    for(int i = 0; i < indices.size(); i++){
        w.at(i) = v.at(indices.at(i));
    }
    return w;
}

//take contiguous slice without having build container object
template<class T>
Ctn<T> slice_range(Ctn<T> v, int start, int end){
    int N = end - start;
    Ctn<T> w(N, 0);
    for(int i = 0; i < N; i++){
        w.at(i) = v.at(start + i);
    }
    return w;
}
    
};
