#ifndef WASS_ALL_H
#define WASS_ALL_H

#include "wassslicer.hh"
#include <vector>
#include <cassert>

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSplit2 : public WassSlicer<T>{
private:
    int num_dists = 2;
public:
    using WassSlicer<T>::WassSlicer;

    T renorm_op(T x, int i_dist){
        if( i_dist == 0 ) return x >= 0 ? x : 0.0;
        else return x <= 0? -x : 0.0;
    }
};

template<class T>
class WassSquare : public WassSlicer<T>{
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;

    T renorm_op(T x, int i_dist){
        return x * x;
    }
};

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassLinExp : public WassSlicer<T>{
protected:
    T c = 1.0;
    T c_inv = 1.0;
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;
    
    void set_sharpness(T x) { this-> c = x; this->c_inv = 1.0 / this->c; }

    T renorm_op(T x, int i_dist){
        if( x >= 0 ) return x + c_inv;
        else return c_inv * exp(c * x);
    }
};

template<class T>
class WassLin : public WassSlicer<T>{
protected:
    valarray<T> c;
    T curr_shift;
    bool c_init = false;
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;

    void set_shift(T x) { this->c = x; }
    void init_c(){
        if(c_init) return;
        valarray<T> tmp(0.0, this->num_rec);
        c = tmp;
        c_init = true;
    }

    bool prescan_hyper(const valarray<T>& f){
        assert( f.size() == this->data.size() );
        //initialize c if it is necessary
        init_c();

        assert( this->c.size() > 0 ); 
        for(int i_r=0; i_r < this->num_rec; i_r++){
            for(int i_t=0; i_t < this->nt; i_t++){
                const int i_gbl = i_t + i_r * this->nt;
                T curr = min<T>(f[i_gbl], this->data[i_gbl]); 
                if( curr < c[i_r] ) c[i_r] = curr;
            } 
        } 
    }
    //define virtual function
    void set_trace_info(int trace_no){
        curr_shift = -1.0 * c[trace_no];
        cerr << "current shift == " << curr_shift << endl;
    }    

    T renorm_op(T x, int i_dist){
        return x + this->curr_shift;
    }
};

template<class T>
class WassExp : public WassSlicer<T>{
protected:
    T c = 1.0;
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;

    void set_sharpness(T x) { this-> c = x; }

    T renorm_op(T x, int i_dist){
        return exp(c*x);
    }
};
#endif
