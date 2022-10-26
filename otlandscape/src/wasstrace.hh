#include "include/cub.hh"
#include "include/wassall.hh"
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include <rsf.h>


float wasstrace(const valarray<float> &f,
    const valarray<float> &g, 
    const valarray<float> &t,
    const valarray<float> &p,
    int mode,
    float c=1.0){

    fprintf(stderr, 
        " CALLING WASSTRACE (%d,%d,%d,%d) ",
        f.size(), g.size(), t.size(), p.size());
    //compute distances
    float value = 0.0;
    if( mode == 0 ){
        for(int i = 0; i < t.size(); i++)
            value += pow(f[i]-g[i], 2.0);
    }
    else if( mode >= 1 ){
       if( mode == 1 ){
           WassSplit<float> my_misfit(g, t, p);
           my_misfit.set_num_splits(2);
           value = my_misfit.eval(f);
       }
       else if( mode == 2 ){
           WassSquare<float> my_misfit(g, t, p);
           value = my_misfit.eval(f);
       }
       else if( mode == 3 ){
           WassLinExp<float> my_misfit(g, t, p);
           my_misfit.set_sharpness(c);
           value = my_misfit.eval(f);
       }
       else if( mode == 4 ){
           //WassLin<float> my_misfit(g, t, p, 1);
           //value = my_misfit.eval(f);
           value = 0.0;
       }
       else if( mode == 5 ){
           WassExp<float> my_misfit(g, t, p);
           value = my_misfit.eval(f);
       }
       else{
           cerr << "Mode " << mode << " not supported in wasstrace.cc\n";
           exit(-2);
       }
    }
    else{
       cerr << "Mode " << mode << " not supported in wasstrace.cc\n";
       exit(-2);
    }
    fprintf(stderr, " RETURNING %f ", value);
    return value;
}
