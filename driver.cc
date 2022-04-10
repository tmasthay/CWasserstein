#include "include/misfit.hh"
#include "include/cub.hh"
#include "include/wass_split.hh"
#include "include/sobolev.hh"
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include <rsf.h>

using T=float;

int main(int argc, char* argv[]){
    sf_init(argc, argv);

    //setup I/O files
    CUB f("in", "i"); f.headin(); //f.report();
    CUB g("data", "i"); g.headin(); //g.report();
    //CUB t("t", "i"); t.headin(); //t.report();
    //CUB p("p", "i"); p.headin(); //p.report();
    CUB output_file("out", "o"); output_file.setup(1);


    //Read axes from f
    sf_axis fa0 = f.getax(0); int nx = sf_n(fa0); T dx = sf_d(fa0);  
    sf_axis fa1 = f.getax(1); int nt = sf_n(fa1); T dt = sf_d(fa1); 
  

    //Read axes for g
    sf_axis ga0 = g.getax(0); int nxg = sf_n(ga0); T dxg = sf_d(ga0);
    sf_axis ga1 = g.getax(1); int ntg = sf_n(ga1); T dtg = sf_d(ga1);

    //Read axes for time t
    //sf_axis ta0 = t.getax(0); int nt_true = sf_n(ta0); T dt_true = sf_d(ta0);

    //Read axes for probability discretization p
    //sf_axis pa0 = t.getax(0); int np = sf_n(pa0); T dp = sf_d(pa0);

    //sanity assertions
    double eps=1e-5;
    assert( nx == nxg ); assert( abs(dx - dxg) < eps );
    assert( nt == ntg ); assert( abs(dt - dtg) < eps );
    //assert( nt == nt_true ); assert( abs(dt - dt_true) < eps );

    //read in data
    valarray<float> f_vec(0.0, nx * nt); f >> f_vec;
    valarray<float> g_vec(0.0, nx * nt); g >> g_vec;
    //valarray<float> t_vec(0.0, nt); t >> t_vec;
    //valarray<float> p_vec(0.0, np); p >> p_vec;

    //WassSplit<float> my_misfit(g_vec, t_vec, p_vec, nx);
    //float value = my_misfit.eval(f_vec);
    Sobolev<float> my_misfit(g_vec, (float) 0.0, nx, nt, dx, dt, (float) 0.0);

    float value = my_misfit.eval(f_vec);
  
    sf_axis output_axis = sf_maxa(1, (float) value, 0.0);
    output_file.putax(0, output_axis);

    cerr << "Final val = " << value << endl;
//    output_file << value;

    return 0;
}
