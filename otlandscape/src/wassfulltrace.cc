#include "include/cub.hh"
#include "wasstrace.hh"
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include <rsf.h>

int main(int argc, char* argv[]){    
    sf_init(argc, argv);

    //setup I/O files
    CUB f("in", "i"); f.headin(); //f.report();
    CUB g("data", "i"); g.headin(); //g.report();
    CUB t("t", "i"); t.headin(); //t.report();
    CUB p("p", "i"); p.headin(); //p.report();

    int verbose, mode;
    float c;
    if(!sf_getint("v", &verbose)) verbose = 0; 
    if(!sf_getint("mode", &mode)) mode = 0;
    if(!sf_getfloat("c", &c)) c = 1.0;

    //Read axes from f
    sf_axis fa0 = f.getax(0); int nt = sf_n(fa0); float dt = sf_d(fa0);
    sf_axis fa1 = f.getax(1); int nx = sf_n(fa1); float dx = sf_d(fa1);

    //Read axes for g
    sf_axis ga0 = g.getax(0); int ntg = sf_n(ga0); float dtg = sf_d(ga0);
    sf_axis ga1 = g.getax(1); int nxg = sf_n(ga1); float dxg = sf_d(ga1); 
 
    //Read axes for time t
    sf_axis ta0 = t.getax(0); int nt_true = sf_n(ta0); float dt_true = sf_d(ta0);

    //Read axes for probability discretization p
    sf_axis pa0 = p.getax(0); int np = sf_n(pa0); float dp = sf_d(pa0);
 
    //sanity assertions
    double eps=1e-5;
    cerr << "nt,ntg=" << nt << "," << ntg << endl;
    assert( nt == ntg ); assert( abs(dt - dtg) < eps );
    assert( nx == nxg ); assert( abs(dx - dxg) < eps );
//    assert( n_cases == ng_cases ); 

    //float** vals;
    //vals = sf_floatalloc2(zs, xs);
    float* vals;
    vals = sf_floatalloc(nx);
    sf_file output_file;
    output_file = sf_output("out");
    sf_putint(output_file, "n1", nx);

    valarray<float> t_vec(0.0, nt); t >> t_vec;
    valarray<float> p_vec(0.0, np); p >> p_vec;

    for(int i = 0; i < nx; i++){
        valarray<float> g_vec(0.0, nt); g >> g_vec;
        valarray<float> f_vec(0.0, nt); f >> f_vec;
        //read in reference data
        vals[i] = wasstrace<float>(f_vec, g_vec, t_vec, p_vec, mode, c);
    }
 
    //vals[0] = value;
    sf_floatwrite(vals, nx, output_file);

    free(vals);
    exit(0);
}
