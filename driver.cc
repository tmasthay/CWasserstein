#include "include/cub.hh"
#include "include/wass_all.hh"
#include "include/sobolev.hh"
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include <rsf.h>


int main(int argc, char* argv[]){
    sf_init(argc, argv);

    using T=float;

    //setup I/O files
    CUB f("in", "i"); f.headin(); //f.report();
    CUB g("data", "i"); g.headin(); //g.report();
    CUB t("t", "i"); t.headin(); //t.report();
    CUB p("p", "i"); p.headin(); //p.report();

    int verbose;
    if(!sf_getint("v", &verbose)) verbose = 0; 
    //get mode of execution
    int mode;
    if(!sf_getint("mode", &mode)) mode = 0;

    //Read axes from f
    sf_axis fa0 = f.getax(0); int nx = sf_n(fa0); T dx = sf_d(fa0);  
    sf_axis fa1 = f.getax(1); int nt = sf_n(fa1); T dt = sf_d(fa1);
    sf_axis fa2 = f.getax(2); int n_cases = sf_n(fa2); T dc = sf_d(fa2); 

    //Read axes for g
    sf_axis ga0 = g.getax(0); int nxg = sf_n(ga0); T dxg = sf_d(ga0);
    sf_axis ga1 = g.getax(1); int ntg = sf_n(ga1); T dtg = sf_d(ga1);
//    sf_axis ga2 = g.getax(2); int ng_cases = sf_n(ga2); T dgc = sf_d(ga2);

    //Read axes for time t
    sf_axis ta0 = t.getax(0); int nt_true = sf_n(ta0); T dt_true = sf_d(ta0);

    //Read axes for probability discretization p
    sf_axis pa0 = p.getax(0); int np = sf_n(pa0); T dp = sf_d(pa0);

    int zs, xs;
    float bz, ez, bx, ex;
    int even_shift = (int) sqrt(n_cases);

    //parse arguments
    if(!sf_getint("zs", &zs)) zs = even_shift;
    if(!sf_getint("xs", &xs)) xs = even_shift;
    if(!sf_getfloat("bz", &bz)) bz = 0.0;
    if(!sf_getfloat("bx", &bx)) bx = 0.0;
    if(!sf_getfloat("ez", &ez)) ez = 1.0;
    if(!sf_getfloat("ex", &ex)) ex = 1.0;
 
    float dz_src = (ez - bz) / ( (float) zs );
    float dx_src = (ex - bx) / ( (float) xs );

    //sanity assertions
    double eps=1e-5;
    assert( nx == nxg ); assert( abs(dx - dxg) < eps );
    assert( nt == ntg ); assert( abs(dt - dtg) < eps );
//    assert( n_cases == ng_cases ); 

    float** vals;
    vals = sf_floatalloc2(zs, xs);
    sf_file output_file;
    output_file = sf_output("out");
    sf_putint(output_file, "n1", zs);
    sf_putint(output_file, "n2", xs);
    sf_putint(output_file, "n3", 1);
    sf_putfloat(output_file, "o1", bz);
    sf_putfloat(output_file, "o2", bx);
    sf_putfloat(output_file, "d1", dz_src);
    sf_putfloat(output_file, "d2", dx_src);    

    //read in reference data
    valarray<float> g_vec(0.0, nx * nt); g >> g_vec;
    valarray<float> f_vec_full(0.0, nx * nt * n_cases); 
    f >> f_vec_full;
    valarray<float> f_vec(0.0, nx*nt);
    const size_t lengths[] = {nx*nt};
    const size_t strides[] = {1};
    const valarray<size_t> lengths_arr(lengths,1);
    const valarray<size_t> strides_arr(strides,1);

    valarray<float> t_vec(0.0, nt); t >> t_vec;
    valarray<float> p_vec(0.0, np); p >> p_vec;

    //compute distances
    for(int iz = 0; iz < zs; iz++){
        for(int ix = 0; ix < xs; ix++){ 
            cerr << ix << "," << iz << endl;
            float value = 0.0;
            int case_no = (ix + iz*xs)*nx*nt;
            gslice curr_slice(case_no,
                lengths_arr,
                strides_arr);
            f_vec = f_vec_full[curr_slice];
            if( mode == 0 ){
                bool do_renormalization = false;
                float s;
                if(!sf_getfloat("s",&s)) s = -1.0;
                Sobolev<float> my_misfit(g_vec, 
                    s, 
                    nx, 
                    nt, 
                    dx, 
                    dt, 
                    (float) 0.0,
                    do_renormalization);
                value = my_misfit.eval(f_vec);
            }
            else if( mode >= 1 ){ 
               if( mode == 1 ){
                   if(verbose) cerr << "Wasserstein splitting\n";
                   WassSplit2<float> my_misfit(g_vec, t_vec, p_vec, nx);
                   my_misfit.set_dists(2);
                   value = my_misfit.eval(f_vec);
               }
               else if( mode == 2 ){
                   if( verbose ) cerr << "Wasserstein squaring\n";
                   WassSquare<float> my_misfit(g_vec, t_vec, p_vec, nx);
                   value = my_misfit.eval(f_vec);
               }
               else if( mode == 3 ){
                   if( verbose ) cerr << "Linexp\n";
                   WassLinExp<float> my_misfit(g_vec, t_vec, p_vec, nx);
                   float c;
                   if(!sf_getfloat("c1", &c)) sf_error("Need c1 for linexp");
                   my_misfit.set_sharpness(c);
                   value = my_misfit.eval(f_vec);
               }
               else if( mode == 4 ){
                   if(verbose) cerr << "Linear\n";
                   WassLin<float> my_misfit(g_vec, t_vec, p_vec, nx);
                   value = my_misfit.eval(f_vec);
               }
               else if( mode == 5 ){
                   if(verbose) cerr << "Exponential\n";
                   WassExp<float> my_misfit(g_vec, t_vec, p_vec, nx);
                   float c;
                   if(!sf_getfloat("c1", &c)) sf_error("Need c1 for exponential"); 
                   value = my_misfit.eval(f_vec);
               }
               else{
                   cerr << "Mode " << mode << " not supported in driver.cc\n";
                   exit(-2);
               }
            }
            else{
                   cerr << "Mode " << mode << " not supported in driver.cc\n";
                   exit(-2);
            }
            vals[iz][ix] = value;
            sf_floatwrite(&vals[iz][ix], 1, output_file);
        }
    }
    cerr << "Wrote correctly\n";

    free(*vals); free(vals);
    exit(0);
}
