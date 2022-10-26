#include "include/cub.hh"
#include "include/wassall.hh"
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
    sf_axis fa0 = f.getax(0); int nt = sf_n(fa0); T dt = sf_d(fa0);  

    //Read axes for g
    sf_axis ga0 = g.getax(0); int ntg = sf_n(ga0); T dtg = sf_d(ga0);

    //Read axes for time t
    sf_axis ta0 = t.getax(0); int nt_true = sf_n(ta0); T dt_true = sf_d(ta0);

    //Read axes for probability discretization p
    sf_axis pa0 = p.getax(0); int np = sf_n(pa0); T dp = sf_d(pa0);
 
    //sanity assertions
    double eps=1e-5;
    cerr << "nt,ntg=" << nt << "," << ntg << endl;
    assert( nt == ntg ); assert( abs(dt - dtg) < eps );
//    assert( n_cases == ng_cases ); 

    //float** vals;
    //vals = sf_floatalloc2(zs, xs);
    float* vals;
    vals = sf_floatalloc(1);
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
    valarray<float> g_vec(0.0, nt); g >> g_vec;
    valarray<float> f_vec(0.0, nt); f >> f_vec; 
    valarray<float> t_vec(0.0, nt); t >> t_vec;
    valarray<float> p_vec(0.0, np); p >> p_vec;

    fprintf(stderr, "nt = %d", nt);

    //compute distances
    float value = 0.0;
    if( mode == 0 ){
        /*
        bool do_renormalization = false;
        float s;
        if(!sf_getfloat("s",&s)) s = -1.0;
        cerr << "HELLO\n";
        Sobolev<float> my_misfit(g_vec, 
            s, 
            nx, 
            nt, 
            dx, 
            dt, 
            (float) 0.0,
            do_renormalization);
        value = my_misfit.eval(f_vec);
        */
        float sum = 0.0;
        for(int i = 0; i < nt; i++)
            sum += pow(f_vec[i]-g_vec[i], 2.0);
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
    vals[0] = value;
    sf_floatwrite(&vals[iz][ix], 1, output_file);

    free(vals);
    exit(0);
}
