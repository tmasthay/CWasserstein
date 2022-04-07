#include "include/misfit.hh"
#include "include/cub.hh"
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>

using T=double;

int main(int argc, char* argv[]){
    sf_init(argc, argv);

    //setup I/O files
    CUB f("in", "i"); f.headin(); f.report();
    CUB g("data", "i"); g.headin(); g.report();
    CUB t("t", "i"); t.headin(); t.report();
    CUB output_file("out", "o"); output_file.setup(1);

    //Read axes from f
    sf_axis fa0 = f.getax(0); int nz = sf_n(fa0); T dz = sf_d(fa0);
    sf_axis fa1 = f.getax(1); int nx = sf_n(fa1); T dx = sf_d(fa1);
    sf_axis fa2 = f.getax(2); int nt = sf_n(fa2); T dt = sf_d(fa2);

    //Read axes for g
    sf_axis ga0 = g.getax(0); int nzg = sf_n(ga0); T dzg = sf_d(ga0);
    sf_axis ga1 = g.getax(1); int nxg = sf_n(ga1); T dxg = sf_d(ga1);
    sf_axis ga2 = g.getax(2); int ntg = sf_n(ga2); T dtg = sf_d(ga2);

    //Read axes for time t
    sf_axis ta0 = t.getax(0); int nt_true = sf_n(ta0); T dt_true = sf_d(ta0);

    //sanity assertions
    double eps=1e-10;
    assert( nz == nzg ); assert( abs(dz - dzg) < eps );
    assert( nx == nxg ); assert( abs(dx - dxg) < eps );
    assert( nt == ntg ); assert( abs(dt - dtg) < eps );
    assert( nt == nt_true ); assert( abs(dt - dt_true) < eps );

    return 0;
}
