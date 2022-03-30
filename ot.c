#include <rsf.h>
#include <unistd.h>

float my_max(float* x, int N){
    int i, max;
    max = x[0];
    for(i = 0; i < N; i++){
        if( x[i] > max ){
            max = x[i];
        }
    }
    return max;
}

float* linspace(float a, float b, int N){
    float *x, delta;
    int i;

    x = sf_floatalloc(N);
    x[0] = a;

    delta = (b-a) / (N-1);

    for(i = 1; i < N-1; i++){
        x[i] = x[i-1] +  delta;
    }
    x[N-1] = b;
    return x;
}

float integrate(float* f, float* x, int N){
    float a=x[0],b=x[N-1],sum=0.0,*gpts,*gwts;
    int order=5,ig,ix;

    float dist = 0.5 * (b-a);
    float mid = 0.5 * (b+a);

    gpts = sf_floatalloc(order);
    gwts = sf_floatalloc(order);

    gpts[0] = -.90618;
    gpts[1] = -0.538469;
    gpts[2] = 0.0;
    gpts[3] = -1.0 * gpts[1];
    gpts[4] = -1.0 * gpts[0];

    gwts[0] = 0.236927;
    gwts[1] = 0.478629;
    gwts[2] = 0.568889;
    gwts[3] = gwts[1];
    gwts[4] = gwts[0];

    for(ig = 0; ig < order; ig++){
        float curr_pt = dist * gpts[ig] + mid;
        //////fprintf(stderr, "ig = %d\n", ig);
        for(ix=1; ix < N; ix++){
            if( x[ix-1] <= curr_pt && curr_pt < x[ix] ){
                float t = (x[ix] - gpts[ig]) / (x[ix] - x[ix-1]);
                sum += gwts[ig] * (t * f[ix-1] + (1-t) * f[ix]);
                break;
            }
        }    
    }
    return dist * sum;
}

void renormalize(float* f, float *x, int N){
    float sum = 0.0;
    int i;
    for(i = 0; i < N - 1; i++){
        sum += 0.5 * (f[i+1] + f[i]) * (x[i+1] - x[i]);
        //sum += f[i];
    }
    if( sum == 0.0 ) {
        //////fprintf(stderr, "Unnormalizable...returning original\n");
        return;
    }
    float reciprocal = 1.0 / sum;
    //fprintf(stderr, "RECIP: %.16e...", reciprocal);
    float max = reciprocal * f[i];
    for( i = 0; i < N; i++ ){
        f[i] = reciprocal * f[i];
        if( f[i] > max ) max = f[i];
    }
    //fprintf(stderr, "max = %.16e\n", max);
}

void softplus_normalize(float *f, float* x, float b, int N){
    int i;
    float C = 0.0;
    for(i = 0; i < N; i++){
        f[i] = log(1.0 + exp(b * f[i]));
    }
    for(i = 0; i < N; i++){
        C += 0.5 * (f[i+1] + f[i]) * (x[i+1] - x[i]);
    }
    for(i = 0; i < N; i++){
       f[i] = f[i] / C;
    }
}

void split_normalize(float* f, 
    float* f_pos, 
    float* f_neg, 
    float* x,
    int N){
    int i;
    for(i = 0; i < N; i++){
        if( f[i] == 0.0 ){
            f_pos[i] = 0.0;
            f_neg[i] = 0.0;
        }
        f_pos[i] = f[i] > 0 ? f[i] : 0.0;
        f_neg[i] = f[i] < 0 ? -f[i] : 0.0;
        //////fprintf(stderr, "BEFORE RENORMALIZATION: (i, pos, neg) = (%d, %.16e, %.16e)\n", i, f_pos[i], f_neg[i]);
    }
    renormalize(f_pos, x, N);
    renormalize(f_neg, x, N);

    for(i = 0; i < N; i++){
        ////fprintf(stderr, "(i, pos, neg) = (%d, %.16e, %.16e)\n", i, f_pos[i], f_neg[i]);
    }
}

void cdf(float* f, float* x, float* F, int N){
    int i;
    for(i = 0; i < N; i++){
        if(i == 0){
            F[i] = 0.0;
        }
        else{
            F[i] = F[i-1] + 0.5 * (x[i]-x[i-1]) * (f[i] + f[i-1]);
        }
    }
}

void quantile(float* p, float* x, float* F, float* G, int Np, int Nx){
    int ip, ix, sx=1;
    float tau = 0.0;
    for(ip = 0; ip < Np; ip++){
        //int found = 0;

        for(ix = sx; ix < Nx; ix++){
            if( F[ix-1] <= p[ip] && p[ip] <= F[ix] ){
                sx = ix;
                if( abs(F[ix] - F[ix-1]) <= tau){
                    //fprintf(stderr, "HELLO\n");
                    G[ip] = 0.5 * (x[ix] + x[ix-1]);
                    //fprintf(stderr, "x[i]=%f, x[i-1]=%f\n", x[ix], x[ix-1]);
                    //////fprintf(stderr, "HERE!\n");
                }
                else{
                    float t = (F[ix] - p[ip]) / (F[ix] - F[ix-1]);
                    //fprintf(stderr, "t = %f\n", t);
                    G[ip] = t * x[ix-1] + (1-t) * x[ix];
                }
               // ////fprintf(stderr, "(ix, sx, p, G) = (%d, %d, %.8f, %.8f)\n", ix, sx, p[ip], G[ip]);
                break;
            }
        }
    }
}

void transport(float* f, float* g, float* T, float* x, int N){
    float *F, *G;

    F = sf_floatalloc(N);
    G = sf_floatalloc(N);

    cdf(f,x,F,N);
    cdf(g,x,G,N);

    quantile(F, x, G, T, N, N);
}

void wass_int(float* dest, float* f, float* g, float* x, float* p, int N, int P){
    float *F, *G, *Finv, *Ginv;
    int ip;

    F = sf_floatalloc(N);
    G = sf_floatalloc(N);
    Finv = sf_floatalloc(P);
    Ginv = sf_floatalloc(P);

    cdf(f,x,F,N);
    cdf(g,x,G,N);

    int i;
    // for(i = 0; i < N; i++){
    //    ////fprintf(stderr, "(i, f,g, F, G) = (%d, %.8f, %.8f, %.8f, %.8f)\n", i,f[i], g[i], F[i], G[i]);
    // }

    quantile(p, x, F, Finv, P, N);
    quantile(p, x, G, Ginv, P, N);

    for(i = 0; i < P; i++){
        //float diff = Finv[i] - Ginv[i];
        ////fprintf(stderr, "diff[%d] = %.15f\n", i, diff);
    }
    for(ip = 0; ip < P; ip++){
        dest[ip] = Finv[ip] - Ginv[ip];
        dest[ip] = dest[ip] * dest[ip];
        //fprintf(stderr, "dest[%d] = %.15e\n", ip, dest[ip]);
    }
}

float wass2(float* f, float* g, float* x, float*p, int N, int P){
    int ip;
    float *y, sum;

    y = sf_floatalloc(P);

    float f_int=0.0, g_int=0.0;
    int ii;
    for(ii = 0; ii < N - 1; ii++){
        float dx = x[ii+1] - x[ii];
        float avgf = 0.5 * (f[ii+1] + f[ii]);
        float avgg = 0.5 * (g[ii+1] + g[ii]);
        f_int += dx * avgf;
        g_int += dx * avgg;
    }
    //fprintf(stderr, "f_int = %f\n", f_int);
    //fprintf(stderr, "g_int = %f\n", g_int);
    float tol = 0.0;
    float support_penalty = 10.0;
    if( f_int <= tol && g_int <= tol ) {
        //fprintf(stderr, "Support complements agree\n");
        return 0.0;
    }
    if( (f_int <= tol && g_int > tol) || (f_int > tol && g_int <= tol) ){
        fprintf(stderr, "Support ISSUES\n");
        return support_penalty;
    }

    //transport(f, g, T, x, N);
    wass_int(y, f, g, x, p, N, P);

    //fprintf(stderr, "max(y) = %.15f", my_max(y, P));
    
    sum = 0.0;
    for(ip = 0; ip < P - 1; ip++){
        //float dp = p[ip+1] - p[ip];
        //float avg = y[ip+1] + y[ip];
        sum = sum + 0.5 * (p[ip+1] - p[ip]) * (y[ip+1] + y[ip]);
//        fprintf(stderr, "sum[%d] = %f\n", ip, sum);
        if( ip == P-2 ){
            //fprintf(stderr, "FINAL = %f\n", p[ip+1]);
        }

        //fprintf(stderr, "(dp, avg, sum) = (%.15f, %.15f, %.15f)\n", dp, avg, sum);
    }
    //sum = integrate(y, linspace(0.0, 1.0, 10*N), 10*N);
    //fprintf(stderr, "\n\n\n\nRETURNING: %f\n\n\n\n", sum);
    return sum;
}

float wass2_split(float* f, float* g, float* x, float* p, int N, int P){
    float *f_pos, *f_neg, *g_pos, *g_neg;
    float pos=0.0, neg=0.0;
    
    f_pos = sf_floatalloc(N);
    f_neg = sf_floatalloc(N);
    g_pos = sf_floatalloc(N);
    g_neg = sf_floatalloc(N);

    //int i;
    split_normalize(f, f_pos, f_neg, x, N);
    split_normalize(g, g_pos, g_neg, x, N);

    int i;
    for(i = 0; i < N; i++){
        //fprintf(stderr, "(i,f,g,f_pos,g_pos,f_neg,g_neg,pos_diff,neg_diff, tot_diff) = (%d, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e)\n", i, f[i], g[i], f_pos[i], g_pos[i], f_neg[i], g_neg[i], f_pos[i] - g_pos[i], f_neg[i] - g_neg[i], f[i] - g[i]);
    }

    pos = wass2(f_pos, g_pos, x, p, N, P);
    neg = wass2(f_neg, g_neg, x, p, N, P);

    //fprintf(stderr, "%.8f...%.8f\n", pos, neg);

    //////fprintf(stderr, "(+w,+mf,+mg,-w,-mf,-mg) = (%.15f,%.15f,%.15f,%.15f,%.15f,%.15f)\n", pos, my_max(f_pos,N), my_max(g_pos,N), neg, my_max(f_neg,N), my_max(g_neg,N));

    if( pos > 100 ) fprintf(stderr, "Pos = %f\n", pos);
    if( neg > 100 ) fprintf(stderr, "Neg = %f\n", neg);
    return pos + neg;
}

float wass2abs(float* f, float* g, float* x, float* p, int N, int P){
    int i;

    for(i = 0; i < N; i++){
        if( f[i] < 0.0 ) f[i] = -1.0 * f[i];
        if( g[i] < 0.0 ) g[i] = -1.0 * g[i];
    }

    renormalize(f,x,N);
    renormalize(g,x,N);

    return wass2(f,g,x,p,N,P);
}

float wass2trace(float*** f, float*** g, float* t, int nz, int nx, int nt, int np){
    //assume n3xn2xn1=(z,x,t) coordinates
    float sum=0.0, *p;
    int iz, ix;

    p = sf_floatalloc(np);
    //p = linspace(0.0, 1.0, np);
    p = linspace(0.01, 0.99, np);

    //fprintf(stderr, "initialized sum = %.16e", sum);
    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            float tmp = wass2_split(f[iz][ix], g[iz][ix], t, p, nt, np);
            sum += tmp;
            if( tmp > 100 ){
                //fprintf(stderr, "(iz,ix,tmp) = (%d, %d, %f)\n", iz,ix,tmp);
                sleep(1);
            }
            //fprintf(stderr, "(sum,curr) = (%f,%f)\n", sum, tmp);
        }
    }
    //fprintf(stderr, "wass2^2 = %.16e\n", sum);
    return sum;
}

float wass2traceabs(float*** f, float*** g, float* t, int nz, int nx, int nt, int np){
        //assume n3xn2xn1=(z,x,t) coordinates
    float sum=0.0, *p;
    int iz, ix;

    p = sf_floatalloc(np);
    p = linspace(0.0, 1.0, np);

    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            sum += wass2abs(f[iz][ix], g[iz][ix], t, p, nt, np);
        }
    }
    return sum;
}

float wass2traceslice(float*** f, float*** g, float* t, int* slice, int ns, int nx, int nt, int np){
    float ***f_sliced, ***g_sliced;
    int is,ix,it;

    f_sliced = sf_floatalloc3(ns, nx, nt);
    g_sliced = sf_floatalloc3(ns, nx, nt);

    for(is = 0; is < ns; is++){
        for(ix = 0; ix < nx; ix++){
            for(it = 0; it < nt; it++){
                ////fprintf(stderr, "(slice[is]-->is,ix,it,f_s, g_s) = (%d-->%d, %d, %d, %.4e, %.4e)", slice[is], is, ix, it, f[slice[is]][ix][it], g[slice[is]][ix][it]);
                f_sliced[is][ix][it] = f[slice[is]][ix][it];
                g_sliced[is][ix][it] = g[slice[is]][ix][it];
                ////fprintf(stderr, "     :::::     (slice[is]-->is,ix,it,f_s, g_s) = (%d-->%d, %d, %d, %.4e, %.4e)\n", slice[is], is, ix, it, f_sliced[is][ix][it], g_sliced[is][ix][it]);
            }
        }
    }
    return wass2trace(f_sliced, g_sliced, t, ns, nx, nt, np);
}

float wass2tracesliceabs(float*** f, float*** g, float* t, int* slice, int ns, int nx, int nt, int np){
    float ***f_sliced, ***g_sliced;
    int is,ix,it;

    f_sliced = sf_floatalloc3(ns, nx, nt);
    g_sliced = sf_floatalloc3(ns, nx, nt);

    for(is = 0; is < ns; is++){
        for(ix = 0; ix < nx; ix++){
            for(it = 0; it < nt; it++){
                f_sliced[is][ix] = f[slice[is]][ix];
                g_sliced[is][ix] = f[slice[is]][ix];
            }
        }
    }
    return wass2traceabs(f_sliced, g_sliced, t, ns, nx, nt, np);
}

float wass2tracesurf(float*** f, float*** g, float* t, int nx, int nt, int np){
    int *slice;

    slice = sf_intalloc(1);
    slice[0] = 0;

    return wass2traceslice(f, g, t, slice, 1, nx, nt, np);
}

float wass2tracesurfabs(float*** f, float*** g, float* t, int nx, int nt, int np){
    int *slice;

    slice = sf_intalloc(1);
    slice[0] = 0;

    return wass2tracesliceabs(f, g, t, slice, 1, nx, nt, np);
}

float wass2tracesurfsp(float*** f, float*** g, float* t, int nx, int nt, int np){
    float b = 200.0;
    float total = 0.0;
    float *p;
    p = sf_floatalloc(np);
    p = linspace(0.0, 1.0, np);

    int* slice;
    slice = sf_intalloc(1);
    slice[0] = 0;
   
    int ix, it;
    for(ix = 0; ix < nx; ix++){
        for(it = 0; it < nt; it++){
           softplus_normalize(f[0][ix], t, b, nt);
           softplus_normalize(g[0][ix], t, b, nt);
           total += wass2(f[0][ix], g[0][ix], t, p, nt, np); 
        }
    }
    return total; 
}

float l2surf(float*** f, float*** g, int nx, int nt){
    float sum=0.0;
    int nb=35;
    int ix,it;
    for(ix = 0; ix < nx; ix++){
        for(it = 0; it < nt; it++){
            float diff = f[0][ix][it] - g[0][ix][it];
            sum += diff * diff;
        }
    }
    return sum;
}

float l2total(float*** f, float*** g, int nz, int nx, int nt){
    float sum=0.0;

    int iz, ix,it;
    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            for(it = 0; it < nt; it++){
                float diff = f[iz][ix][it] - g[iz][ix][it];
                sum = sum + diff * diff;
            }
        }
    }
    return sum;
}

float wass2squaretotal(float*** f, float*** g,
                    float* t, float* p,
                    int nz, int nx, int nt, int np){
    float sum=0.0;
   
    int iz, ix, it;
    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            for(it = 0; it < nt; it++){
                f[iz][ix][it] = f[iz][ix][it] * f[iz][ix][it];
                g[iz][ix][it] = g[iz][ix][it] * g[iz][ix][it];
            }
            renormalize(f[iz][ix], t, nt);
            renormalize(g[iz][ix], t, nt); 
            sum += wass2(f[iz][ix], g[iz][ix], t, p, nt, np);
        }
    }
    return sum;
}

int main(int argc, char* argv[]){
    //assume structured grid
    int nz, nx, nt, nt_true, nz_check, nx_check, nt_check, mode;
    float ***f, ***g, *t, *p;

    sf_file f_file, g_file, t_file, wass_file;

    sf_init(argc, argv);

    //define inputs
    f_file = sf_input("in");
    g_file = sf_input("g");
    t_file = sf_input("t");

    //define outputs
    //T_file = sf_output("out");
    wass_file = sf_output("out");


    //read command-line inputs
    if (!sf_histint(f_file, "n1", &nz)) sf_error("No n3");
    if (!sf_histint(f_file, "n2", &nx)) sf_error("No n2");
    if (!sf_histint(f_file, "n3", &nt)) sf_error("No n1");
    if (!sf_histint(g_file, "n1", &nz_check)) sf_error("No n3: g");
    if (!sf_histint(g_file, "n2", &nx_check)) sf_error("No n2: g");
    if (!sf_histint(g_file, "n3", &nt_check)) sf_error("No n1: g");
    if (!sf_histint(t_file, "n1", &nt_true)) sf_error("No n1: t");
    
    //read command-line literal inputs
    if (!sf_getint("mode", &mode)) mode=1;

    // ////fprintf(stderr, "(nt,nt_check,nt_true) = (%d,%d,%d)\n", nt, nt_check, nt_true);

    //Input validation
    if( nz_check != nz ) sf_error("Dimension mismatch: z, (%d,%d)\n", nz_check, nz);
    if( nx_check != nx ) sf_error("Dimension mismatch: x, (%d,%d)\n", nz_check, nz);
    if( nt_check != nt ) sf_error("Dimension mismatch: t1 yo");
    if( nt != nt_true ) sf_error("Dimension mismatch: t2 yo, nt=%d, ntrue=%d", nt, nt_true);

    //allocate memory
    f = sf_floatalloc3(nz, nx, nt);
    g = sf_floatalloc3(nz, nx, nt);
    t = sf_floatalloc(nt);

    //read in data 
    int iz, ix;
    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            //////fprintf(stderr, "(iz,ix) = (%d, %d)\n", iz, ix);
            sf_floatread(f[iz][ix], nt, f_file);
            sf_floatread(g[iz][ix], nt, g_file);
        }
    }
    sf_floatread(t, nt, t_file);

    // int it;
    // for(iz = 0; iz < nz; iz++){
    //     for(ix = 0; ix < nx; ix++){
    //         for(it = 0; it < nt; it++){
    //             //float tmp = f[iz][ix][it] - g[iz][ix][it];
    //             ////fprintf(stderr, "(iz,ix,it,diff) = (%d,%d,%d,%.8f)\n", iz,ix,it,tmp);
    //         }
    //     }
    // }

    //transport(f,g,T,x,N);
    float distance;
    int np=nt;

    switch( mode ){
        case 1:
            distance = wass2tracesurf(f,g,t,nx,nt,np); break;
        case 2:
            distance = l2surf(f,g,nx,nt); break;
        case 3:
            distance = wass2trace(f,g,t,nz,nx,nt,np); break;
        case 4:
            distance = l2total(f,g,nz,nx,nt); break;
        case 5:
            distance = wass2tracesurfabs(f,g,t,nx,nt,np); break;
        case 6:
            distance = wass2tracesurfsp(f,g,t,nx,nt,np); break;
        case 7:
            distance = wass2squaretotal(f,g,t,p,nz,nx,nt,np); break;
        default:
            fprintf(stderr, "Case no. %d not supported. Exiting. \n", mode);
            exit(-1);
    }

    float* distance_tmp;
    int it;
    distance_tmp = sf_floatalloc(1);
    for(it = 0; it < 1; it++){
        distance_tmp[it] = distance;
    }
    sf_putint(wass_file, "n1", 1);
    sf_putint(wass_file, "n2", 1);
    sf_putint(wass_file, "n3", 1);
    sf_floatwrite(distance_tmp, 1, wass_file);

    int get_detailed_info = 0;
    
    if( get_detailed_info ) {
        fprintf(stderr, "Getting detailed info\n");
        sf_file p_file, cdffpos_file, cdffneg_file, cdfgpos_file, cdfgneg_file, ff_file, fpos_file, fneg_file, gg_file, gpos_file, gneg_file;
        sf_file qfpos_file, qfneg_file, qgpos_file, qgneg_file, integrandpos_file, integrandneg_file;
        float ***F_pos, ***F_neg, ***G_pos, ***G_neg, *f_pos, *f_neg, *g_pos, *g_neg, ***QFpos, ***QFneg, ***QGpos, ***QGneg, ***integrand_pos, ***integrand_neg;

        p_file = sf_output("p");
        cdffpos_file = sf_output("f_cdfpos");
        cdffneg_file = sf_output("f_cdfneg");
        cdfgneg_file = sf_output("g_cdfneg");
        cdfgpos_file = sf_output("g_cdfpos");
        ff_file = sf_output("ff");
        fpos_file = sf_output("fpos");
        fneg_file = sf_output("fneg");
        gg_file = sf_output("gg");
        gpos_file = sf_output("gpos");
        gneg_file = sf_output("gneg");

        
        //quantile files
        qfpos_file = sf_output("qfpos");
        qfneg_file = sf_output("qfneg");
        qgpos_file = sf_output("qgpos");
        qgneg_file = sf_output("qgneg");
        integrandpos_file = sf_output("integrand_pos");
        integrandneg_file = sf_output("integrand_neg");
 
        fprintf(stderr, "YO!!!\n");

        p = sf_floatalloc(np);
        //q = sf_floatalloc3(nz,nx,np);

        F_pos = sf_floatalloc3(nz,nx,nt);
        F_neg = sf_floatalloc3(nz,nx,nt);

        G_pos = sf_floatalloc3(nz,nx,nt);
        G_neg = sf_floatalloc3(nz,nx,nt);

        f_pos = sf_floatalloc(nt);
        f_neg = sf_floatalloc(nt);

        g_pos = sf_floatalloc(nt);
        g_neg = sf_floatalloc(nt);

        p = linspace(0.0, 1.0, np);

        QFpos = sf_floatalloc3(nz,nx,np);
        QFneg = sf_floatalloc3(nz,nx,np);
        QGpos = sf_floatalloc3(nz,nx,np);
        QGneg = sf_floatalloc3(nz,nx,np);

        integrand_pos = sf_floatalloc3(nz,nx,np);
        integrand_neg = sf_floatalloc3(nz, nx, np);


        int jz, jx;

        sf_putint(fpos_file, "n1", nt);
        sf_putint(fpos_file, "n2", nx);
        sf_putint(fpos_file, "n3", nz);
        sf_putint(fneg_file, "n1", nt);
        sf_putint(fneg_file, "n2", nx);
        sf_putint(fneg_file, "n3", nz);
 

        sf_putint(gpos_file, "n1", nt);
        sf_putint(gpos_file, "n2", nx);
        sf_putint(gpos_file, "n3", nz);
        sf_putint(gneg_file, "n1", nt);
        sf_putint(gneg_file, "n2", nx);
        sf_putint(gneg_file, "n3", nz);

        sf_putint(cdffpos_file, "n1", nt);
        sf_putint(cdffpos_file, "n2", nx);
        sf_putint(cdffpos_file, "n3", nz);
            
        sf_putint(cdffneg_file, "n1", nt);
        sf_putint(cdffneg_file, "n2", nx);
        sf_putint(cdffneg_file, "n3", nz);

        sf_putint(cdfgpos_file, "n1", nt);
        sf_putint(cdfgpos_file, "n2", nx);
        sf_putint(cdfgpos_file, "n3", nz);

        sf_putint(cdfgneg_file, "n1", nt);
        sf_putint(cdfgneg_file, "n2", nx);
        sf_putint(cdfgneg_file, "n3", nz);
             
        for(jz = 0; jz < nz; jz++){
            for(jx = 0; jx < nx; jx++){
                split_normalize(f[jz][jx], f_pos, f_neg, t, nt);
                split_normalize(g[jz][jx], g_pos, g_neg, t, nt);
        
                cdf(f_pos, t, F_pos[jz][jx], nt);
                cdf(f_neg, t, F_neg[jz][jx], nt);
        
                cdf(g_pos, t, G_pos[jz][jx], nt);
                cdf(g_neg, t, G_neg[jz][jx], nt);
        
                quantile(p, t, F_pos[jz][jx], QFpos[jz][jx], np, nt);
                quantile(p, t, F_neg[jz][jx], QFneg[jz][jx], np, nt);
                quantile(p, t, G_pos[jz][jx], QGpos[jz][jx], np, nt);
                quantile(p, t, G_neg[jz][jx], QGneg[jz][jx], np, nt);
        
                wass_int(integrand_pos[jz][jx], f_pos, g_pos, t,p, nt, np);
                wass_int(integrand_neg[jz][jx], f_neg, g_neg, t, p, nt, np);
        
                //fprintf(stderr, "HERE!!!\n");
                sf_floatwrite(F_pos[jz][jx], nt, cdffpos_file);
                //fprintf(stderr, "HERE FIRST!!!\n");
                sf_floatwrite(F_neg[jz][jx], nt, cdffneg_file);
                sf_floatwrite(G_pos[jz][jx], nt, cdfgpos_file);
                sf_floatwrite(G_neg[jz][jx], nt, cdfgneg_file);
                //fprintf(stderr, "Cumulative dist. works\n");
                sf_floatwrite(f[jz][jx], nt, ff_file);
                //fprintf(stderr, "Density writing works\n");
                sf_floatwrite(f_pos, nt, fpos_file);
                sf_floatwrite(f_neg, nt, fneg_file);
                sf_floatwrite(g[jz][jx], nt, gg_file);
                sf_floatwrite(g_pos, nt, gpos_file);
                sf_floatwrite(g_neg, nt, gneg_file);
                sf_floatwrite(p, np, p_file);
                sf_floatwrite(QFpos[jz][jx], np, qfpos_file);
                sf_floatwrite(QFneg[jz][jx], np, qfneg_file);
                sf_floatwrite(QGpos[jz][jx], np, qgpos_file);
                sf_floatwrite(QGneg[jz][jx], np, qgneg_file);
                sf_floatwrite(integrand_pos[jz][jx], np, integrandpos_file);
                sf_floatwrite(integrand_neg[jz][jx], np, integrandneg_file);
                fprintf(stderr, "YOOOO!!\n");
        /*
        split_normalize(f[0][0], f_pos, f_neg, t, nt);
        split_normalize(g[0][0], g_pos, g_neg, t, nt);

        sf_putint(cdffpos_file, "n1", nt);
        sf_putint(cdffpos_file, "n2", 1);
        sf_putint(cdffpos_file, "n3", 1);
            
        sf_putint(cdffneg_file, "n1", nt);
        sf_putint(cdffneg_file, "n2", 1);
        sf_putint(cdffneg_file, "n3", 1);

        sf_putint(cdfgpos_file, "n1", nt);
        sf_putint(cdfgpos_file, "n2", 1);
        sf_putint(cdfgpos_file, "n3", 1);

        sf_putint(cdfgneg_file, "n1", nt);
        sf_putint(cdfgneg_file, "n2", 1);
        sf_putint(cdfgneg_file, "n3", 1);

        cdf(f_pos, t, F_pos[0][0], nt);
        cdf(f_neg, t, F_neg[0][0], nt);

        cdf(g_pos, t, G_pos[0][0], nt);
        cdf(g_neg, t, G_neg[0][0], nt);

        quantile(p, t, F_pos[0][0], QFpos[0][0], np, nt);
        quantile(p, t, F_neg[0][0], QFneg[0][0], np, nt);
        quantile(p, t, G_pos[0][0], QGpos[0][0], np, nt);
        quantile(p, t, G_neg[0][0], QGneg[0][0], np, nt);

        wass_int(integrand_pos[0][0], f_pos, g_pos, t,p, nt, np);
        wass_int(integrand_neg[0][0], f_neg, g_neg, t, p, nt, np);

        //fprintf(stderr, "HERE!!!\n");
        sf_floatwrite(F_pos[0][0], nt, cdffpos_file);
        //fprintf(stderr, "HERE FIRST!!!\n");
        sf_floatwrite(F_neg[0][0], nt, cdffneg_file);
        sf_floatwrite(G_pos[0][0], nt, cdfgpos_file);
        sf_floatwrite(G_neg[0][0], nt, cdfgneg_file);
        //fprintf(stderr, "Cumulative dist. works\n");
        sf_floatwrite(f[0][0], nt, ff_file);
        //fprintf(stderr, "Density writing works\n");
        sf_floatwrite(f_pos, nt, fpos_file);
        sf_floatwrite(f_neg, nt, fneg_file);
        sf_floatwrite(g[0][0], nt, gg_file);
        sf_floatwrite(g_pos, nt, gpos_file);
        sf_floatwrite(g_neg, nt, gneg_file);
        sf_floatwrite(p, np, p_file);
        sf_floatwrite(QFpos[0][0], np, qfpos_file);
        sf_floatwrite(QFneg[0][0], np, qfneg_file);
        sf_floatwrite(QGpos[0][0], np, qgpos_file);
        sf_floatwrite(QGneg[0][0], np, qgneg_file);
        sf_floatwrite(integrand_pos[0][0], np, integrandpos_file);
        sf_floatwrite(integrand_neg[0][0], np, integrandneg_file);
        */
        }
    }
}
}

