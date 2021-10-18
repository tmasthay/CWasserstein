#include <rsf.h>

float* linspace(float a, float b, int N){
    float *x, delta;
    int i;

    x = sf_floatalloc(N);
    x[0] = a;

    delta = (b-a) / (N-1);

    for(i = 1; i < N; i++){
        x[i] = x[i-1] +  delta;
    }
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
        fprintf(stderr, "ig = %d\n", ig);
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
    }
    if( sum == 0.0 ) return;
    float reciprocal = 1.0 / sum;
    for( i = 0; i < N; i++ ){
        f[i] = reciprocal * f[i];
    }
}

void split_normalize(float* f, 
    float* f_pos, 
    float* f_neg, 
    float* x,
    int N){
    int i;
    for(i = 0; i < N; i++){
        f_pos[i] = f[i] >= 0 ? f[i] : 0.0;
        f_neg[i] = f[i] < 0 ? -f[i] : 0.0;
    }
    renormalize(f_pos, x, N);
    renormalize(f_neg, x, N);
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
    for(ip = 0; ip < Np; ip++){
        //int found = 0;

        for(ix = sx; ix < Nx; ix++){
            if( F[ix-1] <= p[ip] && p[ip] <= F[ix] ){
                sx = ix;
                if( F[ix] == F[ix-1] ){
                    G[ip] = 0.5 * (x[ix] + x[ix-1]);
                    //fprintf(stderr, "HERE!\n");
                }
                else{
                    float t = (F[ix] - p[ip]) / (F[ix] - F[ix-1]);
                    G[ip] = t * x[ix-1] + (1-t) * x[ix];
                }
               // fprintf(stderr, "(ix, sx, p, G) = (%d, %d, %.8f, %.8f)\n", ix, sx, p[ip], G[ip]);
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

    quantile(p, x, F, Finv, P, N);
    quantile(p, x, G, Ginv, P, N);

    for(ip = 0; ip < P; ip++){
        dest[ip] = Finv[ip] - Ginv[ip];
        dest[ip] = dest[ip] * dest[ip];
    }
}

float wass2(float* f, float* g, float* x, float*p, int N, int P){
    int ix;
    float *y, sum;

    y = sf_floatalloc(P);

    //transport(f, g, T, x, N);
    wass_int(y, f, g, x, p, N, P);
    
    sum = 0.0;
    for(ix = 0; ix < N; ix++){
        sum = sum + 0.5 * (x[ix] - x[ix-1]) * (y[ix] + y[ix-1]);
    }
    //sum = integrate(y, linspace(0.0, 1.0, 10*N), 10*N);
    return sum;
}

float wass2_split(float* f, float* g, float* x, float* p, int N, int P){
    float *f_pos, *f_neg, *g_pos, *g_neg;
    float pos=0.0, neg=0.0;
    
    f_pos = sf_floatalloc(N);
    f_neg = sf_floatalloc(N);
    g_pos = sf_floatalloc(N);
    g_neg = sf_floatalloc(N);

    split_normalize(f, f_pos, f_neg, x, N);
    split_normalize(g, g_pos, g_neg, x, N);

    pos = wass2(f_pos, g_pos, x, p, N, P);
    neg = wass2(f_neg, g_neg, x, p, N, P);

    return pos*pos + neg*neg;
}

float wass2trace(float*** f, float*** g, float* t, int nz, int nx, int nt, int np){
    //assume n3xn2xn1=(z,x,t) coordinates
    float sum=0.0, *p;
    int iz, ix;

    p = sf_floatalloc(np);
    p = linspace(0.0, 1.0, np);

    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            sum += wass2_split(f[iz][ix], g[iz][ix], t, p, nt, np);
        }
    }
    return sum;
}

int main(int argc, char* argv[]){
    //assume structured grid
    int nz, nx, nt, nt_true, nz_check, nx_check, nt_check;
    float ***f, ***g, *t;

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
    if (!sf_histint(f_file, "n3", &nz)) sf_error("No n3");
    if (!sf_histint(f_file, "n2", &nx)) sf_error("No n2");
    if (!sf_histint(f_file, "n1", &nt)) sf_error("No n1");
    if (!sf_histint(g_file, "n3", &nz_check)) sf_error("No n3: g");
    if (!sf_histint(g_file, "n2", &nx_check)) sf_error("No n2: g");
    if (!sf_histint(g_file, "n1", &nt_check)) sf_error("No n1: g");
    if (!sf_histint(t_file, "n1", &nt_true)) sf_error("No n1: t");

    fprintf(stderr, "nz = %d", nz);
    fprintf(stderr, "nx = %d", nx);
    fprintf(stderr, "nt = %d", nt);

    //Input validation
    if( nz_check != nz ) sf_error("Dimension mismatch: z");
    if( nx_check != nx ) sf_error("Dimension mismatch: x");
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
            //fprintf(stderr, "(iz,ix) = (%d, %d)\n", iz, ix);
            sf_floatread(f[iz][ix], nt, f_file);
            sf_floatread(g[iz][ix], nt, g_file);
        }
    }
    sf_floatread(t, nt, t_file);

    //transport(f,g,T,x,N);
    float distance;
    int np=5*nt;
    distance = wass2trace(f,g,t,nz, nx, nt,np);

    float* distance_tmp;
    distance_tmp = sf_floatalloc(1);
    int it;
    for(it = 0; it < 1; it++){
        distance_tmp[it] = distance;
    }
    sf_putint(wass_file, "n1", 1);
    sf_putint(wass_file, "n2", 1);
    sf_putint(wass_file, "n3", 1);
    sf_floatwrite(distance_tmp, 1, wass_file);
}
