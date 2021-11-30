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

float integrate_quad(float* f, float* x, int N){
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
        ////////fprintf(stderr, "ig = %d\n", ig);
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

float integrate(float* f, float* x, int N){
    if( N < 2 ){
        //fprintf(stderr, "Need more than 1 point\n");
	exit(-1);
    }
    int ix;
    //assume regular grid
    float dx = x[1] - x[0];
    float sum = 0.0;
    for(ix = 1; ix < N; ix++){
	sum = sum + 0.5 * (f[ix] - f[ix-1]);
    }
    return sum;
}

float integrate2(float** f, float* x, float* y, int M, int N){
    //fprintf(stderr, "integrate: entering\n");
    float sum=0.0;
    int i;

    if( N < 2 || M < 2 ){
        //fprintf(stderr, "Need more points M=%d, N=%d\n", M,N);
	exit(-1);
    }
    float dx = x[1] - x[0];
    //fprintf(stderr, "integrate: dx evaluated...entering loop\n");
    for(i = 0; i < M; i++){
        sum = sum + 0.5 * dx * integrate(f[i], y, N);
    }
    return sum;
}

void renormalize2(float** f, float* x, float* y, int M, int N){
    //fprintf(stderr, "Entering renormalize\n");
    int i,j;
    float integral = integrate2(f,x,y,M,N);

    float tol = 1e-8;

    //fprintf(stderr, "renormalize: integral evaluated\n");
    if( abs(integral) > tol ){
        for(i = 0; i < M; i++){
           for(j = 0; j < N; j++){
	       f[i][j] = f[i][j] / integral;
	   }
        }
    }
    //fprintf(stderr, "renormalize: returning\n");
}

void split_normalize2(float** f,
    float** f_pos,
    float** f_neg,
    float* x,
    float* y,
    int M,
    int N){
    //fprintf(stderr, "Entering split_normalize2.\n");
    int ix,iy;
    for(ix = 0; ix < M; ix++){
        for(iy = 0; iy < N; iy++){
	    if( f[ix][iy] > 0 ){
	        f_pos[ix][iy] = f[ix][iy];
		f_neg[ix][iy] = 0.0;
	    }
	    else{
	        f_pos[ix][iy] = 0.0;
		f_neg[ix][iy] = abs(f[ix][iy]);
	    }
	}
    }

    //fprintf(stderr, "split_normalize Loop success...renormalizing\n");
    renormalize2(f_pos, x, y, M, N);
    renormalize2(f_neg, x, y, M, N);
    //fprintf(stderr, "split_normalize success\n");
}

float* moment1(float** f, float* x, float* y, int M, int N){
    //fprintf(stderr, "entering moment1\n");

    float* answer;
    answer = sf_floatalloc(2);

    answer[0] = 0.0;
    answer[1] = 0.0;

    int ix, iy;
    float dx=x[1]-x[0];
    float dy=y[1]-y[0];

    //fprintf(stderr, "Moment1 initialization success\n");

    for(ix = 0; ix < M; ix++){
        for(iy = 0; iy < N; iy++){
	    //fprintf(stderr, "moment1 (ix,iy) = (%d,%d)\n", ix, iy);
	    answer[0] += 0.25 * dx * dy * x[ix] * f[ix][iy];
	    answer[1] += 0.25 * dx * dy * y[iy] * f[ix][iy];
	}
    }
    //fprintf(stderr, "moment1 success...returning.\n")
    return answer;
}

float moment2(float** f, float* x, float* y, int M, int N){
    float answer=0.0;

    int ix,iy;
    float dx = x[1] - x[0];
    float dy = y[1] - y[0];
    for(ix = 0; ix < M; ix++){
	float xsq = x[ix] * x[ix];
        for(iy = 0; iy < N; iy++){
	    float euclid2 = xsq + y[iy] * y[iy];
	    float increment = 0.25 * dx * dy * euclid2 * f[ix][iy];
	    if( increment < 0 ){
	        fprintf(stderr, "(inc, dx, dy, euclide2, f[%d][%d]) = (%.15f, %.15f, %.15f, %.15f, %.15f)\n", ix, iy, increment, dx, dy, euclid2, f[ix][iy]);
	    }
	    answer += 0.25 * dx * dy * euclid2 * f[ix][iy];
	}
    }
    return answer;
}

float sw2(float** f, float** g, float* x, float* y, int M, int N){
    //fprintf(stderr, "Entering SW2\n");

    float **f_pos, **f_neg, **g_pos, **g_neg;
    float *fpm1, *fnm1, *gpm1, *gnm1;
    float fpm2, fnm2, gpm2, gnm2;

    f_pos = sf_floatalloc2(M,N);
    f_neg = sf_floatalloc2(M,N);
    g_pos = sf_floatalloc2(M,N);
    g_neg = sf_floatalloc2(M,N);

    fpm1 = sf_floatalloc(2);
    fnm1 = sf_floatalloc(2);
    gpm1 = sf_floatalloc(2);
    gnm1 = sf_floatalloc(2);

    //fprintf(stderr, "SW2 allocations passed.\n");

    split_normalize2(f, f_pos, f_neg, x, y, M, N);
    split_normalize2(g, g_pos, g_neg, x, y, M, N);

    //fprintf(stderr, "SW2 renormalization success\n");

    fpm1 = moment1(f_pos, x, y, M, N);
    fnm1 = moment1(f_neg, x, y, M, N);
    gpm1 = moment1(g_pos, x, y, M, N);
    gnm1 = moment1(g_neg, x, y, M, N);

    //fprintf(stderr, "SW2 first moments evaluated\n");

    fpm2 = moment2(f_pos, x, y, M, N) 
	- pow(fpm1[0], 2.0) 
	- pow(fpm1[1], 2.0);
    fnm2 = moment2(f_neg, x, y, M, N) 
	- pow(fnm1[0], 2.0) 
	- pow(fnm1[1], 2.0);
    gpm2 = moment2(g_pos, x, y, M, N) 
	- pow(gpm1[0], 2.0) 
	- pow(gpm1[1], 2.0);
    gnm2 = moment2(g_neg, x, y, M, N) 
	- pow(gnm1[0], 2.0) 
	- pow(gnm1[1], 2.0);

    //fprintf(stderr, "SW2 second moments evaluated\n");

    float mom1_term = pow(fpm1[0] - gpm1[0], 2.0) +
        pow(fpm1[1] - gpm1[1], 2.0) + 
	pow(fnm1[0] - gnm1[0], 2.0) +
	pow(fnm1[1] - gnm1[1], 2.0);

    float mom2_term = pow(sqrt(fpm2) - sqrt(gpm2), 2.0) +
        pow(sqrt(fnm2) - sqrt(gnm2), 2.0);

    //fprintf(stderr, "SW2 final terms computed\n");
    

    float x_tmp = -0.00000000000;
    return 0.5 * (mom1_term + mom2_term);
}

int main(int argc, char* argv[]){
    //fprintf(stderr, "Main hit\n");

    //assume structured grid
    int nz, nx, nt, nt_true, nz_check, nx_check, nt_check;
    float ***f, ***g, *x, *t;

    sf_file f_file, g_file, x_file, t_file, wass_file;

    sf_init(argc, argv);

    //define inputs
    f_file = sf_input("in");
    g_file = sf_input("g");
    x_file = sf_input("x");
    t_file = sf_input("t");

    //define outputs
    wass_file = sf_output("out");


    //read command-line inputs
    if (!sf_histint(f_file, "n1", &nz)) sf_error("No n3");
    if (!sf_histint(f_file, "n2", &nx)) sf_error("No n2");
    if (!sf_histint(f_file, "n3", &nt)) sf_error("No n1");
    if (!sf_histint(g_file, "n1", &nz_check)) sf_error("No n3: g");
    if (!sf_histint(g_file, "n2", &nx_check)) sf_error("No n2: g");
    if (!sf_histint(g_file, "n3", &nt_check)) sf_error("No n1: g");
    if (!sf_histint(t_file, "n1", &nt_true)) sf_error("No n1: t");

    //Input validation
    if( nz_check != nz ) sf_error("Dimension mismatch: z");
    if( nx_check != nx ) sf_error("Dimension mismatch: x");
    if( nt_check != nt ) sf_error("Dimension mismatch: t1 yo");
    if( nt != nt_true ) sf_error("Dimension mismatch: t2 yo, nt=%d, ntrue=%d", nt, nt_true);

    //allocate memory
    f = sf_floatalloc3(nz, nx, nt);
    g = sf_floatalloc3(nz, nx, nt);
    x = sf_floatalloc(nx);
    t = sf_floatalloc(nt);

    //read in data 
    int iz, ix;
    for(iz = 0; iz < nz; iz++){
        for(ix = 0; ix < nx; ix++){
            sf_floatread(f[iz][ix], nt, f_file);
            sf_floatread(g[iz][ix], nt, g_file);
        }
    }
    sf_floatread(x, nx, x_file);
    sf_floatread(t, nt, t_file);

    
    //fprintf(stderr, "Data read in properly.\n");

    float distance = sw2(f[0], g[0], x, t, nx, nt);

    //fprintf(stderr, "SW distance successfully approximated.\n");

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
}

