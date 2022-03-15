def create_base():
    N = 100
    
    nz = N
    dz = 1.0 / nz
    
    nx = N
    dx = 1.0 / nx
    
    nt = 1000
    dt = 5e-03
    
    d_forward = {
            'case' : 'synthetic',
            'nz' : nz,
            'dz' : dz,
            'nx' : nx,
            'dx' : dx,
            'dt' :  dt,
            'nt' :  nt,
            'sszf' :  0.5,
            'ssxf' :  0.5,
            'eszf' :  0.5,
            'esxf' :  0.5,
            'nsz' :  1,
            'nsx' :  1,
            'nb' : 35,
            'fm' : 5.0,
            'amp': 1000.0,
            'vp': '0.5 + 10*sin(x1+x2)',
            'vs': '0.707 * input',
            'rho': '1.0'
    }
    return d_forward
