import numpy as np
import m8r

def layered_medium(layer_vals, layer_ends, x):
    return np.array([layer_vals[i] for i in np.digitize(x, layer_ends)])

def create_base():
    N = 100
    
    nz = N
    dz = 1.0 / nz
    
    nx = N
    dx = 1.0 / nx
    
    nt = 1000
    dt = 5e-03 

    amp = 1000.0
    snr = 2.0
    
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
            'amp': amp,
            'noise': amp / snr,
            'vp_expr': '0.5 + 10*sin(x1+x2)',
            'vs_expr': '0.707 * input',
            'rho_expr': '1.0'
    }

    @run_once
    def create_vp(the_name):
        layer_vals = [0.25, 0.5, 0.75]
        layer_ends = [0.25, 0.75]
        
        x = np.linspace(0,1,nx)

        layers = layered_medium(layer_vals, layer_ends, x)
        layers_full = np.array([layers for i in range(nz)])
        f = m8r.Output(the_name + '.rsf')

        f.put('d1', dz)
        f.put('d2', dx)
        f.put('n1', nz)
        f.put('n2', nx)

        f.write(layers_full)

    @run_once
    def create_vs(the_name):
        layer_vals = [0.25, 0.5, 0.75]
        layer_ends = [0.25, 0.75]
        
        x = np.linspace(0,1,nx)

        s_factor = 0.707
        layers = s_factor * layered_medium(layer_vals, layer_ends, x)
        layers_full = np.array([layers for i in range(nz)])
        f = m8r.Output(the_name + '.rsf')
        f.put('d1', dz)
        f.put('d2', dx)
        f.put('n1', nz)
        f.put('n2', nx)
        f.write(layers_full)

    @run_once
    def create_rho(the_name):
        Flow(the_name, None, 
            'math output=%s d1=%.f n1=%d d2=%.4f n2=%d'%(
                d_forward['rho_expr'],
                d_forward['dz'], d_forward['nz'], 
                d_forward['dx'], d_forward['nx']))

    d_forward['vp'] = 'vp'
    d_forward['vs'] = 'vs'
    d_forward['rho'] = 'rho'

    create_vp(d_forward['vp'])
    create_vs(d_forward['vs'])
    create_rho(d_forward['rho'])

    return d_forward
