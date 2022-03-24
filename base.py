import numpy as np
import m8r
from seqflow import *

Flow=SeqFlowV

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

    def create_vp(the_name):
        layer_vals = [5, 10, 15]
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

    def create_vs(the_name):
        layer_vals = [5, 10, 15]
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

    d_forward['vp'] = create_vp
    d_forward['vs'] = create_vs
    
    return d_forward
