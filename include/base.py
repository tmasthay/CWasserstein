import numpy as np
import m8r
import copy
from rsf.proj import *
from helper import grey

def layered_medium(layer_vals, layer_ends, x):
    return np.array([layer_vals[i] for i in np.digitize(x, layer_ends)])

def run_once(f):
    def wrapper(s):
        if wrapper.num_calls < 3 :
            wrapper.num_calls += 1  
            return f(s)
    wrapper.num_calls=0
    return wrapper

@run_once
def create_vp(args):
    #get dictionary and target file name
    d = args['dict']
    the_name = d[args['key']]

    #get dof number
    nz = d['nz']
    nx = d['nx']

    #grid step
    dz = d['dz']
    dx = d['dx']

    #grid origins
    oz = d['oz']
    ox = d['ox']

    #get sigmoid tuning
    alpha = args['alpha']

    #change velocities for cfl condition testing
    cfl = 1.0

    #define the layers
    layer_vals = cfl * np.array([6.0, 8.0  , 8.5  , 9.3  , 9.8  , 10.0 , 11.0, 11.8])
    layer_ends = [5.0, 250.0, 450.0, 500.0, 520.0, 575.0, 625.0]

#    layer_vals = [2.0]
#    layer_ends = []
    
    cmd = 'sfmath output="%.8f'%(layer_vals[0])
    for i in range(len(layer_ends)):
        cmd += ' + %.8f / (1.0 + exp(-%.8f * (x1 - %.8f)))'%(
            layer_vals[i+1], 
            alpha, 
            layer_ends[i])
    cmd += '" o1=%.8f o2=%.8f n1=%d n2=%d d1=%.8f d2=%.8f'%(
        oz, ox, nz, nx, dz, dx)

    Flow(the_name, None, cmd)

    plot_cmd = grey('P velocity', color='tokyo')
    Result(the_name, the_name, plot_cmd)
        

    return None

@run_once
def create_vs(args):
    #get dictionary and target file name
    d = args['dict']
    the_name = d[args['key']]

    #get dof number
    nz = d['nz']
    nx = d['nx']

    #grid step
    dz = d['dz']
    dx = d['dx']

    #grid origins
    oz = d['oz']
    ox = d['ox']

    #get sigmoid tuning
    alpha = args['alpha']

    #just have vs=sqrt(2) vp for now
    cfl = 1.0

    #define the layers
    layer_vals = cfl * np.array([3.75, 4.3, 4.6, 5.0, 5.25, 5.5, 6.0, 6.15, 6.3])
    layer_ends = [5.0, 250.0, 450.0, 500.0, 520.0, 575.0, 625.0]
   
    #layer_vals = [2.0*0.707]
    #layer_ends = []
    cmd = 'sfmath output="%.8f'%(layer_vals[0])
    for i in range(len(layer_ends)):
        cmd += ' + %.8f / (1.0 + exp(-%.8f * (x1 - %.8f)))'%(
            layer_vals[i+1], 
            alpha, 
            layer_ends[i])
    cmd += '" o1=%.8f o2=%.8f n1=%d n2=%d d1=%.8f d2=%.8f'%(
        oz, ox, nz, nx, dz, dx)

    Flow(the_name, None, cmd)

    plot_cmd = grey('S velocity', color='tokyo')
    Result(the_name, the_name, plot_cmd)

    return None

@run_once
def create_rho(args):
    d = args['dict']
    the_name = d[args['key']]
    nz = d['nz']
    nx = d['nx']
    dz = d['dz']
    dx = d['dx']
    rho = d['rho_expr']

    cmd = 'math output="%s" d1=%.f n1=%d d2=%.4f n2=%d'%(
            rho,
            dz, nz, 
            dx, nx)
    
    Flow(the_name, None, cmd)

    plot_cmd = grey('Density', color='tokyo')
    Result(the_name, the_name, plot_cmd)

    return the_name, None, cmd

def create_base():
    N = 200
  
    
    az = 0.0
    bz = 700.0 
    nz = N
    dz = (bz-az) / nz
   
    space_between_receivers=70.0
    num_receivers = 10.0 
    ax = 0.0
    dx = space_between_receivers / num_receivers
    nx = N
    bx = dx * nx
    
    nt = 10000
    dt = 0.0071
    
    '''
    #test for scale recovery
    az = 0.0
    bz = 1.0
    nz = N
    dz = 1.0 / N

    ax = 0.0
    bx = 1.0
    nx = N
    dx = 1.0 / N

    nt = 1000
    dt = 0.001
    '''

    amp = 1.0e3
    snr = 2.0
    
    d_forward = {
            'case' : 'synthetic',
            'oz' : az,
            'nz' : nz,
            'dz' : dz,
            'ox' : ax,
            'nx' : nx,
            'dx' : dx,
            'dt' :  dt,
            'nt' :  nt,
            'sszf' :  0.25,
            'ssxf' :  0.25,
            'eszf' :  0.75,
            'esxf' :  0.75,
            'nsz' :  2,
            'nsx' :  2,
            'nb' : 35,
            'fm' : 25.0,
            'amp': amp,
            'noise': amp / snr,
            'vp_expr': '0.5 + 10*sin(x1+x2)',
            'vs_expr': '0.707 * input',
            'rho_expr': '1.0'
    }

    d_forward['vp'] = 'vp'
    d_forward['vs'] = 'vs'
    d_forward['rho'] = 'rho'

    vp_args = {'key': 'vp', 'dict': d_forward, 'alpha': 100}

    vs_args = copy.copy(vp_args)
    vs_args.update({'key': 'vs'})

    rho_args = copy.copy(vp_args)
    rho_args.update({'key': 'rho'})

    create_vp(vp_args)
    create_vs(vs_args)
    create_rho(rho_args)

    return d_forward
