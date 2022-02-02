import sys
import numpy as np
from unwrap import unwrap
import copy
from seqflow import *
import matplotlib.pyplot as plt
from subprocess import check_output as co

Flow=SeqFlow
Plot=SeqPlot

def create_r(p):
    '''
    Params p: 
    srcz -- source z coordinate
    srcx -- source x coordinate

    oz,dz,nz -- origin z, delta_z, number of z 
    ox,dx,nx -- origin x, delta_x, number of x
    '''
    sz = p['sz']
    sx = p['sx']
    oz = p['oz']
    dz = p['dz']
    nz = p['nz']
    ox = p['ox']
    dx = p['dx']
    nx = p['nx']

    y1 = "(x1-%.10e)*(x1-%.10e)"%(sz, sz)
    y2 = "(x2-%.10e)*(x2-%.10e)"%(sx,sx)
    Flow("r_%.4f_%.4f"%(sz,sx), None, 
        '''
        math output="sqrt(%s + %s)"
            o1=%.10e d1=%.10e n1=%d
            o2=%.10e d2=%.10e n2=%d
        '''%(y1,y2,oz,dz,nz,ox,dx,nx))

def get_ot(f1, f2, suffix, ot_mode):
    Flow('tmp_transp', f1, 'transp plane=13')
    Flow('tmp_transp2', f2, 'transp plane=13')
    Flow('diff_%s'%suffix, 'tmp_transp tmp_transp2',
        '''
        ./ot.exe g=${SOURCES[1]} t=t.rsf mode=%d
        '''%ot_mode)

    res = co('sfdisfil < diff_%s.rsf'%suffix, shell=True).decode('utf-8')
    return float(res.split(':')[-1].replace(' ',''))


def create_cases(p):
    '''
    Params p:
    k -- characteristic wavenumber
    
    srcz -- source z coordinate
    srcx -- source x coordinate

    oz,dz,nz -- origin z, delta_z, number of z 
    ox,dx,nx -- origin x, delta_x, number of x
    '''
    k = p['k']
    srcz = p['srcz']
    srcx = p['srcx']
    ref = p['ref']
    ref_name = p['ref_name']
    ot_mode = p['ot_mode']
    ox = p['ox']
    dx = p['dx']
    nx = p['nx']

    dist = np.zeros((len(srcz), len(srcx)))

    for (iz,sz) in enumerate(p['srcz']):
        for (ix,sx) in enumerate(p['srcx']):
            d = copy.copy(p)
            d['sz'] = sz
            d['sx'] = sx
            create_r(d)

            suffix = '%.4f_%.4f'%(sz, sx)
            Flow('green_real_%s'%suffix,
                'r_%s'%suffix,
                '''
                math output="cos(%.10e * input) / input" 
                '''%k)
            Flow('green_imag_%s'%suffix,
                'r_%s'%suffix,
                '''
                math output="sin(%.10e * input) / input"
                '''%k)

            Flow('green_real_top_%s'%suffix,
                 'green_real_%s'%suffix,
                 'window n1=1 f1=0')

            Flow('green_imag_top_%s'%suffix,
                 'green_imag_%s'%suffix,
                 'window n1=1 f1=0')

            Plot('green_real_top_%s'%suffix,None,'graph')
            Plot('green_imag_top_%s'%suffix,None,'graph')

            if( not ref ):
                dist1 = get_ot('green_real_top_%s'%suffix,
                    'green_real_top_%s'%ref_name,
                    'real_' + suffix,
                    ot_mode)

                dist2 = get_ot('green_imag_top_%s'%suffix,
                    'green_imag_top_%s'%ref_name,
                    'imag_' + suffix,
                    ot_mode)
                dist[iz][ix] = dist1 + dist2
                
            else:
                Flow('t',None,
                    '''
                    math output="x1" o1=%.10e d1=%.10e n1=%d
                    '''%(ox,dx,nx))
    return dist
               
def go():
    nz = 50
    nx = 50
    d = {
        'k':2.0,
        'srcz': np.linspace(1.0, 10.0, nz),
        'srcx': np.linspace(-10.0, 10.0, nx),
        'oz': 0.0,
        'ox': -10.0,
        'nz':200,
        'nx': 200,
        'dz':0.1,
        'dx':0.1,
        'ref': False,
        'ot_mode': int(sys.argv[1])
    }

    sz_ref = d['srcz'][int(nz / 2)]
    sx_ref = d['srcx'][int(nx / 2)]

    d_ref = copy.copy(d)
    d_ref['srcz'] = [sz_ref]
    d_ref['srcx'] = [sx_ref]
    d_ref['ref'] = True
    d_ref['ref_name'] = '%.4f_%.4f'%(sz_ref, sx_ref)

    d['ref_name'] = d_ref['ref_name']

    ignore_val = create_cases(d_ref)
    dist = create_cases(d)

    np.save('dist-%d.npy'%d['ot_mode'], dist)
    np.save('z-%d.npy'%d['ot_mode'], d['srcz'])
    np.save('x-%d.npy'%d['ot_mode'], d['srcx'])

    Z,X = np.meshgrid(d['srcz'],d['srcx'])
    fig,ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(Z,X,dist)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig('ot_convexity-%d.png'%d['ot_mode'])

go()
 
