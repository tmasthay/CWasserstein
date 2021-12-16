import os
from seqflow import *

def get_dict(name, sz, sx):  
    nz=100
    nx=100
    nt=1000
    dz=1.0 / nz
    dx=1.0 / nx
    dt=5.0 / nt
    
    nb=30
    fm=5.0
    vp='0.5'
    vs='0.707 * input'
    rho='1.0'
    
    if( not os.path.exists('t.rsf') ):
         SeqFlowV('t.rsf', None, 
             '''
             math output="x1" d1=%.4f n1=%d 
             '''%(dt, nt))

         SeqFlowV('x.rsf', None, 
             '''
             math output="x1" d1=%.4f n1=%d 
             '''%(dx, nx))
    
    d1 = {
            'case' : name,
            'nz' : nz,
            'dz' : dz,
            'nx' : nx,
            'dx' : dx,
            'dt' :  dt,
            'nt' :  nt,
            'sszf' :  sz,
            'ssxf' :  sx,
            'eszf' :  sz,
            'esxf' :  sx,
            'nsz' :  1,
            'nsx' :  1,
            'nb' : nb,
            'fm' : fm,
            'vp': vp,
            'vs': vs,
            'rho': rho
    }
    
    return d1
 
    
