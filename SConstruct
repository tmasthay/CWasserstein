from rsf.proj import *
from rtest import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from time import time
import m8r

test = False

if( not test ):
    t = time()
    
    nz = 20
    nx = 20
    z = np.linspace(0.2, 0.8, nz)
    x = np.linspace(0.2, 0.8, nx)
    
    pts = np.array([yy for yy in product(z,x)])
    
    misfits = run_test(pts)
    Z,X = np.meshgrid(z,x)
    misfits = misfits.reshape(Z.shape)
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(Z,X,misfits)
    plt.title('W2 Hyperbolic Radon')
    plt.savefig('W2-hyperbolictrace.png')
    
    print('Total time: %.4f'%(time() - t))
else:
    Flow('test-file',None,'math output="1.0" n1=100 d1=0.01')



