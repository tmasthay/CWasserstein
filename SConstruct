from rsf.proj import *
from rtest import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from time import time
import sys

landscape = True
inversion = False

if( landscape ):
    t = time()
    
    nz = 25
    nx = 25
    z = np.linspace(0.2, 0.8, nz)
    x = np.linspace(0.2, 0.8, nx)
    
    pts = np.array([yy for yy in product(z,x)])
    
    set_mode(1 if len(sys.argv) == 1 else int(sys.argv[1]))
    misfits = run_test(pts)
    Z,X = np.meshgrid(z,x)
    misfits = misfits.reshape(Z.shape)
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(Z,X,misfits)
    plt.title('W2 Hyperbolic Radon')
    plt.savefig('W2-hyperbolictrace.png')
    
    print('Total time: %.4f'%(time() - t))

if( inversion ):
   os.system('python inversion.py')
    


