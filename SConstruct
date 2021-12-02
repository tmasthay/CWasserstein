from rsf.proj import *
from rtest import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from time import time
import sys

landscape = True
inversion = False

def mode_to_str(mode):
    if( mode == 1 ):
        return 'W2 trace split-renormalize surface'
    elif( mode == 2 ):
        return 'L2 surface'
    elif( mode == 3 ):
        return 'W2 trace split-renormalize global'
    elif( mode == 4 ):
        return 'L2 global'
    elif( mode == 5 ):
        return 'W2 trace abs-renormalize surface'

def go():
    if( landscape ):
        t = time()
        
        nz = 25
        nx = 25
        z = np.linspace(0.2, 0.8, nz)
        x = np.linspace(0.2, 0.8, nx)
        
        pts = np.array([yy for yy in product(z,x)])
        

        mode = 1 if len(sys.argv) == 1 else int(sys.argv[1])
        set_mode(mode)

        misfits = run_test(pts)
        Z,X = np.meshgrid(z,x)
        misfits = misfits.reshape(Z.shape)
        
        plt_title = mode_to_str(mode)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(Z,X,misfits)
        plt.title(plt_title)
        plt.savefig('_'.join(plt_title.split(' ')) + '.png')
        
        print('Total time: %.4f'%(time() - t))
    
    if( inversion ):
       os.system('python inversion.py')

go()
    


