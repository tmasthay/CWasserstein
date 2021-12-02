from rsf.proj import *
from rtest import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from time import time
import sys
from subprocess import check_output as co
from mode_to_str import *

landscape = True
inversion = False

def setup_output_directory(case_name):
    fig_dir = '../figures'
    date_info = co('date',shell=True) \
        .decode('utf-8') \
        .replace('\n','') \
        .split(' ')
    
    date_form = ''.join(date_info[1:4]) + '_' + date_info[-1]
    top_dir = fig_dir + '/' + date_form

    if( not os.path.exists( top_dir ) ):
        os.system('mkdir %s'%top_dir)

    top_dir = top_dir + '/' + '_'.join(case_name.split(' '))
    if( not os.path.exists( top_dir ) ):
        os.system('mkdir %s'%top_dir)

    top_dir = top_dir + '/' + date_info[4]
    if( os.path.exists( top_dir ) ):
        print('FATAL ERROR: two jobs sent off within second')
        print('CHECK TIME ZONE SETTINGS ON COMPUTER')
        exit(-1)
    else:
        os.system('mkdir %s'%top_dir)

    f = open(top_dir + '/GITHASH.txt', 'w')
    f.write(co('git log', shell=True) \
        .decode('utf-8') \
        .split('commit')[1] \
        .split()[0])

    return top_dir

def run_mode(mode):
    if( landscape ):
        t = time()
        
        nz = 10
        nx = 10
        z = np.linspace(0.2, 0.8, nz)
        x = np.linspace(0.2, 0.8, nx)
        
        pts = np.array([yy for yy in product(z,x)])
        
        set_mode(mode)

        misfits = run_test(pts)
        Z,X = np.meshgrid(z,x)
        misfits = misfits.reshape(Z.shape)
        
        plt_title = mode_to_str(mode)
        dir = setup_output_directory(plt_title)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(Z,X,misfits)
        plt.title(plt_title)
        plt.savefig(dir + '/' + '_'.join(plt_title.split(' ')) + '.png')
        
        print('Total time: %.4f'%(time() - t))
    
    if( inversion ):
       os.system('python inversion.py')

def go():
    if( len(sys.argv) == 1 ):
        run_mode(1)
    else:
        y = [int(xx) for xx in sys.argv.split(' ')[1:]]
        for yy in y:
            run_mode(yy)
    


