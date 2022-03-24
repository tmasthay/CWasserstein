from rsf.proj import * 
from rtest import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from time import time
import sys
from subprocess import check_output as co
from mode_to_str import *
import signal

global TOP_DIR

#Usage
# python SConstruct modes=[] grid=[nz,nx] threshold=0
#ignore me...ssh git test

landscape = True
inversion = False

def parse(s):
    is_present = eval(' or '.join([str('%s'%s in t) for t in sys.argv])) 
    if( is_present ):
        return eval(':'.join(sys.argv).split('%s='%s)[1].split(':')[0])
    else:
        return None

def setup_output_directory(case_name):
    global TOP_DIR
    fig_dir = '../figures'
    date_info = co('date',shell=True) \
        .decode('utf-8') \
        .replace('\n','') \
        .split(' ')
    
    date_form = ''.join(date_info[1:3]) + '_' + date_info[-1]
    top_dir = fig_dir + '/' + date_form

    if( not os.path.exists( top_dir ) ):
        os.system('mkdir %s'%top_dir)

    top_dir = top_dir + '/' + date_info[3]
    if( os.path.exists( top_dir ) ):
        print('FATAL ERROR: two jobs sent off within second')
        print('CHECK TIME ZONE SETTINGS ON COMPUTER')
        exit(-1)
    else:
        os.system('mkdir %s'%top_dir)

    top_dir = top_dir + '/' + '_'.join(case_name.split(' '))
    if( not os.path.exists( top_dir ) ):
        os.system('mkdir %s'%top_dir)

    f = open(top_dir + '/GITHASH.txt', 'w')
    f.write(co('git log', shell=True) \
        .decode('utf-8') \
        .split('commit')[1] \
        .split()[0])

    if( not os.path.exists( top_dir + '/vpl_files' ) ):
        os.system('mkdir %s'%(top_dir + '/vpl_files'))

    if( not os.path.exists( top_dir + '/rsf_files' ) ):
        os.system('mkdir %s'%(top_dir + '/rsf_files'))

    return top_dir

def move_intermediate_files(top_dir):
    print('MOVING INTERMEDIATE RSF AND VPL FILES')
    os.system('mv *.vpl %s/vpl_files'%top_dir)
    os.system('mv *.rsf %s/rsf_files'%top_dir)

def handler(signum, frame):
    global TOP_DIR
    input('CTRL+C PRESSED: EXITING AFTER MOVING FILES')
    move_intermediate_files(TOP_DIR)
    exit(1)
    
def run_mode(mode):
    if( landscape ):
        global TOP_DIR
        signal.signal(signal.SIGINT, handler)
        t = time()

        print('yo1')        
        zx_dims = parse('grid')
        if( type(zx_dims) == type(None) ):
            nz = 10
            nx = 10
        else:
            nz = zx_dims[0]
            nx = zx_dims[1]

        z = np.linspace(0.2, 0.8, nz)
        x = np.linspace(0.2, 0.8, nx)
        
        pts = np.array([yy for yy in product(z,x)])
        
        set_mode(mode)

        misfits = run_test(pts)
        threshold = parse('threshold')
        if( type(threshold) != type(None) ):
            for i in range(len(misfits)):
                misfits[i] = min(threshold, misfits[i])

        Z,X = np.meshgrid(z,x)
        misfits = misfits.reshape(Z.shape)

        print('Z:\n%s'%Z)
        print('X:\n%s'%X)
        print('misfits:\n%s'%misfits)

        plt_title = mode_to_str(mode)
        my_dir = setup_output_directory(plt_title)
        fig,ax = plt.subplots()
        img = plt.imshow(misfits, \
            extent=[np.min(Z), np.max(Z), np.min(X), np.max(X)])
        fig.colorbar(img, ax=ax)
        plt.xlabel('Z')
        plt.ylabel('X')
        plt.title(plt_title)
        case_dir = my_dir + '/' + '_'.join(plt_title.split(' '))
        plt.savefig(case_dir + '.png')
        np.save(case_dir + 'Z.npy', Z)
        np.save(case_dir + 'X.npy', X)
        np.save(case_dir + 'misfits.npy', misfits)


        os.system('cp %s.png .'%(case_dir))
        os.system('mv %s.png 1.png'%(case_dir.split('/')[-1]))
        print(case_dir)
        #move_intermediate_files(TOP_DIR)
        
        print('Total time: %.4f'%(time() - t))
    
    if( inversion ):
       os.system('python inversion.py')

def go():
    modes = parse('modes')
    if( type(modes) != type(None) ):
        for m in modes:
            run_mode(m)
    else:
        run_mode(1)

signal.signal(signal.SIGINT, handler)
signal.signal(signal.SIGALRM, handler)

go()
