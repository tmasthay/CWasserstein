from seqflow import *
from newtest import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from time import time
import sys
from subprocess import check_output as co
from mode_to_str import *
import signal
from purge import purge

Flow=SeqFlowV
Plot=SeqPlot

global TOP_DIR
global l2_mode

l2_mode = True

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

def sfl2(src):
    s = co('sfattr < %s.rsf'%src, shell=True) \
        .decode('utf-8') \
        .split('\n')[3] \
        .split('=')[-1]
    try:
        return float(s)
    except: 
        print('Error at %s.rsf'%src)
        exit(-1)
        

def square_normalize(dest, src, alpha=0.0):
    C = 1.0 / sfl2(src)
    Flow(dest, src, 'math output="%.16f*input*input + %.4f"'%(C,alpha))

def hs(f1, f2, n1, n2, del_char='z', alpha=0.0, pad=1, s=-1.0):
    g1 = f1 + '_norm'
    g2 = f2 + '_norm'
    square_normalize(g1, f1, alpha)
    square_normalize(g2, f2, alpha)

    h1 = f1 + '_fft'
    h2 = f2 + '_fft'
    Flow(h1, g1, 
        '''
        dd type=complex | fft3 axis=1 | fft3
        ''')
    Flow(h2, g2, 
        '''
        dd type=complex | fft3 axis=1 | fft3
        ''')
    
    Flow('tmp', '%s %s'%(h1,h2),
       '''
       math output="abs(input - y)"
           y=${SOURCES[1]} | dd type=float
       ''')
    if( not os.path.exists('kernel.rsf') ):
        Flow('kernel', 'tmp',
            '''
            math output="1.0 / (1.0 + x1^2 + x2^2)"
            ''')
    Flow('tmp_integrand', 'tmp kernel',
        '''
        math output="input * y" y=${SOURCES[1]}
        ''')
    final_result = sfl2('tmp_integrand')
    #os.system('sfrm wav%s_[0-9].*rsf'%del_char)
    #purge()
    return final_result

def l2diff(f1, f2, del_char='z'):
    Flow('diff', '%s %s'%(f1,f2),
        '''
        math output="input - y" y=${SOURCES[1]}
        ''')
    return_val = sfl2('diff')
    #os.system('sfrm wav%s_[0-9].*[0-9].rsf'%del_char)
    #purge()
    return return_val
    
 
def run_mode(mode):
    if( landscape ):
        global TOP_DIR
        global l2_mode
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
        
        misfits = np.zeros(len(pts))
        alpha = 0.0
        for (i,pt) in enumerate(pts):
            zz = pt[0]
            xx = pt[1]
            ztop, xtop = create_case(zz,xx)
            suff = d_forward['case'] + '_top'
            if( l2_mode ):
                misfits[i] = l2diff(ztop, 'wavz_' + suff, 'z')
                misfits[i] += l2diff(ztop, 'wavx_' + suff, 'x')
            else:
                misfits[i] = hs(ztop, 'wavz_' + suff, nz, nx, 'z', 
                    alpha=alpha)
                misfits[i] += hs(xtop, 'wavx_' + suff, nz, nx, 'x', 
                    alpha=alpha)
             
        min_threshold = parse('min_threshold')
        max_threshold = parse('max_threshold')
        if( type(min_threshold) != type(None) ):
            for i in range(len(misfits)):
                misfits[i] = min(min_threshold, misfits[i])

        if( type(max_threshold) != type(None) ):
            for i in range(len(misfits)):
                misfits[i] = min(max_threshold, misfits[i])

        Z,X = np.meshgrid(z,x)
        misfits = misfits.reshape(Z.shape)

        print('Z:\n%s'%Z)
        print('X:\n%s'%X)
        print('misfits:\n%s'%misfits)

        plt_title = mode_to_str(mode)
        if( l2_mode ):
            plt_title = "L2_surf_Fomel"
        else:
            plt_title = "Hminus1"
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

will_run=False

if(will_run):
    go()
