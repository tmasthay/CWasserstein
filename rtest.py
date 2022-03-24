#from rsf.proj import *
from seqflow import *
from forward import *
from os import system
from subprocess import check_output as co
import copy
import numpy as np
import signal
from base import create_base

SeqFlowLcl = SeqFlow
SeqPlotLcl = SeqPlot

global d_forward
global mode
global do_radon

do_radon = False

d_forward = create_base()

def set_mode(the_mode):
     global mode
     mode = the_mode

def get_val(s):
    try:
        return float(s.split(':')[-1].replace('\n','').replace(' ',''))
    except:
        return 1e10

def process(x, d):
    global do_radon

    t = time()
    print('Entering process')
    SeqFlowLcl(x + '_top', x,
            '''
            window n1=1 f1=0 | 
            put d2=%.4f unit2=s label2=Time
            '''%(d['dt']))

    SeqPlotLcl(x + '_top', x + '_top',
        '''
        grey color=seismic.csv scalebar=y title="%s_top"
        '''%x)

    if( do_radon ):
        SeqFlow(x + '_top_lrt', x + '_top',
            '''
            sftransp plane=12 | radon adj=y dp=4 p0=-800 np=400
            ''')
    
        SeqPlotLcl(x + '_top', x + '_top', 
                '''
                grey color=seismic.csv scalebar=y title="%s_top"
                '''%x)
        SeqPlotLcl(x + '_top_lrt', x + '_top_lrt', 
                '''
                grey color=seismic.csv scalebar=y title="%s_top_lrt"
                '''%x)
    t4 = time()
    print('Exiting process...total time: %.4f seconds'%(t4 - t))

def combine(x,y,d):
    #retrieve global variables
    global mode
    global do_radon

    #get output file name
    comparison = d['output']


    #get surface data input file name
    t0 = x + '_top'
    t1 = y + '_top'

    if(do_radon):
        b0 = x + '_top_lrt'
        b1 = y + '_top_lrt'
        inputs = '.vpl '.join([t0,t1,b0,b1, ''])[:-1]
        system('vppen gridnum=2,2 %s > %s.vpl'%(inputs, comparison))

    #define command for Wasserstein distance
    t_mode = 't.rsf'
    cmd='./ot.exe g=${SOURCES[1]} t=%s mode=%d'%(t_mode, mode)


    #take Wasserstein distance between radon transforms 
    if( do_radon ):
        b0 = x + '_top_lrt'
        b1 = y + '_top_lrt'
        c0 = t0 + '_trans'
        c1 = t1 + '_trans'
        plane=23
        SeqFlowLcl(c0, b0, 'transp plane=%d'%plane)
        SeqFlowLcl(c1, b1, 'transp plane=%d'%plane)
        SeqFlowLcl('dist', '%s.rsf %s.rsf'%(c0, c1), cmd)
    #else take Wasserstein distance between original signals
    else:
        c0 = t0 + '_trans'
        c1 = t1 + '_trans'
        plane=23
        SeqFlowLcl(c0, t0, 'transp plane=%d'%plane)
        SeqFlowLcl(c1, t1, 'transp plane=%d'%plane)
        SeqFlowLcl('dist', '%s.rsf %s.rsf'%(c0, c1), cmd)

    return get_val(co('sfdisfil < dist.rsf', shell=True).decode('utf-8'))

def run_case(d, reference_name):
    global d_veltran

    forward(d)
    namez = 'wavz_' + d['case']
    namex = 'wavx_' + d['case']

#    process(namez, d_veltran)
#    process(namex, d_veltran)
    process(namez, d)
    process(namex, d)

    d_output = {'output': 'compare_' + d['case']}
    zc = combine(namez, 'wavz_' + reference_name, {'output': 'cmp_z_' + d['case']})
    xc = combine(namex, 'wavx_' + reference_name, {'output': 'cmp_x_' + d['case']})
    if( type(zc) == type(None) or type(xc) == type(None) ):
        distance = 1e10
    else:
        distance = zc + xc

    os.system('rm /var/tmp/%s*.rsf@'%namez)
    os.system('rm /var/tmp/%s*.rsf@'%namex)

    return distance

def run_case_short(src_z, src_x):
    global d_forward
    global mode

    e = copy.copy(d_forward)
    e['case'] = str(src_z) + '_' + str(src_x) + '_' + str(mode)
    e['sszf'] = src_z
    e['ssxf'] = src_x

    #input('run_case_short_args = \n%s\n%s'%(e,d_forward['case']))
    return run_case(e, d_forward['case'])

def run_test(my_pts):
    misfits = []
    N = len(my_pts)
    for (i,pt) in enumerate(my_pts):
        print(i)
        t = time()
        misfits.append(run_case_short(pt[0], pt[1]))
        print('Iteration %d of %d took %.4f seconds'%(i,N,time() - t))
    return np.array(misfits)

if( not os.path.exists('x.rsf') ):
    SeqFlow('x',None,
        '''math output="x1" n1=%d d1=%f'''%(d_forward['nx'],d_forward['dx']))

    SeqFlow('t',None,
       '''math output="x1" n1=%d d1=%f'''%(d_forward['nt'], d_forward['dt']))

if( not os.path.exists('wavz_' + d_forward['case'] + '.rsf') ):
    forward(d_forward)
    process('wavz_' + d_forward['case'], d_forward)
    process('wavx_' + d_forward['case'], d_forward)
    print('wavz_' + d_forward['case'] + '_top.rsf should exist')
