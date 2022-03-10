#from rsf.proj import *
from seqflow import *
from forward import *
from os import system
from subprocess import check_output as co
import copy
import numpy as np
import signal

SeqFlowLcl = SeqFlow
SeqPlotLcl = SeqPlot

global d_forward
global d_veltran
global mode

N = 100

nz = N
dz = 1.0 / nz

nx = N
dx = 1.0 / nx

nt = 400
dt = 5e-03

d_forward = {
        'case' : 'synthetic',
        'nz' : nz,
        'dz' : dz,
        'nx' : nx,
        'dx' : dx,
        'dt' :  dt,
        'nt' :  nt,
        'sszf' :  0.5,
        'ssxf' :  0.5,
        'eszf' :  0.5,
        'esxf' :  0.5,
        'nsz' :  1,
        'nsx' :  1,
        'nb' : 35,
        'fm' : 5.0,
        'vp': '0.5',
        'vs': '0.707 * input',
        'rho': '1.0'
}

d_veltran = {
        'v0': 0.5,
        'dv': 0.01,
        'adj': 'n',
        'nx': 100,
        'dx': 0.01,
        'x0': 0
}

def set_mode(the_mode):
     global mode
     mode = the_mode

def get_val(s):
    try:
        return float(s.split(':')[-1].replace('\n','').replace(' ',''))
    except:
        print('Sending alarm after failed to convert float')
        signal.alarm(5)

def process(x, d):
    t = time()
    print('Entering process')
    SeqFlowLcl(x + '_top', x,
            '''
            window n1=1 f1=0 | 
            put v0=%f dv=%f
            '''%(d['v0'], d['dv']))
    t1 = time()
    print('top time: %.4f'%(t1 - t))
    #input(x)

#    SeqFlowLcl(x + '_top_hrt', x + '_top', 
#            '''
#            veltran adj=%s nx=%d dx=%f x0=%f
#            '''%(d['adj'], d['nx'], d['dx'], d['x0']))

    SeqFlow(x + '_top_lrt', x + '_top',
        '''
        radon adj=y dp=0.1 p0=-5.0 np=100
        ''')
    t2 = time()
    print('radon time: %.4f'%(t2 - t1))

    SeqPlotLcl(x + '_top', x + '_top', 
            '''
            grey color=seismic.csv scalebar=y title="%s_top"
            '''%x)
    t3 = time()
    print('top plot time: %.4f'%(t3 - t2))
#    SeqPlotLcl(x + '_top_hrt', x + '_top_hrt', 
#            '''
#            grey color=seismic.csv scalebar=y title="%s_top_hrt"
#            '''%x)
    SeqPlotLcl(x + '_top_lrt', x + '_top_lrt', 
            '''
            grey color=seismic.csv scalebar=y title="%s_top_lrt"
            '''%x)
    t4 = time()
    print('lrt plot time: %.4f'%(t4-t3))
    ##input('yo4')
    print('Exiting process...total time: %.4f seconds'%(t4 - t))

def combine(x,y,d):
    global mode

    comparison = d['output']
    t0 = x + '_top'
    t1 = y + '_top'
    b0 = x + '_top_lrt'
    b1 = y + '_top_lrt'
    inputs = '.vpl '.join([t0,t1,b0,b1, ''])[:-1]
    system('vppen gridnum=2,2 %s > %s.vpl'%(inputs, comparison))

    t_mode='x.rsf'
    #t_mode = 't.rsf'
    cmd='./ot.exe g=${SOURCES[1]} t=%s mode=%d'%(t_mode, mode)

    c0 = t0 + '_trans'
    c1 = t1 + '_trans'
    plane=13 if t_mode[0] == 'x' else 23
#    SeqFlowLcl(c0, t0, 'transp plane=23 | transp plane=12')
#    SeqFlowLcl(c1, t1, 'transp plane=23 | transp plane=12')
    SeqFlowLcl(c0, t0, 'transp plane=%d'%plane)
    SeqFlowLcl(c1, t1, 'transp plane=%d'%plane)
    SeqFlowLcl('dist', '%s.rsf %s.rsf'%(c0, c1), cmd)
    #SeqFlowLcl('dist', '%s.rsf %s.rsf'%(c0, c1), cmd)
    #SeqFlowLcl('dist', '%s.rsf %s.rsf'%(t0, t1), cmd.replace('x.rsf','t.rsf'))

    return get_val(co('sfdisfil < dist.rsf', shell=True).decode('utf-8'))

def run_case(d, reference_name):
    global d_veltran

    forward(d)
    namez = 'wavz_' + d['case']
    namex = 'wavx_' + d['case']

    process(namez, d_veltran)
    process(namex, d_veltran)

    d_output = {'output': 'compare_' + d['case']}
    zc = combine(namez, 'wavz_' + reference_name, {'output': 'cmp_z_' + d['case']})
    xc = combine(namex, 'wavx_' + reference_name, {'output': 'cmp_x_' + d['case']})
    distance = zc + xc

    #os.system('rm /var/tmp/%s*.rsf@'%namez)
    #os.system('rm /var/tmp/%s*.rsf@'%namex)

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
    process('wavz_' + d_forward['case'], d_veltran)
    process('wavx_' + d_forward['case'], d_veltran)
