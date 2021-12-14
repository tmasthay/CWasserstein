from seqflow import *
from forward import *
from os import system
from input_dict import get_dict
import numpy as np
from subprocess import check_output as co
import matplotlib.pyplot as plt
def create_top(the_name):
     SeqFlowV('%s_top.rsf'%the_name, '%s.rsf'%the_name, 'window f1=0 n1=1')

def create_idv_plots(the_name):
     SeqFlowV('%s.vpl'%the_name, '%s.rsf'%the_name, 'grey')

def create_idv(name, sz, sx):
    tru_name = 'wavz_' + name
    d = get_dict(name, sz, sx)
    forward(d)
    create_top(tru_name)
    create_idv_plots(tru_name)
    create_idv_plots(tru_name + '_top')
    SeqFlowV(tru_name + '_top_trans.rsf', tru_name + '_top.rsf', 'transp plane=23')

def create_comp_plots(name1, name2):
    name1 = 'wavz_' + name1
    name2 = 'wavz_' + name2
    names_vpl = '%s.vpl %s.vpl'%(name1, name2)
    system('vppen gridnum=2,1 %s > compare_%s.vpl'%(names_vpl, name2))

def w2(name1, name2):
    name1 = 'wavz_' + name1 + '_trans'
    name2 = 'wavz_' + name2 + '_trans'
    cmd = './ot.exe g=%s.rsf t=t.rsf mode=1 < %s.rsf > dist.rsf'%(name1,name2)
    print('cmd = %s'%cmd)
    system(cmd)
    s = co('sfdisfil < dist.rsf', shell=True) \
        .decode('utf-8') \
        .replace('\n', '') \
        .replace(' ','') \
        .split(':')[1]
    return float(s)

def run_script():
    create_idv('ref', 0.5, 0.5)

    w2_lcl = lambda x : w2('ref_top', x)

    N = 20
    t0 = 0.05
    t1 = 0.95 
    v_tests = np.linspace(t0, t1 ,N)
    h_tests = np.linspace(t0, t1,N)

    vw2 = np.zeros(len(v_tests))
    hw2 = np.zeros(len(h_tests))

    for (i,v) in enumerate(v_tests):
        name = 'v_%.2f'%v
        create_idv(name, v, 0.5)
        create_comp_plots('ref_top', name + '_top')
        vw2[i] = w2_lcl(name + '_top')
   
    for (i,h) in enumerate(h_tests):
        name = 'h_%.2f'%h
        create_idv(name, 0.5, h)
        create_comp_plots('ref_top', name + '_top') 
        hw2[i] = w2_lcl(name + '_top')

    plt.plot(v_tests, vw2)
    plt.title('Vertical Test')
    plt.xlabel('Shift')
    plt.ylabel('Misfit')
    plt.savefig('vertical.png')
    plt.clf()

    plt.plot(h_tests, hw2)
    plt.title('Horizontal Test')
    plt.xlabel('Shift')
    plt.ylabel('Misfit')
    plt.savefig('horizontal.png')
    
run_script()
