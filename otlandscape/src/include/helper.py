import numpy as np
from itertools import product
import copy
from forward import forward
from subprocess import check_output as co
from rsf.proj import *
import re
import base
     
def attach(s, suf):
    t = re.sub('[0-9]*', '', s)
    return s[:len(t)] + suf + s[len(t):]
    
def dict_to_str(d):
    ignore_keys=[]
    s = ''
    for key in d.keys():
        if( key not in ignore_keys ):
            s += '%s=%s '%(key,str(d[key]))
    return s
    
def get_mode(branch):
    branch = branch.lower()
    if( branch == 'sobolev' ):
        return 0
    elif( branch == 'w2split' ):
        return 1
    elif( branch == 'w2square' ):
        return 2
    elif( branch == 'w2linexp' ):
        return 3
    elif( branch == 'w2linear' ):
        return 4
    elif( branch == 'w2exp' ):
        return 5
    return -1
    
def get_scons_field(the_name):
    return co('''grep "^%s[ ]*=" SConstruct'''%the_name, shell=True) \
    .decode('utf-8') \
    .replace('\n', '') \
    .split('=')[-1]
     
def create_reference_data():
    d = base.create_base()
    output_files, input_files, fc = forward(d)
    Flow(output_files, input_files, fc)
    
    wavz_syn = output_files.split(' ')[0].replace('stacked', '')
    wavx_syn = output_files.split(' ')[1].replace('stacked', '')
    
    Flow(wavz_syn, wavz_syn + 'stacked', 'window f3=0 n3=1')
    Flow(wavx_syn, wavx_syn + 'stacked', 'window f3=0 n3=1')
    
    wavz_syn_fourier = attach(wavz_syn, 'fft')
    wavx_syn_fourier = attach(wavx_syn, 'fft')
    
    fourier_cmd = 'fft1 | fft3 sign=1 | real'
    
    Flow(wavz_syn_fourier, wavz_syn, fourier_cmd)
    Flow(wavx_syn_fourier, wavx_syn, fourier_cmd)
    
    output_files = re.sub('_synthetic', '', output_files)
    
    return output_files, input_files, fc, d
 
def create_synthetic_data(z, x, modes, ps, np_factor=1.0):
    w2_created = False
    sobolev_created = False
    
    output_files, input_files, fc, d = create_reference_data()
    
    #create p
    n_probs = int(np_factor * d['nt'])
    Flow('p', None, 
        'math output="x1" n1=%d d1=%.8f o1=0.0'%(n_probs, 1.0 / n_probs))
    
    #create t
    Flow('t', None, 
        'math output="x1" n1=%d d1=%.8f o1=0'%(d['nt'], d['dt']))
    
    #define associated lists for w2 and sobolev cases
    allz_w2 = ''
    allx_w2 = ''
    
    allz_sobolev = ''
    allx_sobolev = ''
  
    z = np.linspace(z[0], z[1], z[2])
    x = np.linspace(x[0], x[1], x[2]) 
    for (i,zz) in enumerate(z):
        for (j,xx) in enumerate(x):
            d.update({'sszf': 'sszf_%d_%d.rsf'%(i,j), \
                'ssxf': 'ssxf_%d_%d.rsf'%(i,j), \
                'eszf': 'eszf_%d_%d.rsf'%(i,j), \
                'esxf': 'esxf_%d_%d.rsf'%(i,j), \
                'sszf_val': zz, 'eszf_val': zz, 'numz_comp': 1, \
                'ssxf_val': xx, 'esxf_val': xx, 'numx_comp': 1, \
                'case': 'shifted_%d_%d'%(i,j)})
            
            output_files, input_files, fc = forward(d) 
            Flow(output_files, input_files, fc)
            for l in output_files.split():
                Result(l, l,
                    ''' 
                    transp plane=12 | 
                    grey scalebar=y color=i
                    ''')

def ricker_shift(wass_exe):
    def get_ricker(f0, t0, sigma, amp, domain='input'):
        if( len(f0) == 1 ):
            shift_term = '(%s - %.4f)'%(domain, t0[0])
            shift_term = '%s * %s'%(shift_term, shift_term)
            sig_inv_squared = 1.0 / sigma[0]**2
            s = '%.4f * (1 - (%.4f * %s)) * exp(-%.4f * %s)'%(
                amp[0], sig_inv_squared, shift_term,
                0.5 * sig_inv_squared, shift_term)
            return s
        else:
            return '(%s) + (%s)'%(
                 get_ricker(f0[:1], t0[:1], sigma[:1], amp[:1], domain),
                 get_ricker(f0[1:], t0[1:], sigma[1:], amp[1:], domain))

    f0 = [20.0, 20.0]
    t0 = np.array([0.0, 2.0])
    sigma = [1.0, 0.5]
    amp = [1.0, 1.0]

    t_min = -5.0
    t_max = 5.0
    nt = 1000
    dt = (t_max - t_min) / nt
   
    Flow('tricker', None, 'math output="x1" n1=%d d1=%.4f'%(nt, dt))
    Flow('ricker_ref', 'tricker', 
        '''
        math output="%s"
            | put n2=1
            | transp plane=12'''%(get_ricker(f0,t0,sigma,amp)))

    shift_min = -3
    shift_max = 3
    num_shifts = 25
    shifts = np.linspace(shift_min, shift_max, num_shifts)

    k = 0
    for i in range(len(shifts)):
        for j in range(len(shifts)):
            k += 1
            t_curr = np.zeros(len(t0))
            t_curr[0] = t0[0] + shifts[i] 
            t_curr[1] = t0[1] + shifts[j]

            Flow('ricker%d'%k, 'tricker',
            '''
            math output="%s"
                | put n2=1
                | transp plane=12'''%(get_ricker(f0,t_curr,sigma,amp)))

    num_cases = len(shifts)**2
    #all_rickers = ' '.join(['ricker%d'%i for i in range(1,num_cases+1)])
    all_rickers = ' '.join(['ricker%d.rsf'%i for i in range(1,num_cases+1)])
    srcs = ' '.join(['${SOURCES[%d]}'%i for i in range(1,num_cases+1)])

    

#    Flow('ricker_shift', all_rickers, 'cat %s'%srcs)
    Flow('ricker_shift', None, 'cat %s'%all_rickers)

    np_factor = 10.0
    prob_samples = np_factor * nt
    Flow('pricker', None, 
        'math output="x1" d1=%.4f n1=%d'%(1.0 / prob_samples,
            prob_samples))

    wc = lambda m : ''' 
    ./${SOURCES[4]} 
    data=${SOURCES[1]} 
    t=${SOURCES[2]} 
    p=${SOURCES[3]}
    bz=%.4f
    ez=%.4f
    bx=%.4f
    ex=%.4f
    zs=%d
    xs=%d
    mode=%s
    '''%(shifts[0], shifts[-1], shifts[0], shifts[-1], len(shifts), len(shifts), m)

    modes = ['w2split', 'sobolev']
    mode_nos = [get_mode(m) for m in modes]
    ps = [{}, {'s': 0.0}]
    
    for (i, mode) in enumerate(modes):
        cmd = wc(mode_nos[i]) + ' ' + dict_to_str(ps[i])
        cmd = cmd[:-1]
        Flow('districker%d'%i, 
            'ricker_shift ricker_ref tricker pricker %s'%wass_exe,
            cmd)
        Result('districker%d'%i, 'districker%d'%i, 
            '''
            grey scalebar=y label1=Pwaveshift label2=Swaveshift title="%s"
            '''%(mode + ' ' + dict_to_str(ps[i])))
            
        
   
