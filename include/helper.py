import numpy as np
from itertools import product
import copy
from forward import forward
from subprocess import check_output as co
from rsf.proj import *
import re
import base

def grey(title, 
    color='seismic', 
    scalebar='y', 
    label1='Time', 
    unit1='s', 
    label2='RecLoc', 
    unit2='km',
    minval='',
    maxval=''):
    s = '''grey color=%s scalebar=%s title="%s" 
        label1=%s unit1=%s label2=%s unit2=%s
        wheretitle=t wherexlabel=b
        '''%(
        color, scalebar,title, label1, unit1, label2, unit2)

    if( len(minval) > 0 ):
        s += ' minval=%s'%(minval)
    if( len(maxval) > 0 ):
        s += ' maxval=%s'%(maxval)
    return s

def snapshot_movie(base_name, ref_name, times,
    axis=3,
    title_pre='',
    title_suf='', 
    color='seismic', 
    scalebar='y', 
    label1='Depth', 
    unit1='km', 
    label2='RecLoc', 
    unit2='km',
    minval='',
    maxval=''):
    for t in times:
        Flow('%s%d'%(base_name, t), ref_name,
            '''
            window n%d=1 f%d=%d
            '''%(axis, axis, t))
        Result('%s%d'%(base_name,t), '%s%d'%(base_name,t), 
            grey('%s%d%s'%(title_pre, t, title_suf),
                color,
                scalebar,
                label1,
                unit1,
                label2,
                unit2,
                minval,
                maxval))

def surface_obs(name,
    color='seismic', 
    scalebar='y', 
    label1='Time', 
    unit1='s', 
    label2='X', 
    unit2='km',
    minval='',
    maxval=''):
    without_numbers = re.sub('[0-9][0-9]*', '', name)
    number = name[len(without_numbers):]
    final_name = without_numbers + 'top' + number
    Flow(final_name, name,
        '''
        window n1=1 f1=0 |
        put label1=%s label2=%s unit1=%s unit2=%s
        ''')
    Result(final_name, final_name,
        grey('Surface obs for case %s'%number,
            color,
            scalebar,
            label1,
            unit1,
            label2,
            unit2,
            minval,
            maxval))
    return final_name
   
def landscape(d_ref, d_perturb):
    s = ' '.join(list(d_perturb[0].keys()))
    for (i,d) in enumerate(d_perturb):
        #later on check for d_ref =''...meaning we don't want to compare to synthetic data
        d_curr = copy.copy(d_ref)
        d_curr.update(d)
        forward(d_curr)
    return s

def proto_landscape(d):
    d_ref = d['ref']
    z = np.linspace(d['zs'], d['ze'], d['nz'])
    x = np.linspace(d['xs'], d['xe'], d['nx'])
    total = [[(zz,xx) for xx in x] for zz in z]
    d_perturb = []
    i = 0
    for zz in z:
        for xx in x:
            i += 1
            d_perturb.append({'sszf': zz, 
                'eszf': zz,
                'ssxf': xx,
                'esxf': xx,
                'name': '%d'%i})

    return landscape(d_ref, d_perturb)

def get_val(file_name):
   val = co('sfdisfil < %s'%file_name, shell=True) \
       .decode('utf-8') \
       .split(':')[-1] \
       .replace(' ','')
   return float(val)

def get_norm(file_name):
    return float(co('sfattr < %s.rsf'%file_name.replace('.rsf',''), 
        shell=True) \
            .decode('utf-8') \
            .split('\n')[3] \
            .split('=')[-1] \
            .replace(' ',''))

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

def modify_src(zz, xx, cmd):
    cmd1 = re.sub('sszf=.*? ', 'sszf=%f '%zz, cmd)
    cmd2 = re.sub('eszf=.*? ', 'eszf=%f '%zz, cmd1)
    cmd3 = re.sub('ssxf=.*? ', 'ssxf=%f '%xx, cmd2)
    cmd4 = re.sub('esxf=.*? ', 'esxf=%f '%xx, cmd3)
    return cmd4

def modify_of(i, of):
   of1 = re.sub('wavz', 'wavz%d'%i, of)
   return re.sub('wavx', 'wavx%d'%i, of1)

def create_reference_data():
    d = base.create_base()
    output_files, input_files, fc = forward(d)
    Flow(output_files, input_files, fc)

    wavz_syn = output_files.split(' ')[0]
    wavx_syn = output_files.split(' ')[1]

    wavz_syn_fourier = attach(wavz_syn, 'fft')
    wavx_syn_fourier = attach(wavx_syn, 'fft')

    fourier_cmd = 'fft1 | fft3 sign=1 | real'

    Flow(wavz_syn_fourier, wavz_syn, fourier_cmd)
    Flow(wavx_syn_fourier, wavx_syn, fourier_cmd)

    output_files = re.sub('_synthetic', '', output_files)

    return output_files, input_files, fc, d

def create_synthetic_data(z, x, modes, np_factor=1.0):
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

    for mode in modes:
        i = 0
        for zz in z:
            for xx in x:
                w2_redun = ('w2' in mode and w2_created)
                l2_redun = ('sobolev' == mode and s == 0.0 and w2_created)
                sob_redun = ('sobolev' == mode and sobolev_created)
 
                if(w2_redun or l2_redun or sob_redun):
                    continue

                if(w2_redun or l2_redun or sob_redun):
                    print('''
                        SHOULD NEVER SEE THIS TEXT: USE OF "continue" 
                        KEYWORD IS FLAWED
                        ''')
                #increment case number
                i += 1
                 
                #do string manipulation to create new file
                new_fc = modify_src(zz, xx, fc) 
                new_output_files = modify_of(i, output_files)
        
                #create new synthetic data
                Flow(new_output_files, input_files, new_fc)
        
                #get surface observations
                wavz_curr = new_output_files.split(' ')[0]
                wavx_curr = new_output_files.split(' ')[1]
        
                if( not sobolev_created \
                    and mode == 'sobolev' \
                    and sobolev_norm > 0.0 ):
                    wavz_orig = wavz_curr
                    wavx_orig = wavx_curr
                    wavz_curr = attach(wavz_curr, 'fft')
                    wavx_curr = attach(wavx_curr, 'fft')
         
                    Flow(wavz_curr, wavz_orig, fourier_command)
                    Flow(wavx_curr, wavx_orig, fourier_command)

                    allz_sobolev += wavz_curr + '.rsf'
                    allx_sobolev += wavx_curr + '.rsf'
 
                    sobolev_created = True 
        
                if( not w2_created \
                    and ('w2' in mode or \
                        ('sobolev' == mode and sobolev_norm == 0.0))): 
                    allz_w2 += wavz_curr + '.rsf '
                    allx_w2 += wavx_curr + '.rsf '
    if( len(allz_w2) > 0 ):
        Flow('allz_w2', None, 'sfcat %s'%allz_w2)
        Flow('allx_w2', None, 'sfcat %s'%allx_w2)
    if( len(allx_w2) > 0 ):
        Flow('allz_sobolev', None, 'sfcat %s'%allz_sobolev)
        Flow('allx_sobolev', None, 'sfcat %s'%allx_sobolev)

