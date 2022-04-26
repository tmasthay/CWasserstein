import numpy as np
from itertools import product
import copy
from forward import forward
from subprocess import check_output as co
from rsf.proj import *
import re

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
            s += '%s=%s\n'%(key,str(d[key]))
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
    
