import numpy as np
import copy
import sys
sys.path.append('./include')

from forward import forward
import base
from helper import *
proj = Project()

#define libs and cflags
libs=['rsf++']
cflags='-I. -I./include -w'
#cflags=''

#add to compilation command
proj.Prepend(LIBS=libs)
proj.Replace(CCFLAGS=cflags)

#compile programs
#elas = proj.Program('./include/myelastic.cc')
elas = proj.Program('./include/elastic_top.c')
driv = proj.Program('driver.cc')
wass_exe = str(driv[0])

"""
#create convexity cases by string subs on command
i = 0

#all cases 
allz = ''
allx = ''

for zz in z:
    for xx in x:
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

        if( mode == 'sobolev' and sobolev_norm > 0.0 ):
            wavz_orig = wavz_curr
            wavx_orig = wavx_curr
            wavz_curr = attach(wavz_curr, 'fft')
            wavx_curr = attach(wavx_curr, 'fft')
 
            Flow(wavz_curr, wavz_orig, fourier_command)
            Flow(wavx_curr, wavx_orig, fourier_command)            

        #plot_str = 'put d1=%.4f d2=%.4f'%(d['dx'], d['dt'])
        plot_str = 'transp plane=12 | '
        plot_str += grey('%.2f,%.2f'%(zz,xx), color='i')
        #Result(wavz_curr, wavz_curr, plot_str)

        allz += wavz_curr + '.rsf '
        allx += wavx_curr + '.rsf '

allz = allz[:-1]
allx = allx[:-1]

Flow('allz', None, 'sfcat %s'%allz)
Flow('allx', None, 'sfcat %s'%allx)
"""

modes = ['w2square']
mode_nos = [get_mode(m) for m in modes]
ps = [dict()]

#add other param for setting nz,nx simultaneously
n_all = 25
#create z for convexity
zs = 0.4
ze = 0.6
nz = n_all
z = np.linspace(zs, ze, nz)

#create x for convexity
xs = 0.2
xe = 0.8
nx = n_all
x = np.linspace(xs, xe, nx)

#define probability space granularity
np_factor = 1.0
create_synthetic_data(z, x, modes, np_factor)

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
    '''%(z[0], z[-1], x[0], x[-1], nz, nx, m)

all_colors = '''cubeyf1 gist_earth izoaz linearlfb lb montag owb rwb seismic spectral
 viridis acton bamako batlow berlin bilbao broc brocO buda cork corkO 
 davos devon grayC hawaii imola inferno lajolla lapaz lisbon magma nuuk 
 oleron oslo plasma roma romaO tofino tokyo turku vik vikO'''.split(' ')

wavz_syn = 'wavz_synthetic'
wavx_syn = 'wavx_synthetic'
wavz_syn_fourier = wavz_syn + 'fft'
wavx_syn_fourier = wavx_syn + 'fft'

for (i,mode) in enumerate(modes):
    cmd = wc(mode_nos[i]) + ' ' + dict_to_str(ps[i])
    cmd = cmd[:-1]
    if( mode == 'sobolev' ):
        Flow('distz_%d'%i, 
            'allz_sobolev %s t p %s'%(wavz_syn_fourier, wass_exe),
            cmd)
        Flow('distx_%d'%i, 
            'allx_sobolev %s t p %s'%(wavx_syn_fourier, wass_exe),
            cmd)
    elif( 'w2' in mode ):
        Flow('distz_%d'%i, 'allz_w2 %s t p %s'%(wavz_syn, wass_exe), cmd)
        Flow('distx_%d'%i, 'allx_w2 %s t p %s'%(wavx_syn, wass_exe), cmd)
    Flow('fulldist_%d'%i, 'distz_%d distx_%d'%(i,i), 'add ${SOURCES[1]}')
    clr='viridis'
    Result('fulldist%d'%i, 'fulldist%d'%i, 
        '''
        grey color=%s scalebar=y title="%s" 
            label1=Zshift label2=Xshift
        '''%(clr, mode))




proj.End()
