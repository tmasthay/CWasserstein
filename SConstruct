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

#define mode of execution
#mode='sobolev' if len(sys.argv) == 1 else sys.argv[1].split('mode=')[-1]
mode='w2square'
mode_no = get_mode(mode)
sobolev_norm = -1.0

#define parameters
cvals = 10.0**np.array(range(-5,5))
params = dict({'c1': cvals[9]})
if( 'exp' not in mode ): params = dict()
param_string = dict_to_str(params)

#add to compilation command
proj.Prepend(LIBS=libs)
proj.Replace(CCFLAGS=cflags)


#compile programs
#elas = proj.Program('./include/myelastic.cc')
elas = proj.Program('./include/elastic_top.c')
driv = proj.Program('driver.cc')
wass_exe = str(driv[0])

#create base case
d = base.create_base()

#get forward command and synthetic output files
output_files, input_files, fc = forward(d)

#create p
n_probs = 100
Flow('p', None, 'math output="x1" n1=%d d1=%.8f o1=0.0'%(n_probs, 1.0 / n_probs))

#create t
Flow('t', None, 'math output="x1" n1=%d d1=%.8f o1=0'%(d['nt'], d['dt']))

#create reference data then format command for synthetic
Flow(output_files, input_files, fc)

#get surface observations
wavz_syn = output_files.split(' ')[0]
wavx_syn = output_files.split(' ')[1]

wavz_syn_fourier = attach(wavz_syn, 'fft')
wavx_syn_fourier = attach(wavx_syn, 'fft')

fourier_command = 'fft1 | fft3 sign=1 | real'

Flow(wavz_syn_fourier, wavz_syn, fourier_command)
Flow(wavx_syn_fourier, wavx_syn, fourier_command)

#pull out "synthetic" for future computations
output_files = re.sub('_synthetic', '', output_files)

#create z for convexity
zs = 0.4
ze = 0.6
nz = 4
z = np.linspace(zs, ze, nz)

#create x for convexity
xs = 0.2
xe = 0.8
nx = 4
x = np.linspace(xs, xe, nx)

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

        """
        wass_cmd = '''
            ./${SOURCES[4]} data=${SOURCES[1]} t=${SOURCES[2]}
            p= 
            '''
        if( mode == 'sobolev' ):
            sobolev_cmd = wass_cmd + ' mode=%d s=%.4f'%(mode_no,
                sobolev_norm)
            #get distance
	    Flow('zdist%d'%i, 
                '%s %s t p %s'%(
                    wavz_curr if sobolev_norm == 0.0 \
                        else wavz_curr_fourier, 
                    wavz_syn if sobolev_norm == 0.0 \
                        else wavz_syn_fourier, 
                    wass_exe),
                sobolev_cmd)

            Flow('xdist%d'%i, 
                '%s %s t p %s'%(
                    wavx_curr if sobolev_norm == 0.0 \
                        else wavz_curr_fourier, 
                    wavx_syn if sobolev_norm == 0.0 \
                        else wavx_syn_fourier, 
                    wass_exe),
                 sobolev_cmd)
        elif(mode[0:2] == 'w2'):
             w2_cmd = wass_cmd + ' mode=%d'%mode_no
             w2_cmd = w2_cmd + ' %s'%param_string
             Flow('zdist%d'%i,
                 '%s %s t p %s'%(
                     wavz_curr,
                     wavz_syn,
                     wass_exe),
                 w2_cmd)
             Flow('xdist%d'%i,
                 '%s %s t p %s'%(
                     wavx_curr,
                     wavx_syn,
                     wass_exe),
                 w2_cmd)"""

allz = allz[:-1]
allx = allx[:-1]

Flow('allz', None, 'sfcat %s'%allz)
Flow('allx', None, 'sfcat %s'%allx)

wass_cmd = '''
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
    '''%(z[0], z[-1], x[0], x[-1], nz, nx, mode_no)

if( mode == 'sobolev' ):
    Flow('distz', 'allz %s t p %s'%(wavz_syn_fourier, wass_exe),
        wass_cmd + ' s=%f'%sobolev_norm)
    Flow('distx', 'allx %s t p %s'%(wavx_syn_fourier, wass_exe),
        wass_cmd + ' s=%f'%sobolev_norm)
elif( 'w2' in mode ):
    wass_cmd = wass_cmd + ' ' + param_string
    Flow('distz', 'allz %s t p %s'%(wavz_syn, wass_exe), wass_cmd)
    Flow('distx', 'allx %s t p %s'%(wavx_syn, wass_exe), wass_cmd)

Flow('fulldist', 'distz distx', 'add ${SOURCES[1]}')


all_colors = '''cubeyf1 gist_earth izoaz linearlfb lb montag owb rwb seismic spectral
 viridis acton bamako batlow berlin bilbao broc brocO buda cork corkO 
 davos devon grayC hawaii imola inferno lajolla lapaz lisbon magma nuuk 
 oleron oslo plasma roma romaO tofino tokyo turku vik vikO'''.split(' ')

clr='viridis'
Result('fulldist', 'fulldist', 'grey color=%s scalebar=y title="%s" label1=Zshift label2=Xshift'%(clr, mode))

proj.End()
