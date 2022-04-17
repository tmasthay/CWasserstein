from rsf.proj import *
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
mode='w2split'
mode_no = get_mode(mode)
sobolev_norm = 0.0

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

"""
z_syn_top = surface_obs(wavz_syn)
x_syn_top = surface_obs(wavx_syn)
"""


#pull out "synthetic" for future computations
output_files = re.sub('_synthetic', '', output_files)

#create z for convexity
zs = 0.4
ze = 0.6
nz = 25
z = np.linspace(zs, ze, nz)

#create x for convexity
xs = 0.2
xe = 0.8
nx = 25
x = np.linspace(xs, xe, nx)

#creat convexity cases by string subs on command
i = 0

for zz in z:
    for xx in x:
        #increment case number
        i += 1
  
        #do string manipulation to create new file
        fc1 = re.sub('sszf=.*? ', 'sszf=%f '%zz, fc)
        fc2 = re.sub('eszf=.*? ', 'eszf=%f '%zz, fc1)
        fc3 = re.sub('ssxf=.*? ', 'ssxf=%f '%xx, fc2)
        new_fc = re.sub('esxf=.*? ', 'esxf=%f '%xx, fc3)
        of1 = re.sub('wavz', 'wavz%d'%i, output_files)
        new_output_files = re.sub('wavx', 'wavx%d'%i, of1)

        #create new synthetic data
        Flow(new_output_files, input_files, new_fc)

        #get surface observations
        wavz_curr = new_output_files.split(' ')[0]
        wavx_curr = new_output_files.split(' ')[1]

        if( mode == 'sobolev' ):
            wavz_curr_fourier = attach(wavz_curr, 'fft')
            wavx_curr_fourier = attach(wavx_curr, 'fft')
 
            Flow(wavz_curr_fourier, wavz_curr, fourier_command)
            Flow(wavx_curr_fourier, wavx_curr, fourier_command)

        #plot_str = 'put d1=%.4f d2=%.4f'%(d['dx'], d['dt'])
        plot_str = 'transp plane=12 | '
        plot_str += grey('%.2f,%.2f'%(zz,xx), color='i')
        Result(wavz_curr, wavz_curr, plot_str)
        """
        z_curr_top = surface_obs(wavz_curr)
        x_curr_top = surface_obs(wavx_curr)
        """

        wass_cmd = '''
            ./${SOURCES[4]} data=${SOURCES[1]} t=${SOURCES[2]} p=${SOURCES[3]} 
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
                 w2_cmd)

proj.End()
