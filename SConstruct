from rsf.proj import *
from subprocess import check_output as co
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
elas = proj.Program('./include/myelastic.cc')
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
z_syn_top = surface_obs(wavz_syn)
x_syn_top = surface_obs(wavx_syn)


#pull out "synthetic" for future computations
output_files = re.sub('_synthetic', '', output_files)

#create z for convexity
zs = 0.2
ze = 0.8
nz = 5
z = np.linspace(zs, ze, nz)

#create x for convexity
xs = 0.2
xe = 0.8
nx = 5
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
        z_curr_top = surface_obs(wavz_curr)
        x_curr_top = surface_obs(wavx_curr)

        wass_cmd = '''
            ./${SOURCES[4]} data=${SOURCES[1]} t=${SOURCES[2]} p=${SOURCES[3]} 
            '''
        #get distance
	Flow('zdist%d'%i, '%s %s t p %s'%(z_curr_top, z_syn_top, wass_exe),
            wass_cmd)

        Flow('xdist%d'%i, '%s %s t p %s'%(x_curr_top, x_syn_top, wass_exe),
            wass_cmd)

proj.End()
