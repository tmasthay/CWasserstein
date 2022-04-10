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
cflags=''

#add to compilation command
proj.Prepend(LIBS=libs)
proj.Replace(CCFLAGS=cflags)


prog = proj.Program('driver.cc')

elas = proj.Program('./include/myelastic.c')

d = base.create_base()

output_files, input_files, fc = forward(d)

Flow(output_files, input_files, fc)

output_files = re.sub('_synthetic', '', output_files)

zs = 0.2
ze = 0.8
nz = 20
z = np.linspace(zs, ze, nz)

xs = 0.2
xe = 0.8
nx = 20
x = np.linspace(xs, xe, nx)

i = 0
for zz in z:
    for xx in x:
        i += 1
        fc1 = re.sub('sszf=.*? ', 'sszf=%f '%zz, fc)
        fc2 = re.sub('eszf=.*? ', 'eszf=%f '%zz, fc1)
        fc3 = re.sub('ssxf=.*? ', 'ssxf=%f '%xx, fc2)
        new_fc = re.sub('esxf=.*? ', 'esxf=%f '%xx, fc3)
        of1 = re.sub('wavz', 'wavz%d'%i, output_files)
        new_output_files = re.sub('wavx', 'wavx%d'%i, of1)
        Flow(new_output_files, input_files, new_fc)


'''
num_snaps = 200
step_int = int(max(num_snaps, d['nt']) / num_snaps)
snapshot_movie('wavz', 'wavz_synthetic', \
    step_int * np.array(range(d['nt'] / step_int)))
'''
 
proj.End()
