import numpy as np
from subprocess import check_output as co
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys
sys.path.append('./include')

from helper import *

plt.rcParams['text.usetex'] = True

def get_title(s, d=dict()):
    if('sobolev' in s):
        val = float(s.split('-')[-1])
        if( val == 0.0 ): return '$L^2$'
        else: return "$H^{%.2f}"%val
    else:
        rescale = s[2:]
        param_string = dict_to_str(d)
        return '$W_2$ with %s normalization\n%s'%(rescale, param_string)


if( len(sys.argv) == 1 ):
    threshold = 1e60
else:
    threshold = float(sys.argv[1])
    print('threshold = %f'%threshold)

x_files = co('ls -t xdist*.rsf | tac', shell=True).decode('utf-8').split()
z_files = co('ls -t zdist*.rsf | tac', shell=True).decode('utf-8').split()

misfits = np.zeros(len(x_files))

mode = eval(co('''grep "^mode=" SConstruct''', shell=True) \
    .decode('utf-8') \
    .replace('\n','') \
    .split('=')[-1])

sobolev_norm = ''
if( mode == 'sobolev' ):
    sobolev_norm = eval(co('''grep "^sobolev_norm[ ]*=" SConstruct''', shell=True) \
        .decode('utf-8') \
        .replace('\n', '') \
        .split('=')[-1])
    mode='%s-%.2f'%(mode, float(sobolev_norm))

for i in range(len(x_files)):
    misfits[i] = get_val(x_files[i]) + get_val(z_files[i])
    misfits[i] = min(misfits[i], threshold)

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

Z,X = np.meshgrid(z,x)
misfits = misfits.reshape(Z.shape)

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(Z,X,misfits)
ax.set_title(r'Convexity plot for %s'%get_title(mode))
ax.set_xlabel('Z')
ax.set_ylabel('X')

fig.savefig('Fig/convexity-%s.png'%mode)

plt.clf()

fig, ax = plt.subplots(1,1)
cb = ax.contourf(Z,X,misfits)
fig.colorbar(cb)
ax.set_title(r'Convexity plot for %s'%get_title(mode))
ax.set_xlabel('Z')
ax.set_ylabel('X')

fig.savefig('Fig/convexity-contour-%s.png'%mode) 

np.save('Fig/misfits.npy', misfits)
np.save('Fig/Z.npy', Z)
np.save('Fig/X.npy', X)
