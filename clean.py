from seqflow import *
from subprocess import check_output as co
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from h1 import *


Flow=SeqFlowV
Plot=SeqPlotV

def transp_data(dest, src):
    Flow(dest, src,
        '''
        transp plane=13 | transp plane=23
        ''')

def get_val(s):
    try:
        return float(s.split(':')[-1].replace('\n','').replace(' ',''))
    except:
        return 1e10

def gv_cmd(s):
    return get_val(co(s, shell=True).decode('utf-8'))

mode=sys.argv[1]

nz=50
nx=50

wavz_files = co('ls wavz_*_*_top.rsf', shell=True).decode('utf-8').split()
wavx_files = co('ls wavx_*_*_top.rsf', shell=True).decode('utf-8').split()

z_ref = 'wavz_synthetic_top.rsf'
x_ref = 'wavx_synthetic_top.rsf'

transp_data('z_ref', z_ref)
transp_data('x_ref', x_ref)

z_ref = 'z_ref'
x_ref = 'x_ref'

cmd = './ot.exe mode=3 t=t.rsf g=${SOURCES[1]}'

z = np.zeros(len(wavz_files))
x = np.zeros(len(wavx_files))
misfits = np.zeros(len(wavz_files))

for i in range(len(wavz_files)):
    transp_data(wavz_files[i].replace('.rsf', 'trans.rsf'), wavz_files[i])
    transp_data(wavx_files[i].replace('.rsf', 'trans.rsf'), wavx_files[i])

    pt = [float(ss) for ss in wavz_files[i].split('_')[1:3]]
    z[i] = pt[0]
    x[i] = pt[1]

    if( mode[0] == '1' ):
        Flow('diff1',
            '%s %s'%(z_ref, wavz_files[i].replace('.rsf','trans.rsf'))
            , cmd)
        Flow('diff2', 
            '%s %s'%(x_ref, wavx_files[i].replace('.rsf','trans.rsf'))
            , cmd)
        misfits[i] = gv_cmd('sfdisfil < diff1.rsf') \
            + gv_cmd('sfdisfil < diff2.rsf')
    elif( mode[0] == '3' ):
        cur_z = wavz_files[i].replace('_top','')
        cur_x = wavx_files[i].replace('_top','')
        ref_z = 'wavz_synthetic'
        ref_x = 'wavx_synthetic'
        Flow('diffz', '%s  %s'%(cur_z, ref_z), cmd)
        Flow('diffx', '%s %s'%(cur_x, ref_x), cmd)
        misfits[i] = gv_cmd('sfdisfil < diffz.rsf') + \
            gv_cmd('sfdisfil < diffx.rsf') 
    elif( mode[0] == 'h' ):
        misfits[i] = hs('wavz_synthetic_top', 
            wavz_files[i].replace('.rsf',''), 
            1,
            1)
        misfits[i] += hs('wavx_synthetic_top', 
            wavx_files[i].replace('.rsf',''), 
            1,
            1)
    else:
        misfits[i] = l2diff('wavz_synthetic_top', wavz_files[i])
        misfits[i] += l2diff('wavx_synthetic_top', wavx_files[i])

version = co('ls [0-9]*.png', shell=True).decode('utf-8').split()
version_no = np.max([int(e.split('.')[0]) for e in version]) + 1

Z = z.reshape((nz,nx))
X = x.reshape((nz,nx))
misfits = misfits.reshape((nz,nx))


np.save('Z-%d.npy'%version_no, Z)
np.save('X-%d.npy'%version_no, X)
np.save('misfits-%d.npy'%version_no, misfits)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(Z,X,misfits)
plt.savefig('%d.png'%version_no)

os.system('mkdir case-%d'%version_no)
os.system('mv X-%d.png case-%d'%(version_no, version_no))
os.system('mv Z-%d.png case-%d'%(version_no, version_no))
os.system('mv misfits-%d.png case-%d'%(version_no, version_no))
os.system('mv %d.png case-%d'%(version_no, version_no))

    



