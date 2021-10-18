from rsf.proj import *
from subprocess import check_output as co
import m8r
import matplotlib.pyplot as plt
import os
import numpy as np
from forward import forward
import renormalize as rn
from scipy.optimize import minimize
import copy
from seqflow import *
from itertools import product as cartesian_product
from numpy import inf,nan

global dist_hist
dist_hist = 0

d = {
        'case' : 'synthetic',
        'nz' : 100,
        'nx' : 100,
        'dt' :  5.0e-4,
        'nt' :  100,
        'sszf' :  0.2,
        'ssxf' :  0.5,
        'eszf' :  0.2,
        'esxf' :  0.5,
        'nsz' :  1,
        'nsx' :  1,
        'nb' : 30,
        'fm' : 20.0
}

d['vp'] = '2.0'
d['vs'] = '0.707 * input'
d['rho'] = '1.0'

forward(d)

ls_out = co('ls', shell=True).decode('utf-8')

t = np.linspace(0.0, d['dt']*(d['nt']-1), d['nt'])

SeqFlow('t',None,
        '''math output="x1" n1=%s d1=%.8f'''%(d['nt'], d['dt']))

if( 'wavz_synthetic.rsf' in ls_out \
        and 'wavx_synthetic.rsf' in ls_out ):
        # wavz_syn_input = m8r.Input('wavz_synthetic.rsf')
        # # wavz_syn = wavz_syn_input.read(shape=(d['nz'],d['nx'],d['nt']))
        # wavz_syn = wavz_syn_input.read(shape=(d['nt'], d['nx'], d['nz']))
        # wavz_syn = wavz_syn.transpose((2,1,0))
        # wavx_syn_input = m8r.Input('wavx_synthetic.rsf')
        # # wavx_syn = wavx_syn_input.read(shape=(d['nz'],d['nx'],d['nt']))
        # wavx_syn = wavx_syn_input.read(shape=())

        def misfit(p, s='C'):
                global dist_hist
                z0 = p[0]
                x0 = p[1]
                e = copy.copy(d)
                e['case'] = '%.4f_%.4f'%(z0,x0)
                e['sszf'] = z0
                e['ssxf'] = x0
                e['eszf'] = z0
                e['esxf'] = x0
                forward(e)

                #define tmp variables with strings
                wxs = 'wavx_' + e['case'] + '.rsf'
                wzs = 'wavz_' + e['case'] + '.rsf'
                wxs_ref = 'wavx_synthetic.rsf'
                wzs_ref = 'wavz_synthetic.rsf'

                #compile program and name executable
                prog = Program('ot.c')
                wass = str(prog[0])

                #build commands: note command for z and x is the same
                x_name = 'distx_%s.rsf'%dist_hist
                z_name = 'distz_%s.rsf'%dist_hist
                x_inputs = '%s %s %s'%(wxs, wxs_ref, wass)
                z_inputs = '%s %s %s'%(wzs, wzs_ref, wass)
                cmd = '''./${SOURCES[2]} g=${SOURCES[1]} t=t.rsf'''
                
                SeqFlow(x_name, x_inputs, cmd)
                SeqFlow(z_name, z_inputs, cmd)

                #parse out wasserstein distance
                x_contrib_str = co('sfdisfil < %s'%x_name, shell=True)
                x_contrib_str = x_contrib_str.decode('utf-8')
                z_contrib_str = co('sfdisfil < %s'%z_name, shell=True)
                z_contrib_str = z_contrib_str.decode('utf-8')

                # input(x_contrib_str)
                # input(z_contrib_str)
                #define tmp function to get floating value from string
                get_val = lambda x : float(eval(x.split(':')[-1]))

                #increment dist_hist for next call
                dist_hist += 1

                #get contributions and return them
                # if( 'C' in s ):
                #         return get_val(x_contrib_str) + get_val(z_contrib_str)
                return get_val(x_contrib_str) + get_val(z_contrib_str)

                # input('Why are you executing?')
                # wavx_input = m8r.Input('wavx_' + e['case'] + '.rsf')
                # wavx = wavx_input.read(shape=(d['nt'],d['nz'],d['nx']))
                # wavx = np.transpose(wavx, (2,1,0))
                
                # wavz_input = m8r.Input('wavz_' + e['case'] + '.rsf')
                # wavz = wavz_input.read(shape=(d['nt'],d['nz'],d['nx']))
                # wavz = np.transpose(wavz, (2,1,0))

                # total_sum = 0.0
                # for i in range(len(wavx)):
                #         for j in range(len(wavx[0])):
                #                 total_sum += rn.wass_poly_disc(wavx_syn[i][j], \
                #                         wavx[i][j], t, rn.better_split_normalize) \
                #                         + rn.wass_poly_disc(wavz_syn[i][j], \
                #                                 wavz[i][j], t, rn.better_split_normalize)
                #                 print('(i,j) = (%s,%s,%s)'%(i,j,total_sum))

                # return total_sum[0]

        global callback_i
        global curr_dict
        callback_i = 0

        def my_callback(xk):
                global callback_i
                global curr_dict
                curr_method = curr_dict['method']
                curr_ini = curr_dict['init']
                s = ''.join(80*['*']) + '\nIteration: %d\n'%callback_i
                s = s + 'METHOD: %s\n'%curr_method
                s = s + 'ini=%s, om=%f, TARGET=%s\n'%(curr_ini, misfit(curr_ini), [d['eszf'], d['esxf']])
                s = s + '''xk = %s, misfit=%f\n'''%(xk, misfit(xk))
                s = s + ''.join(80*['*']) + '\n'
                print(s)
                callback_i += 1

        # methods = ['Nelder-Mead']
        # ini = [[0.5, 1.0]]

        # cases = []
        # for meth in methods:
        #         cases.append([])
        #         for ini_curr in ini:
        #                 cases[-1].append({'method': meth, 
        #                         'init': ini_curr,
        #                         'hist': [], 
        #                         'misfit_hist': []})

        # for (m, meth) in enumerate(methods):
        #         for (i, ini_curr) in enumerate(ini):
        #                 curr_dict = cases[m][i]
        #                 minimize(misfit, 
        #                         ini_curr, 
        #                         method=meth, 
        #                         callback=my_callback,
        #                         tol=1e-10)
        #                 cases[m][i] = curr_dict
        #                 callback_i = 0
        nz = 5
        nx = 5
        sz = 0.1
        ez = 0.9
        sx = 0.1
        ex = 0.9
        pz = np.linspace(sz, ez, nz)
        px = np.linspace(sx,ex,nx)
        Z,X = np.meshgrid(pz,px)
        F = np.zeros(Z.shape)
        for i in range(len(pz)):
                for j in range(len(px)):
                        curr = [pz[i], px[j]]
                        F[j][i] = misfit(curr)
                        print('(curr, misfit) = (%s,%s)'%(curr, F[j][i]))

        input(F)

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.contour3D(Z,X,F,cmap='viridis')
        plt.title('Wasserstein Landscape')
        plt.savefig('radial.png')
        #print('m in C = %.15f'%misfit([0.5, 0.5], 'C'))
        #print('m in Python = %f'%misfit([0.5,0.5], 'P'))

else:
        print('Synthetic data not generated yet. Try running "scons" again.')

End()