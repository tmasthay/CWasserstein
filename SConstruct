from rsf.proj import *
from subprocess import check_output as co
import m8r
import matplotlib.pyplot as plt
import os
import numpy as np
from forward import forward
from scipy.optimize import minimize
from seqflow import *
from itertools import product as cartesian_product
from numpy import inf,nan
from time import time
from tmp import *
from copy import *
from gradients import *

global dist_hist
dist_hist = 0

global d
d = {
        'case' : 'synthetic',
        'nz' : 100,
        'dz' : 0.01,
        'nx' : 100,
        'dx' : 0.01,
        'dt' :  5e-03,
        'nt' :  500,
        'sszf' :  0.5,
        'ssxf' :  0.5,
        'eszf' :  0.5,
        'esxf' :  0.5,
        'nsz' :  1,
        'nsx' :  1,
        'nb' : 30,
        'fm' : 5.0
}

d['vp'] = '0.5'
d['vs'] = '0.707 * input'
d['rho'] = '1.0'

forward(d)

ls_out = co('ls', shell=True).decode('utf-8')

t = np.linspace(0.0, d['dt']*(d['nt']-1), d['nt'])

SeqFlow('t',None,
        '''math output="x1" n1=%s d1=%.8f'''%(d['nt'], d['dt']))

SeqFlow('x',None,
        '''math output="x1" n1=%s d1=%.8f'''%(d['nx'], d['dx']))

if( 'wavz_synthetic.rsf' in ls_out \
        and 'wavx_synthetic.rsf' in ls_out ):
        # wavz_syn_input = m8r.Input('wavz_synthetic.rsf')
        # # wavz_syn = wavz_syn_input.read(shape=(d['nz'],d['nx'],d['nt']))
        # wavz_syn = wavz_syn_input.read(shape=(d['nt'], d['nx'], d['nz']))
        # wavz_syn = wavz_syn.transpose((2,1,0))
        # wavx_syn_input = m8r.Input('wavx_synthetic.rsf')
        # # wavx_syn = wavx_syn_input.read(shape=(d['nz'],d['nx'],d['nt']))
        # wavx_syn = wavx_syn_input.read(shape=())

        def misfit(p, s='W2'):
                misfit_exec_time_start = time()
                global dist_hist
                global d
                t_big = time()
                z0 = p[0]
                x0 = p[1]
                q = deepcopy(p)
                e = deepcopy(d)
                e['case'] = '%.4f_%.4f'%(z0,x0)
                e['sszf'] = z0
                e['ssxf'] = x0
                e['eszf'] = z0
                e['esxf'] = x0

                if( not os.path.exists('wavz_' + e['case'] + '.rsf') ):
                        forward(e)

                #define tmp variables with strings
                wxs = 'wavx_' + e['case'] + '.rsf'
                wzs = 'wavz_' + e['case'] + '.rsf'
                wxs_ref = 'wavx_synthetic.rsf'
                wzs_ref = 'wavz_synthetic.rsf'

                #compile program and name executable
                prog = Program('ot.c')
                wass = str(prog[0])

                if( 'W' in s ):
                        #build commands: note command for z and x is the same
                        x_name = 'distx_%s.rsf'%dist_hist
                        z_name = 'distz_%s.rsf'%dist_hist
                        x_inputs = '%s %s %s'%(wxs, wxs_ref, wass)
                        z_inputs = '%s %s %s'%(wzs, wzs_ref, wass)
                        cmd = '''./${SOURCES[2]} g=${SOURCES[1]} t=t.rsf'''
                
                        SeqFlow(x_name, x_inputs, cmd)
                        SeqFlow(z_name, z_inputs, cmd)
                        
                        SeqFlow('trans_' + wxs, wxs, 'sftransp plane=23')
                        SeqFlow('trans_' + wxs_ref, wxs_ref, 'sftransp plane=23')
                        SeqFlow('trans_' + wzs, wzs, 'sftransp plane=23')
                        SeqFlow('trans_' + wzs_ref, wzs_ref, 'sftransp plane=23')

                        x_inputs_trans = '%s %s %s'%('trans_' + wxs, 'trans_' + wxs_ref, wass)
                        z_inputs_trans = '%s %s %s'%('trans_' + wzs, 'trans_' + wzs_ref, wass)
                        cmd2 = cmd.replace('t=t.rsf', 't=x.rsf')
                        SeqFlow('trans_' + x_name, x_inputs_trans, cmd2)
                        SeqFlow('trans_' + z_name, z_inputs_trans, cmd2)

                        #parse out wasserstein distance
                        x_contrib_str = co('sfdisfil < %s'%x_name, shell=True)
                        x_contrib_str = x_contrib_str.decode('utf-8')
                        z_contrib_str = co('sfdisfil < %s'%z_name, shell=True)
                        z_contrib_str = z_contrib_str.decode('utf-8')

                        xct = co('sfdisfil < %s'%('trans_' + x_name), shell=True).decode('utf-8')
                        zct = co('sfdisfil < %s'%('trans_' + z_name), shell=True).decode('utf-8')

                        # input(x_contrib_str)
                        # input(z_contrib_str)
                        #define tmp function to get floating value from string
                        #get_val = lambda x : float(eval(x.split(':')[-1]))
                        def get_val(x):
                                #print('x=%s'%x)
                                return float(eval(x.split(':')[-1]))

                        #increment dist_hist for next call
                        dist_hist += 1

                        #get contributions and return them
                        # if( 'C' in s ):
                        #         return get_val(x_contrib_str) + get_val(z_contrib_str)
                        #return get_val(xct) + get_val(zct)
                        print('%s\nMISFIT EXEC TIME: %f'%(''.join(80*['*']), time() - misfit_exec_time_start))
                        option = 'time'
                        if( option == 'both' ):
                                return get_val(x_contrib_str) + get_val(z_contrib_str) + get_val(xct) + get_val(zct)
                        elif( option == 'time' ):
                                return get_val(x_contrib_str) + get_val(z_contrib_str)
                        else:
                                return get_val(xct) + get_val(zct)
                else:
                        wxss = wxs.replace('.rsf','')
                        wzss = wzs.replace('.rsf','')
                        wxss_ref = wxs_ref.replace('.rsf','')
                        wzss_ref = wzs_ref.replace('.rsf','')
                        get_top = lambda x : SeqFlow('%s_top.rsf'%x, x + '.rsf', 'window n1=1 f1=0 | put o3=0', True)
                        extract = lambda x : float([y for y in x.split('\n') if 'norm' in y][0].split('=')[1])
                        get_top(wxss)
                        get_top(wzss)
                        get_top(wxss_ref)
                        get_top(wzss_ref)
                        SeqFlow('tmpx', '%s_top.rsf %s_top.rsf'%(wxss, wxss_ref), 'math output="input-y" y=${SOURCES[1]}')
                        SeqFlow('tmpz', '%s_top.rsf %s_top.rsf'%(wzss, wzss_ref), 'math output="input-y" y=${SOURCES[1]}')
                        sx = extract(co('sfattr < tmpx.rsf', shell=True).decode('utf-8'))
                        sz = extract(co('sfattr < tmpz.rsf', shell=True).decode('utf-8'))
                        return sx**2 + sz**2

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
        
        actual = True
        if( actual ):
                t_big = time()
                global callback_i
                global curr_dict
                global point_guesses
                global misfit_guesses
                callback_i = 0

                point_guesses = []
                misfit_guesses = []

                def my_callback(xk):
                        start_time = time()
                        global callback_i
                        global curr_dict
                        global point_guesses
                        global misfit_guesses
                        curr_method = curr_dict['method']
                        curr_ini = curr_dict['init']
                        #point_guesses.append(xk)
                        #misfit_guesses.append(misfit(xk))
                        s = 'Iteration: %d\n'%callback_i
                        #s = s + 'METHOD: %s\n'%curr_method
                        #s = s + 'ini=%s, om=%f, TARGET=%s\n'%(curr_ini, misfit(curr_ini), [d['eszf'], d['esxf']])
                        #s = s + '''xk = %s, misfit=%.15f\n'''%(xk, misfit_guesses[-1])
                        s = s + 'xk=%s\n'%xk
                        s = s + 'CALLBACK_TIME=%f\n'%(time() - start_time)
                        s = s + ''.join(80*['*']) + '\n'
                        print(s)
                        callback_i += 1


                z_ini = np.linspace(0.3,0.7,5)
                x_ini = np.linspace(0.3,0.7,5)
                
                ini = [np.array(guess) for guess in cartesian_product(z_ini, x_ini)]
                deriv_info = [[lambda x: numerical_deriv(x,misfit,0.05), lambda x: numerical_hessian(x,misfit,0.05)]]

                methods = ['Nelder-Mead']

                final_answers = []

                point_guesses.append(ini[0])
                misfit_guesses.append(misfit(ini[0]))
                cases = []
                y = 1.0
                for meth in methods:
                        cases.append([])
                        for ini_curr in ini:
                                cases[-1].append({'method': meth, 
                                        'init': ini_curr,
                                        'hist': [], 
                                        'misfit_hist': []})

                last_time = time()
                for (m, meth) in enumerate(methods):
                        for (i, ini_curr) in enumerate(ini):
                                curr_dict = cases[m][i]
                                money_signs = ''.join(80 * ['$'])
                                s = '%s\nSTARTING ITERATION: (%d, %d) of (%d, %d)\n'%(money_signs, m,i,len(methods), len(ini))
                                s = s + 'LAST ITERATION TOOK %f seconds.\n%s'%(time() - last_time, money_signs)
                                print(s)
                                last_time = time()
                                y = minimize(misfit, 
                                        ini_curr, 
                                        method=meth, 
                                        callback=my_callback,
                                        tol=0.0,
                                        jac=deriv_info[i][0],
                                        hess=deriv_info[i][1],
                                        options={'maxiter': 20})
                                cases[m][i] = curr_dict
                                final_answers[i] = y['x']
                                callback_i = 0

                print('Total time: %.4f'%(time() - t_big))
                print('y = %s'%y)
                #print('m in C = %.15f'%misfit([0.1, 0.9], 'C'))
                #print('m in Python = %f'%misfit([0.5,0.5], 'P'))

                xx = [p[0] for p in point_guesses]
                yy = [p[1] for p in point_guesses]
                blue_channel = np.linspace(0.0, 1.0, len(xx))
                red_channel = np.linspace(1.0, 0.0, len(xx))
                C = [[red_channel[i], 0.0, blue_channel[i]] for i in range(len(blue_channel))]
                print(C)

                ref_pt = [d['eszf'], d['esxf']]
                far_off = [np.linalg.norm(ref_pt - guess) for guess in final_answers]
                
                Z,X = np.meshgrid(z_ini,x_ini)

                fig3 = plt.figure()
                ax = plt.axes(projection='3d')
                ax.plot_surface(Z,X,far_off.reshape((len(z_ini),len(x_ini))),cmap='viridis')
                plt.title('Error versus Initial Guess')
                plt.xlabel('X_initial')
                plt.ylabel('Y_initial')
                plt.zlabel('Final Error')
                plt.savefig('guesses.png')
                plt.clf()

                # plt.scatter(xx,yy, color=C)
                # plt.scatter(0.5,0.5,color=[0.0,1.0,0.0])
                # plt.title('Surface-only W2  (red first guess, blue last guess, green ideal)')
                # plt.savefig('scatter.png')
                # plt.clf()
                # plt.plot(np.array(range(len(point_guesses))), misfit_guesses)
                # plt.title('Surface-only W2')
                # plt.ylabel('Misfit')
                # plt.xlabel('Iteration Number')
                # plt.savefig('misfit.png')

                # plt.clf()
                # plt.figure()
        else:
                nz = 25
                nx = 25
                sz = 0.3
                ez = 0.7 
                sx = 0.3
                ex = 0.7
                pz = np.linspace(sz, ez, nz)
                px = np.linspace(sx,ex,nx)
                Z,X = np.meshgrid(pz,px)
                F = np.zeros(Z.shape)
                G = np.zeros(Z.shape)
                for i in range(len(pz)):
                        for j in range(len(px)):
                                t = time()
                                curr = [pz[i], px[j]]
                                # F[j][i] = misfit(curr)
                                #F[j][i] = misfit(curr, 'L2')
                                print('(j,i,x) = (%d,%d,%s)'%(j,i,curr), end='')
                                G[j][i] = misfit(curr, 'W2')
                                print('...misfit computed')
                                t = time() - t
                                print('(curr, L2misfit, W2, exec time) = (%s,%s,%s,%.4f)'%(curr, F[j][i], G[j][i], t))

                G = [[x if x not in [inf, nan] else 200.0 for x in e] for e in G]
                G = [[x if x <= 199.0 else 200.0 for x in e] for e in G]
                G = np.array(G)
                #input(F)

                # fig = plt.figure()
                # ax = plt.axes(projection='3d')
                # ax.contour3D(Z,X,F,cmap='viridis')
                # plt.title('L2 landscape')
                # plt.savefig('L2landscape-global-contour.png')

                # fig2 = plt.figure()
                # ax = plt.axes(projection='3d')
                # ax.contour3D(Z,X,G,cmap='viridis')
                # plt.title('W2 landscape')
                # plt.savefig('W2landscape-surfaceOnly-contour.png')

                fig3 = plt.figure()
                ax = plt.axes(projection='3d')
                ax.plot_surface(Z,X,F,cmap='viridis')
                plt.title('L2 landscape')
                plt.savefig('L2landscape-global.png')

                fig4 = plt.figure()
                ax = plt.axes(projection='3d')
                ax.plot_surface(Z,X,G,cmap='viridis')
                plt.title('W2 landscape')
                plt.savefig('W2landscape-surfaceOnly.png')
else:
        print('Synthetic data not generated yet. Try running "scons" again.')

#windowify(d, [11])
#print(y)

End()