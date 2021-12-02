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
from copy import *
from gradients import *
from plotall import plot_all

def SeqFlowLocal(output_files, input_files, cmd):
    verbose = False
    SeqFlow(output_files, input_files, cmd, verbose)

def get_val(s):
    try:
        return float(s.split(':')[-1].replace('\n','').replace(' ',''))
    except:
        input('get_val error: Press enter to dump s before exiting')
        print(s)
        exit(-1)

def get_val_full(s):
    try:
        t = co('sfdisfil < %s'%s, shell=True).decode('utf-8')
        return float(t.split(':')[-1].replace('\n','').replace(' ',''))
    except:
        input('get_val_full error: Press enter to dump s before exiting')
        print(t)
        exit(-1)

global dist_hist
dist_hist = 0

global wass_hist
wass_hist = []

global d
d = {
    'case' : 'synthetic',
    'nz' : 200,
    'dz' : 0.005,
    'nx' : 200,
    'dx' : 0.005,
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

SeqFlowLocal('t',None,
    '''math output="x1" n1=%s d1=%.8f'''%(d['nt'], d['dt']))

SeqFlowLocal('x',None,
    '''math output="x1" n1=%s d1=%.8f'''%(d['nx'], d['dx']))

os.system('cp t.rsf tmpt.txt')
os.system('cp x.rsf tmpx.txt')
os.system('cp /var/tmp/t.rsf@ /var/tmp/tmpt.txt')
os.system('cp /var/tmp/x.rsf@ /var/tmp/tmpx.txt')

if( 'wavz_synthetic.rsf' in ls_out \
    and 'wavx_synthetic.rsf' in ls_out ):

    def misfit(p, s='W2'):
        if( min(p) <= 0.0 or max(p) >= 1.0 ):
            return 1e10
        misfit_exec_time_start = time()
        global dist_hist
        global d
        global wass_hist
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
        #prog = Program('sw.c')
        prog = Program('ot.c')
        wass = str(prog[0])

        #cmd = './${SOURCES[2]} g=${SOURCES[1]} t=t.rsf x=x.rsf'
        cmd = './${SOURCES[2]} g=${SOURCES[1]} t=t.rsf'

        SeqFlowLocal('distx.rsf', '%s %s %s'%(wxs_ref, wxs, wass), cmd)
        SeqFlowLocal('distz.rsf', '%s %s %s'%(wzs_ref, wzs, wass), cmd)

        return get_val_full('distx.rsf') + get_val_full('distz.rsf')
 
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
            point_guesses.append(xk)
            misfit_guesses.append(misfit(xk))
            s = 'Iteration: %d\n'%callback_i
            #s = s + 'METHOD: %s\n'%curr_method
            #s = s + 'ini=%s, om=%f, TARGET=%s\n'%(curr_ini, misfit(curr_ini), [d['eszf'], d['esxf']])
            s = s + '''xk = %s, misfit=%.15f\n'''%(xk, misfit_guesses[-1])
            s = s + 'xk=%s\n'%xk
            s = s + 'CALLBACK_TIME=%f\n'%(time() - start_time)
            s = s + ''.join(80*['*']) + '\n'
            print(s)
            callback_i += 1


        nz_ini = 15
        nx_ini = 15
        num_iter = 20
        tau = 0.0
        z_ini = np.linspace(0.3,0.7,nz_ini)
        x_ini = np.linspace(0.3,0.7,nx_ini)
        
        ini = [np.array(guess) for guess in cartesian_product(z_ini, x_ini)]
        converge_error = False
        if( not converge_error ):
            ini = [[0.6, 0.6]]

        deriv_info = len(ini) * ([[lambda x: numerical_deriv(x,misfit,0.1), lambda x: numerical_hessian(x,misfit,0.1)]])
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
                curr_res = minimize(misfit, 
                    ini_curr, 
                    method=meth, 
                    callback=my_callback,
                    tol=tau,
                    jac=deriv_info[i][0],
                    hess=deriv_info[i][1],
                    options={'maxiter': num_iter})
                cases[m][i] = curr_dict
                final_answers.append(curr_res['x'])
                callback_i = 0
                os.system('''rm $(ls *.rsf | grep -v "synthetic")''')
                os.system('''rm $(ls /var/tmp/[a-zA-Z][a-zA-Z]*.rsf@ | grep -v "synthetic")''')
                os.system('''rm $(ls /var/tmp/[v]*.rsf@ | grep -v "synthetic")''')
                os.system('cp tmpt.txt t.rsf')
                os.system('cp tmpx.txt x.rsf')
                os.system('cp /var/tmp/tmpt.txt /var/tmp/t.rsf@')
                os.system('cp /var/tmp/tmpx.txt /var/tmp/x.rsf@')

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

        plot_all(wass_hist)
        os.system('rm *.rsf')

        if( converge_error ):
            print('yo')
            # ref_pt = [d['eszf'], d['esxf']]
            # far_off = [np.linalg.norm(ref_pt - guess) for guess in final_answers]
            # far_off = np.array(far_off)
            
            # Z,X = np.meshgrid(z_ini,x_ini)

            # fig3 = plt.figure()
            # ax = plt.axes(projection='3d')
            # ax.plot_surface(Z,X,far_off.reshape((len(z_ini),len(x_ini))),cmap='viridis')
            # plt.title('Error versus Initial Guess: L2 surface')
            # plt.xlabel('X_initial')
            # plt.ylabel('Y_initial')
            # #plt.zlabel('Final Error')
            # plt.savefig('guesses.png')
            # plt.clf()

        else:
            print('Creating scatter plot')
            plt.scatter(xx,yy, color=C)
            plt.scatter(0.5,0.5,color=[0.0,1.0,0.0])
            plt.title('sw  (red first guess, blue last guess, green ideal)')
            plt.savefig('scatter.png')
            plt.clf()
            plt.plot(np.array(range(len(point_guesses))), misfit_guesses)
            plt.title('Surface-only SW2')
            plt.ylabel('Misfit')
            plt.xlabel('Iteration Number')
            plt.savefig('misfit.png')

        # plt.clf()
        # plt.figure()
else:
    print('Synthetic data not generated yet. Try running "scons" again.')

#windowify(d, [11])
#print(y)

End()
