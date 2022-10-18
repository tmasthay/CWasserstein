import os
from time import time

def clean_files(files, ext='.rsf', check_none=False):
    if( check_none and type(files) == type(None) ):
        return []
    elif( type(files) == str ):
        files = files.split(' ')
    if( type(files) != list ):
        print('FATAL: Invalid file format type for clean_files')
        exit(-1)
    return [e if e[-4:] in ['.rsf', '.vpl', '.exe'] else e + ext for e in files]

def clean_cmd(output_files, input_files, cmd):
    for i in range(len(output_files)):
        #input('TARGETS[%d] -> %s'%(i, output_files[i]))
        cmd = cmd.replace('${TARGETS[%d]}'%i, output_files[i])
    
    for i in range(len(input_files)):
        #input('SOURCES[%d] -> %s'%(i, input_files[i]))
        cmd = cmd.replace('${SOURCES[%d]}'%i, input_files[i])

    y = [e.strip() for e in cmd.split('|')]
    x = [e \
        if e.split(' ')[0][:2] in ['sf', './'] \
        else 'sf' + e for e in y]
        
    if( len(input_files) > 0 ):
        x[0] = '< ' + input_files[0] + ' ' + x[0]
    x[-1] = x[-1] + ' > ' + output_files[0]
    return ' | '.join(x).replace('\n', '')
    
def SeqFlow(output_files, input_files, cmd, verbose=False):
    out_clean = clean_files(output_files)
    inp_clean = clean_files(input_files, check_none=True)
    cmd = clean_cmd(out_clean, inp_clean, cmd)
    if( verbose ):
        print(''''\'\'\'%s\'\'\'...time='''%cmd)
    t = time()
    y = os.system(cmd)
    if( verbose ):
        print(time() - t)

def SeqPlot(output_files, input_files, cmd, verbose=False):
    out_clean = clean_files(output_files, ext='.vpl')
    inp_clean = clean_files(input_files, check_none=True)
    cmd = clean_cmd(out_clean, inp_clean, cmd)
    if( verbose ):
        print(''''\'\'\'%s\'\'\'...time='''%cmd)
    t = time()
    y = os.system(cmd)
    if( verbose ):
        print(time() - t)

def SeqFlowV(output_files, input_files, cmd):
    SeqFlow(output_files, input_files, cmd, verbose=True)

def SeqPlotV(output_files, input_files, cmd):
    SeqPlot(output_files, input_files, cmd, verbose=True)
