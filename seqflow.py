import os

def clean_files(files, default_ext='.rsf', check_none=False):
    if( check_none and type(files) == type(None) ):
        return []
    elif( type(files) == str ):
        files = files.split(' ')
    if( type(files) != list ):
        print('FATAL: Invalid file format type for clean_files')
        exit(-1)
    return [e if '.' in e else e + '.rsf' for e in files]

def clean_cmd(output_files, input_files, cmd):
    for i in range(len(output_files)):
        cmd = cmd.replace('${TARGETS[%d]}'%i, output_files[i])
    
    for i in range(len(input_files)):
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
        print(''''\'\'\'%s\'\'\''''%cmd)
    y = os.system(cmd)
    