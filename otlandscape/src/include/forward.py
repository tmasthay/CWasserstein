from rsf.proj import *
from subprocess import check_output as co

def  forward(d):
    #case name tag
    case_name = d['case']

    #append tag to fields
    def get_field(field_name, ref=False):
        if(ref):
            return field_name + '_synthetic'
        else:
            return field_name + '_' + case_name

    #discretization info
    n1 = d['nz']
    n2 = d['nx']
    n3 = d['nt']
    d1 = d['dz']
    d2 = d['dx']
    d3 = d['dt']

    #elastic physics
    vp = d['vp']
    vs = d['vs']
    rho = d['rho']

    #source info
    #"start source fractional z,x"
    sszf = d['sszf']
    ssxf = d['ssxf']

    #"end source fractional z,x"
    eszf = d['eszf']
    esxf = d['esxf']

    #"number of sources z,x"
    nsx = d['nsx']
    nsz = d['nsz']

    #PML size, number of ghost boundary cells
    nb = d['nb']

    #Ricker frequency
    fm = d['fm']

    #amplitude
    amp = d['amp']

    #noise level
    noise = d['noise']

    #verbosity
    verb = d['verb']

    delta_x = (d['esxf_val'] - d['ssxf_val']) / d['numx_comp']
    delta_z = (d['eszf_val'] - d['sszf_val']) / d['numz_comp']

    Flow(d['sszf'], None,
        '''
        math output="x1" n1=%d n2=%d
            o1=%.4f o2=%.4f
            d1=%.4f d2=%.4f
        '''%(d['numz_comp'], d['numx_comp'],
            d['sszf_val'], d['ssxf_val'],
            delta_z, delta_x))

    Flow(d['eszf'], d['sszf'], 'math output="input"')

    print('''
        math output="x2" n1=%d n2=%d
            o1=%.4f o2=%.4f
            d1=%.4f d2=%.4f
        '''%(d['numz_comp'], d['numx_comp'],
            d['sszf_val'], d['ssxf_val'],
            delta_z, delta_x))
    Flow(d['ssxf'], None,
        '''
        math output="x2" n1=%d n2=%d
            o1=%.4f o2=%.4f
            d1=%.4f d2=%.4f
        '''%(d['numz_comp'], d['numx_comp'],
            d['sszf_val'], d['ssxf_val'],
            delta_z, delta_x))
    Flow(d['esxf'], d['ssxf'], 'math output="input"')
 
    def forward_command(nbt, fmt, dtt, ssxft, sszft, esxft, eszft, nxft, nzft,kt, verbb):
        s = 'nb=%d fm=%d dt=%.15f nt=%s kt=%s'%(nbt,fmt,dtt,kt,kt)
        s = s + ' nxf=%s nzf=%s'%(nxft, nzft)
        s = s + ' amp=%.4f verb=%d'%(amp, verbb)
        t = './${SOURCES[7]} vs=${SOURCES[1]} rho=${SOURCES[2]} wavx=${TARGETS[1]}'
        t = t + ' sszf=${SOURCES[3]} eszf=${SOURCES[4]}'
        t = t + ' ssxf=${SOURCES[5]} esxf=${SOURCES[6]}'
        return t + ' ' + s

    fc = forward_command(nb, fm, d3, ssxf, sszf, esxf, eszf, nsx, nsz, n3, verb)

    #elas_exe = './include/elastictop.exe' if type(d['sszf']) == type(1.0) \
    #    else './include/elastictopstack.exe'
    elas_exe = './include/elastictopstack.exe'

#    Flow('%s %s'%(get_field('wavz'), get_field('wavx')), \
#         '%s %s %s %s'%(vp, vs, rho, elas_exe), \
#         fc)
    output_files = '%s %s'%(get_field('wavz'), get_field('wavx'))
    input_files = '%s %s %s'%(vp, vs, rho)
    input_files = input_files + ' %s %s %s %s'%(sszf, eszf, ssxf, esxf)
    input_files = input_files + ' %s'%elas_exe
    return output_files, input_files, fc
    
