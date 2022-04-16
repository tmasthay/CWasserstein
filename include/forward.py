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
 
    def forward_command(nbt, fmt, dtt, ssxft, sszft, esxft, eszft, nxft, nzft,kt, verbb):
        s = 'nb=%d fm=%d dt=%.15f nt=%s kt=%s ssxf=%s sszf=%s'%(nbt,fmt,dtt,kt,kt,ssxft,sszft)
        s = s + ' esxf=%s eszf=%s nxf=%s nzf=%s'%(esxft, eszft, nxft, nzft)
        s = s + ' amp=%.4f verb=%d'%(amp, verbb)
        t = './${SOURCES[3]} vs=${SOURCES[1]} rho=${SOURCES[2]} wavx=${TARGETS[1]}'
        return t + ' ' + s

    fc = forward_command(nb, fm, d3, ssxf, sszf, esxf, eszf, nsx, nsz, n3, verb)

    elas_exe = './include/elastic_top.exe'

#    Flow('%s %s'%(get_field('wavz'), get_field('wavx')), \
#         '%s %s %s %s'%(vp, vs, rho, elas_exe), \
#         fc)
    output_files = '%s %s'%(get_field('wavz'), get_field('wavx'))
    input_files = '%s %s %s %s'%(vp, vs, rho, elas_exe)
    return output_files, input_files, fc
    
