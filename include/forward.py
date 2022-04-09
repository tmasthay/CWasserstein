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
 
    def forward_command(nbt, fmt, dtt, ssxft, sszft, esxft, eszft, nxft, nzft,kt):
        s = 'nb=%d fm=%d dt=%.15f nt=%s kt=%s ssxf=%s sszf=%s'%(nbt,fmt,dtt,kt,kt,ssxft,sszft)
        s = s + ' esxf=%s eszf=%s nxf=%s nzf=%s'%(esxft, eszft, nxft, nzft)
        s = s + ' amp=%.4f'%amp
        t = './${SOURCES[3]} vs=${SOURCES[1]} rho=${SOURCES[2]} wavx=${TARGETS[1]}'
        return t + ' ' + s
 
    path = 'include/'
    src = path + 'myelastic.c'
    elas_exe = src.replace('.c', '.exe')

    fc = forward_command(nb, fm, d3, ssxf, sszf, esxf, eszf, nsx, nsz, n3)

    Flow('%s %s'%(get_field('wavz'), get_field('wavx')), \
         '%s %s %s %s'%(vp, vs, rho, elas_exe), \
         fc)
