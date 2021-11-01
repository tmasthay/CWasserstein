#from rsf.proj import *
from seqflow import *

def  forward(d):
    #case name tag
    case_name = d['case']

    #append tag to fields
    def get_field(field_name):
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

    def create_field(fexpr, field_name, input_file=None):
        if( input_file != None ):
            input_file = input_file.replace('.rsf', '') + '.rsf'

        output_file = get_field(field_name).replace('.rsf','')
        output_file += '.rsf'
        SeqFlow(output_file,
            input_file,
            '''
            math output="%s" n1=%s n2=%s d1=%.15f d2=%.15f
            label1=x1 unit1=km label2=x2 unit2=km
            '''%(fexpr, n1, n2, d1, d2))

    def forward_command(nbt, fmt, dtt, ssxft, sszft, esxft, eszft, nxft, nzft,kt):
        s = 'nb=%d fm=%d dt=%.15f nt=%s kt=%s ssxf=%s sszf=%s'%(nbt,fmt,dtt,kt,kt,ssxft,sszft)
        s = s + ' esxf=%s eszf=%s nxf=%s nzf=%s'%(esxft, eszft, nxft, nzft)
        t = './${SOURCES[3]} vs=${SOURCES[1]} rho=${SOURCES[2]} wavx=${TARGETS[1]}'
        return t + ' ' + s

    def combine(arr):
        a = [get_field(s).replace('.rsf','') + '.rsf' for s in arr]
        a = ' '.join(a)
        return a
    
    #prog = Program('myelastic.c')
    #elas = str(prog[0])
    elas = 'myelastic.exe'

    fc = forward_command(nb, fm, d3, ssxf, sszf, esxf, eszf, nsx, nsz, n3)

    create_field(vp, 'vp')
    create_field(vs, 'vs', get_field('vp'))
    create_field(rho, 'rho')
    SeqFlow(combine(['wavz','wavx']), combine(['vp', 'vs', 'rho']) + ' ' + elas,
        fc)

def gauss_test(mu,sig, time_shifts, nz, nx, nt, dz, dx, dt):
    muz = mu[0]
    mux = mu[1]
    mut = mu[2]

    sigz = sig[0]
    sigx = sig[1]
    sigt = sig[2]

    def my_gauss(var_name, mean, stddev):
        C = sqrt(2) * stddev
        C = 1.0 / C
        D = mean * C
        s = '%.15e * %s - %.15e'%(C,var_name,D)
        s = '%s * %s'%(s,s)
        return 'exp(%s)'%s

    w2 = []

    basezx = '%s * %s'%(my_gauss('x1', muz, sigz), my_gauss('x2', mux, sigx))
    tail_cmd = 'n1=%d n2=%d n3=%d d1=%.15e d2=%.15e d3=%.15e'%(nz,nx,nt,dz,dx,dt)

    def get_curr(the_shift):
        return '%s * %s'%(basezx, my_gauss('x3', mut-the_shift, sigt))

    Flow('t_test', None, 'math output="x1" n1=%d d1=%.15e'%(nt, dt))
    Flow('ref_t_test', None, 'math output="%s" %s'%(get_curr(0.0), tail_cmd))

    for t in time_shifts:
        curr = get_curr(t)
        output_name = 'test_%.4e'%t
        Flow(output_name, None, 'math output="%s" %s'%(curr,tail_cmd))

    