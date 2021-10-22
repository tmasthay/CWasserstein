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
    d1 = 1.0 / n1
    d2 = 1.0 / n2
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
        fc, True)

    