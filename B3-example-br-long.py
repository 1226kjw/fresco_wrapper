from io import IncrementalNewlineDecoder
from wfresco import *

f = Wfresco('B3-example-br-long.frin', 'alpha+c12 -> alpha+c12* @ 100 MeV; nuc def')

f.set_parameters(
hcm=0.01,
    rmatch=-60.000,
    rintp=0.15,
    rsp=0.0,
    rasym=1000.00,
    accrcy=0.0010000,
    jtmin=0.0,
    jtmax=9000.0,
    absend=-50.0000,
    jump=[1,10,50,200],
    jbord=[0.0,200.0,300.0,1000.0,9000.0],
    theta_range=[0, 20, 0.05],
    cutr=-20.00,
    ips=0.0000,
    it0=0,
    iter=0,
    iblock=21,
    nnu=24,
    smallchan=1.00E-12,
    smallcoup=1.00E-12,
    chans=1,
    smats=2,
    xstabl=1,
    cdcc=1,
    elab=656.0000,
    pel=1,
    exl=1,
    lab=1,
    lin=1,
    lex=1
)

f.set_partition(
    qval=0.1370, nex=21, pwf=True,
    proj='8B', target='208Pb',
    states=[
        {
            'proj':["1.5-", 0],
            'target':["0.0+", 0],
            'cpot':1
        },
        {'proj':["0.5+", 0.1583], 'target':1, 'cpot':1},
        {'proj':["0.5+", 0.2180], 'target':1, 'cpot':1},
        {'proj':["0.5+", 0.3260], 'target':1, 'cpot':1},
        {'proj':["0.5+", 0.4830], 'target':1, 'cpot':1},
        {'proj':["0.5+", 0.6889], 'target':1, 'cpot':1},
        {'proj':["0.5+", 0.9438], 'target':1, 'cpot':1},
        {'proj':["0.5+", 1.2478], 'target':1, 'cpot':1},
        {'proj':["0.5+", 1.6007], 'target':1, 'cpot':1},
        {'proj':["0.5+", 2.0027], 'target':1, 'cpot':1},
        {'proj':["0.5+", 2.4536], 'target':1, 'cpot':1},
        {'proj':["0.5+", 2.9536], 'target':1, 'cpot':1},
        {'proj':["0.5+", 3.5025], 'target':1, 'cpot':1},
        {'proj':["0.5+", 4.1005], 'target':1, 'cpot':1},
        {'proj':["0.5+", 4.7474], 'target':1, 'cpot':1},
        {'proj':["0.5+", 5.4434], 'target':1, 'cpot':1},
        {'proj':["0.5+", 6.1884], 'target':1, 'cpot':1},
        {'proj':["0.5+", 6.9824], 'target':1, 'cpot':1},
        {'proj':["0.5+", 7.8253], 'target':1, 'cpot':1},
        {'proj':["0.5+", 8.7173], 'target':1, 'cpot':1},
        {'proj':["0.5+", 9.6583], 'target':1, 'cpot':1},
    ]
)

f.set_partition(
    qval=0.0000, nex=-1, pwf=True,
    proj='7Be', target='208Pb + p', masst=209, zt=83,
    states=[
        {
            'proj':['0.0+',0.0],
            'target':['0.0+',0.0],
            'cpot':2
        }
    ]
)

f.set_pot(1,
    pots=[
        {
            'type':0,
            'shape':0,
            'p':[1, 0, 2.65]
        },
    ]
)
f.set_pot(2,
    pots=[
        {
            'type':0,
            'shape':0,
            'p':[208, 0, 1.3]
        },
        {
            'type':1,
            'shape':0,
            'p':[114.2000, 1.2860, 0.8530, 9.4400, 1.7390, 0.8090, 0.0000]
        }
    ]
)
f.set_pot(3,
    pots=[
        {
            'type':0,
            'shape':0,
            'p':[208, 0, 1.3]
        },
        {
            'type':1,
            'shape':0,
            'p':[34.8190, 1.1700, 0.7500, 15.3400, 1.3200, 0.6010, 0.0000]
        }
    ]
)
f.set_pot(4,
    pots=[
        {
            'type':0,
            'shape':0,
            'p':[1, 0, 2.3910]
        },
        {
            'type':1,
            'shape':0,
            'p':[44.6750, 2.3910, 0.4800]
        },
        {
            'type':3,
            'shape':0,
            'p':[4.8980, 2.3910, 0.4800]
        }
    ]
)

f.set_overlap(
    kn1=1, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=1, l=1, lmax=0, sn=0.5, j=1.5, nam=1, ampl=1.0000,
    kbpot=4, be=0.1370, isc=1, ipc=0,
)
f.set_overlap(
    kn1=2, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-0.0182, isc=12, ipc=2, nk=20, er=-0.0344
)
f.set_overlap(
    kn1=3, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-0.0771, isc=12, ipc=2, nk=20, er=-0.0834
)
f.set_overlap(
    kn1=4, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-0.1850, isc=12, ipc=2, nk=20, er=-0.1324
)
f.set_overlap(
    kn1=5, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-0.3419, isc=12, ipc=2, nk=20, er=-0.1814
)
f.set_overlap(
    kn1=6, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-0.5479, isc=12, ipc=2, nk=20, er=-0.2304
)
f.set_overlap(
    kn1=7, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-0.8028, isc=12, ipc=2, nk=20, er=-0.2794
)
f.set_overlap(
    kn1=8, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-1.1067, isc=12, ipc=2, nk=20, er=-0.3284
)
f.set_overlap(
    kn1=9, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-1.4596, isc=12, ipc=2, nk=20, er=-0.3774
)
f.set_overlap(
    kn1=10, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-1.8616, isc=12, ipc=2, nk=20, er=-0.4264
)
f.set_overlap(
    kn1=11, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-2.3125, isc=12, ipc=2, nk=20, er=-0.4754
)
f.set_overlap(
    kn1=12, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-2.8125, isc=12, ipc=2, nk=20, er=-0.5245
)
f.set_overlap(
    kn1=13, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-3.3614, isc=12, ipc=2, nk=20, er=-0.5735
)
f.set_overlap(
    kn1=14, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-3.9594, isc=12, ipc=2, nk=20, er=-0.6225
)
f.set_overlap(
    kn1=15, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-4.6064, isc=12, ipc=2, nk=20, er=-0.6715
)
f.set_overlap(
    kn1=16, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-5.3023, isc=12, ipc=2, nk=20, er=-0.7205
)
f.set_overlap(
    kn1=17, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-6.0473, isc=12, ipc=2, nk=20, er=-0.7695
)
f.set_overlap(
    kn1=18, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-6.8413, isc=12, ipc=2, nk=20, er=-0.8185
)
f.set_overlap(
    kn1=19, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-7.6843, isc=12, ipc=2, nk=20, er=-0.8675
)
f.set_overlap(
    kn1=20, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-8.5763, isc=12, ipc=2, nk=20, er=-0.9165
)
f.set_overlap(
    kn1=21, kn2=0, ic1=1, ic2=2, _in=1,
    kind=0, nn=0, l=0, lmax=0, sn=0.5, j=0.5, nam=1, ampl=1.0000,
    kbpot=4, be=-9.5173, isc=12, ipc=2, nk=20, er=-0.9655
)

f.set_coupling(icto=1, icfrom=2, kind=3, ip1=2, ip2=0, ip3=0, p1=3.00, p2=2.00,
    cfp=[
        {'in':1, 'ib':1, 'ia':1, 'kn':1, 'a':1},
        {'in':1, 'ib':2, 'ia':1, 'kn':2, 'a':1},
        {'in':1, 'ib':3, 'ia':1, 'kn':3, 'a':1},
        {'in':1, 'ib':4, 'ia':1, 'kn':4, 'a':1},
        {'in':1, 'ib':5, 'ia':1, 'kn':5, 'a':1},
        {'in':1, 'ib':6, 'ia':1, 'kn':6, 'a':1},
        {'in':1, 'ib':7, 'ia':1, 'kn':7, 'a':1},
        {'in':1, 'ib':8, 'ia':1, 'kn':8, 'a':1},
        {'in':1, 'ib':9, 'ia':1, 'kn':9, 'a':1},
        {'in':1, 'ib':10, 'ia':1, 'kn':10, 'a':1},
        {'in':1, 'ib':11, 'ia':1, 'kn':11, 'a':1},
        {'in':1, 'ib':12, 'ia':1, 'kn':12, 'a':1},
        {'in':1, 'ib':13, 'ia':1, 'kn':13, 'a':1},
        {'in':1, 'ib':14, 'ia':1, 'kn':14, 'a':1},
        {'in':1, 'ib':15, 'ia':1, 'kn':15, 'a':1},
        {'in':1, 'ib':16, 'ia':1, 'kn':16, 'a':1},
        {'in':1, 'ib':17, 'ia':1, 'kn':17, 'a':1},
        {'in':1, 'ib':18, 'ia':1, 'kn':18, 'a':1},
        {'in':1, 'ib':19, 'ia':1, 'kn':19, 'a':1},
        {'in':1, 'ib':20, 'ia':1, 'kn':20, 'a':1},
        {'in':1, 'ib':21, 'ia':1, 'kn':21, 'a':1},
    ]
)