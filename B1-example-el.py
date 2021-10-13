from wfresco import *

f = Wfresco('B1-example-el.frin', 'p+Ni78')

f.set_parameters(
    hcm=0.1,
    rmatch=60,
    jtmin=0.0,
    jtmax=50,
    absend=0.0010,
    theta_range=[0,180,1],
    chans=1,
    smats=2,
    xstabl=1,
    elab=[6.9, 11.0, 49.350],
    nlab=[1,1]
)

f.set_partition(
    qval=-0.000, nex=1,
    proj='p', massp=1, zp=1, target='Ni78',
    states=[
        {
            'proj':["0.5+", 0],
            'target':["0.0+", 0],
            'cpot':1
        }
    ]
)

f.set_pot(1, ap=1, at=78, rc=1.2,
    pots=[
        {
            'type':CENTRAL_POTENTIAL_VOLUME,
            'p':[40, 1.2, 0.65, 10.0, 1.2, 0.5]
        }
    ]
)