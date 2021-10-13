from wfresco import *

f = Wfresco('B2-example-inel.frin', 'alpha+c12 -> alpha+c12* @ 100 MeV; nuc def')

f.set_parameters(
    hcm=0.05,
    rmatch=20,
    jtmin=0.0,
    jtmax=40,
    absend=0.01,
    theta_range=[0,180,1],
    iter=1,
    ips=0,
    iblock=0,
    chans=1,
    smats=2,
    xstabl=1,
    elab=100
)

f.set_partition(
    qval=0.000, nex=2,
    proj='alpha', massp=4, zp=2, target='12C',
    states=[
        {
            'proj':["0.0+", 0],
            'target':["0.0+", 0],
            'cpot':1
        },
        {
            'proj':1,
            'target':["2.0+", 4.43],
            'cpot':1
        }
    ]
)

f.set_pot(1, ap=4, at=12, rc=1.2,
    pots=[
        {
            'type':CENTRAL_POTENTIAL_VOLUME,
            'p':[40, 1.2, 0.65, 10.0, 1.2, 0.5]
        },
        {
            'type':DEFORMED_TARGET,
            'p':[1.3]
        },
        {
            'type':CENTRAL_POTENTIAL_DERIVATIVE,
            'p':[0.0, 1.2, 0.65, 6.0, 1.2, 0.5]
        },
        {
            'type':DEFORMED_TARGET,
            'p':[1.3]
        },
    ]
)