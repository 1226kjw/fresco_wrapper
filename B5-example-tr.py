from wfresco import *

f = Wfresco('B5-example-tr.frin', 'n14(f17,ne18)c13 @ 170 MeV')
f.set_parameters(
    hcm         = 0.03,
    rmatch      = 40,
    rintp       = 0.20,
	hnl			= 0.1,
	rnl			= 5.0,
	centre		= 0.0,
    jtmin       = 0.0,
    jtmax       = 120.0,
    absend      = -1.0,
    theta_range = [0.0, 140.0, 0.10],
    iter		= 1,
	nnu			= 36,
	chans		= 1,
	xstabl		= 1,
    elab		= 27.9
)

f.set_partition(nex=1,
    proj='f17', target='n14',
    states=[
        {
            'proj':["2.5+", 0],
            'target':["1.0+", 0],
            'cpot':1
        }
    ]
)
f.set_partition(nex=1, qval=3.6286,
    proj='ne18', target='c13',
    states=[
        {
            'proj':["0+", 0],
            'target':["0.5+", 0],
            'cpot':2
        }
    ]
)
f.set_pot(1, ap=17, at=14, rc=1.3,
	pots=[
		{'type':CENTRAL_POTENTIAL_VOLUME, 'p':[37.2, 1.2, 0.6, 21.6, 1.2, 0.69]},
	]
)
f.set_pot(2, ap=18, at=13, rc=1.3,
	pots=[
		{'type':CENTRAL_POTENTIAL_VOLUME, 'p':[37.2, 1.2, 0.6, 21.6, 1.2, 0.69]},
	]
)
f.set_pot(3, at=17, rc=1.2,
	pots=[
		{'type':CENTRAL_POTENTIAL_VOLUME, 'p':[50, 1.2, 0.65]},
		{'type':SPIN_ORBIT_PROJECTILE, 'p':[6, 1.2, 0.65]},
	]
)
f.set_pot(4, at=13, rc=1.2,
	pots=[
		{'type':CENTRAL_POTENTIAL_VOLUME, 'p':[50, 1.2, 0.65]},
		{'type':SPIN_ORBIT_PROJECTILE, 'p':[6, 1.2, 0.65]},
	]
)
f.set_pot(5, ap=17, at=14, rc=1.3,
	pots=[
		{'type':CENTRAL_POTENTIAL_VOLUME, 'p':[37.2, 1.2, 0.6, 21.6, 1.2, 0.69]},
	]
)

f.set_overlap(kn1=1, ic1=1, ic2=2, _in=1, kind=0, nn=1, l=2, sn=0.5, j=2.5, kbpot=3, be=3.922, isc=1, ipc=0)
f.set_overlap(kn1=2, ic1=2, ic2=1, _in=2, kind=3, nn=1, l=1, sn=0.5, ia=1, ib=1, j=1.0, kbpot=4, be=7.5506, isc=1, ipc=0)

f.set_coupling(icto=-2, icfrom=1, kind=7, ip1=0, ip2=-1, ip3=5,
	cfp=[
		{'in':1, 'ib':1, 'ia':1, 'kn':1, 'a':1},
		{'in':2, 'ib':1, 'ia':1, 'kn':2, 'a':1}
	]
)