from wfresco import *

f = Wfresco('B6-example-capture.frin', '14C(n,g)15C E1 only')

f.set_parameters(
	hcm = 0.1,
	rmatch=100,
	jtmin=0,
	jtmax=4.5,
	absend=-1,
	theta_range=[0,0],
	iter=1,
	elab=[0.005, 4.005],
	nlab=50
)

f.set_partition(nex=1,
	proj='neutron', massp=1.0087, zp=0, target='14c',
	states=[
        {
            'proj':["0.5+", 0],
            'target':["0.0+", 0],
            'cpot':1
        }
	]
)
f.set_partition(nex=1, qval=1.218,
	proj='Gamma', massp=0, zp=0, target='15c',
	states=[
        {
            'proj':["1.0+", 0],
            'target':["0.5+", 0],
            'cpot':3
        }
	]
)

f.set_pot(1,
	pots=[
		{'type':0, 'shape':0, 'p':[14, 0, 1.3]},
		{'type':1, 'shape':0, 'p':[57, 1.7, 0.7]},
		{'type':3, 'shape':0, 'p':[0, 1.7, 0.5]}
	]
)
f.set_pot(2,
	pots=[
		{'type':0, 'shape':0, 'p':[14, 0, 1.2]},
		{'type':1, 'shape':0, 'p':[55.77, 1.223, 0.5]},
		{'type':3, 'shape':0, 'p':[0.5, 1.223, 0.5]}
	]
)

f.set_overlap(kn1=1, ic1=1, ic2=2, _in=-2, kind=0, nn=2, l=0, sn=0.5, j=0.5, kbpot=2, be=1.218, isc=1)
f.set_coupling(icto=2, icfrom=1, kind=2, ip1=-1, ip2= 1,
	cfp=[
		{'in':2, 'ib':1, 'ia':1, 'kn':1, 'a':1.000}
	]
)