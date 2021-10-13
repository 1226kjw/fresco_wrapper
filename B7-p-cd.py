from wfresco import *

f = Wfresco('B7-p-cd.frin', 'p + 112Cd elastic')

f.set_parameters(
	hcm=0.010,
	rmatch=20,
	jtmin=0,
	jtmax=200,
	theta_range=[0, 180, 2],
	xstabl=1,
	elab=27.9
)

f.set_partition(nex=1, qval=0,
	proj='Projton', massp=1, zp=1, target='112cd',
	states=[
        {
            'proj':["0.5+", 0],
            'target':["0.0+", 0],
            'cpot':1
        }
	]
)

f.set_pot(1,
	pots=[
		{'type':0, 'p':[112, 0, 1.2]},
		{'type':1, 'p':[52.5, 1.17, 0.75, 3.5, 1.32, 0.61]},
		{'type':2, 'p':[0, 0, 0, 8.5, 1.32, 0.61]},
		{'type':3, 'p':[6.2, 1.01, 0.75]},
	]
)