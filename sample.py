from wfresco import *

fresco = Input('sample.frin', '11Be + 197Au')
fresco.set_parameters(
	# hcm=0.1, rmatch=60, jtmin=0, jtmax=50, theta_range=[0,180,1],
	# chans=1, smats=2, xstabl=1, elab=[6.9,11,49.350]
	hcm         = 0.040,
	rmatch      = 50.000,
	rintp       = 0.24,
	rasym       = 340.00,
	accrcy      = 0.0010,
	jtmin       = 0.0,
	jtmax       = 1500.0,
	absend      = -0.0010,
	jump        = [10, 20, 0, 0, 0, 0],
	jbord       = [300.0, 1300.0, 0, 0, 0, 0],
	theta_range = [1.0, 180.0, 1.00],
	ips         = 0.0500,
	iblock      = 2,
	chans       = 1,
	smats       = 2,
	xstabl      = 1,
	nlpl        = 0,
	elab        = 42,
)

fresco.set_projectile('11BE')
fresco.set_target('197Au')
fresco.set_partition(qval=0.000, pwf='T', nex=2)

fresco.set_state(proj=["0.5+", 0], target=["0.0+",0.0], cpot=1)
fresco.set_state(proj=["0.5-", 0.3200], target=1, cpot=1)

# fresco.set_pot(1,
# 	{'ap':1.000, 'at':78.000, 'rc':1.2},
# 	{'type':1, 'p':[40.00, 1.2, 0.65, 10.0, 1.2, 0.500]}
# )

fresco.set_pot(1,
	{'type':COULOMB_POTENTIAL, 'p':[197.000, 11.000, 0.0010]},
	{'type':PROJECTILE_COUPLED, 'p':[197.000, 11.000, 0.0010], 'step':[{'ib':1,'ia':2,'k':1,'str':0.4817},{'ib':-2,'ia':1,'k':1,'str':0.4817}]},
	{'type':PROJECTILE_COUPLED, 'p':[0.497, 0.0000, 0.0000], 'itt':'F', 'shape':12,'step':[{'ib':1,'ia':2,'k':1,'str':0.4968},{'ib':-2,'ia':1,'k':1,'str':0.4968}]},
	{'type':CENTRAL_POTENTIAL_VOLUME, 'p':[40.000, 1.2290, 0.6120, 15.0000, 1.2290, 0.6120, 0.0000]},
)

fresco.write_input()
fresco.run()
#fresco.ls()
#fresco.ls(16)
fresco.show(16)
