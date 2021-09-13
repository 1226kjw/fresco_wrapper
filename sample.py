from wfresco import *
#import wfresco

f = Wfresco('sample.frin', '11Be + 197Au')
f.set_parameters(
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

f.set_partition(qval=0.000, pwf='T', nex=2)
f.set_projectile('11BE')
f.set_target('197Au')

f.set_state(proj=["0.5+", 0], target=["0.0+",0.0], cpot=1)
f.set_state(proj=["0.5-", 0.3200], target=1, cpot=1)

f.set_pot(1,
    #fresco.coulomb_potential(p=[197.000, 11.000, 0.0010], shape=12, step=[{'ib':1,'ia':2,'k':1,'str':0.4817},{'ib':-2,'ia':1,'k':1,'str':0.4817}]),
    {'type':COULOMB_POTENTIAL, 'p':[197.000, 11.000, 0.0010]},
    {'type':PROJECTILE_COUPLED, 'p':[197.000, 11.000, 0.0010], 'step':[{'ib':1,'ia':2,'k':1,'str':0.4817},{'ib':-2,'ia':1,'k':1,'str':0.4817}] },
    {'type':CENTRAL_POTENTIAL_VOLUME, 'p':[40.000, 1.2290, 0.6120, 15.0000, 1.2290, 0.6120, 0.0000]},
    {'type':PROJECTILE_COUPLED, 'p':[0.497, 0.0000, 0.0000], 'itt':False, 'shape':12,'step':[{'ib':1,'ia':2,'k':1,'str':0.4968},{'ib':-2,'ia':1,'k':1,'str':0.4968}]},
)

f.write_input()
f.run()

#table = fresco.get_table(16)

f.ls()
f.ls(16)
#f.show(16)
table = f.get_table(16)
for i in table:
    print(i)
    for j in table[i]:
        print(j)

f.show(16)