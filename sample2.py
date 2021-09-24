from wfresco import *

f = Wfresco('sample2.frin', '1-ADCC d + 208Pb')
f.set_parameters(
    hcm         = 0.050,
    rmatch      = 30.000,
    rintp       = 0.40,
    rasym       = 300.00,
    accrcy      = 0.0010,
    hnl         = 0.025,
    rnl         = 1.50,
    centre      = 0.00,
    jtmin       = 0.0,
    jtmax       = 40.0,
    absend      = -1.0010,
    jump        = [4, 0, 0, 0, 0, 0],
    jbord       = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    theta_range = [20.0, 0.0, 20.00],
    cutl        = -1.60,
    cutr        = -10.00,
    cutc        = 0.00,
    nnu         = 24,
    mtmin       = 8,
    nrbases     = 10,
    nrbmin      = 5,
    buttle      = 0,
    pralpha     = 'F',
    pcon        = 0,
    rmatr       = 20.00,
    chans       = 1,
    listcc      = 1,
    smats       = 4,
    xstabl      = 0,
    veff        = 1,
    elab        = 50,
)

f.set_partition(qval=0.000, nex=11,proj='Deuteron',massp=2,zp=1,target='208Pb',
    state=[
        {
            #proj=["0.0+", 0], target=["0.0+",0.0], cpot=1)
        }
    ]
)
f.set_projectile('Deuteron', mass=2,z=1)
f.set_target('208Pb', mass=208, z=82)

f.set_state(proj=["0.0+", 0], target=["0.0+",0.0], cpot=1)
f.set_state(proj=["1.0+", 0], target=1, cpot=1)
f.set_state(proj=["1.0+", 0], target=1, cpot=1)
f.set_state(proj=["1.0+", 0], target=1, cpot=1)
f.set_state(proj=["1.0+", 0], target=1, cpot=1)
f.set_state(proj=["1.0+", 0], target=1, cpot=1)



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
#f.ls(16)
#table = f.get_table(16)
# for i in table:
#     print(i)
#     for j in table[i]:
#         print(j)

f.show(16)
