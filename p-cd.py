from wfresco import *

f = Wfresco('sample2.frin', '1-ADCC d + 208Pb')
f.set_parameters(
    hcm         = 0.10,
    rmatch      = 20.000,
    rintp       = 1,
    jtmin       = 0.0,
    jtmax       = 200.0,
    absend      = 0.0100,
    theta_range = [0.0, 180.0, 2.00],
    iter=1,
    chans=1,
    listcc=0,
    treneg=0,
    cdetr=0,
    smats=2,
    xstabl=1,
    elab=27.9
)

#f.set_partition(qval=0.000, pwf='T', nex=2)

# f.set_projectile('11Be')
# f.set_target('197Au')

# f.set_state(proj=["0.5+", 0], target=["0.0+",0.0], cpot=1)
# f.set_state(proj=["0.5-", 0.3200], target=1, cpot=1)

f.set_partition(qval=0.000, pwf=True, nex=1,
    proj='Proton', massp=1, zp=1,target='112Cd',
    state=[
        {
            'proj':["0.5+", 0],
            'target':["0.0+", 0],
            'cpot':1
        }
    ]
)

f.set_pot(1,
    {'type':COULOMB_POTENTIAL, 'p':[112, 0, 1.2]},
    {'type':CENTRAL_POTENTIAL_VOLUME, 'p':[52.5, 1.17, 0.75, 3.5, 1.32, 0.61, 0]},
    {'type':CENTRAL_POTENTIAL_DERIVATIVE, 'p':[0, 0, 0, 8.5, 1.32, 0.61, 0]},
    {'type':SPIN_ORBIT_PROJECTILE, 'p':[6.2, 1.01, 0.75]}
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

#f.show(16)
