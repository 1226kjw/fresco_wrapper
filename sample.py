import wfresco

fresco = wfresco.Input('sample.frin', '11Be + 197Au')
fresco.set_parameters(
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

fresco.set_projectile("11Be")
fresco.set_target("197Au")

fresco.set_state({
    'jp'        : 0.5,
    'ptyp'      : 1,
    'ep'        : 0.0000,
    'cpot'      : 1,
    'jt'        : 0.0,
    'ptyt'      : 1,
    'et'        : 0.0000,
})
fresco.set_state({
    'jp'        : 0.5,
    'ptyp'      : -1,
    'ep'        : 0.3200,
    'cpot'      : 1,
    'copyt'     : 1,
})

fresco.set_pot({
    'kp'        : 1,
    'type'      : wfresco.COULOMB_POTENTIAL,
    'p'         : [197.000, 11.0000, 0.0010],
})
fresco.set_pot({
    'kp'        : 1,
    'type'      : wfresco.PROJECTILE_COUPLED,
    'p'         : [0.482, 0.0000, 0.0000]},
    step=[
    {
        'ib'    : 1,
        'ia'    : 2,
        'k'     : 1,
        'str'   : 0.4817
    },
    {
        'ib'    : -2,
        'ia'    : 1,
        'k'     : 1,
        'str'   : 0.4817
    }]
)

fresco.set_pot({
    'kp'        : 1,
    'type'      : wfresco.CENTRAL_POTENTIAL_VOLUME,
    'p'         : [40.000, 1.2290, 0.6120, 15.0000, 1.2290, 0.6120, 0.0000],
})
fresco.set_pot({
    'kp'        : 1,
    'type'      : wfresco.PROJECTILE_COUPLED,
    'itt'       : 'F',
    'shape'     : 12,
    'p'         : [0.497, 0.0000, 0.0000]},
    step=[
    {
        'ib'    : 1,
        'ia'    : 2,
        'k'     : 1,
        'str'   : 0.4968
    },
    {
        'ib'    : -2,
        'ia'    : 1,
        'k'     : 1,
        'str'   : 0.4968
    }]
)

fresco.write_input()
fresco.run()

fresco.show(16,201)
