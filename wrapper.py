import pandas as pd
from fractions import Fraction
import subprocess
import os

def write_pt(f, data, pt: str):
    if pt == 'p':
        f.write("  namep='%8s'" % data['Name'][:8])
        f.write("  massp=%8.4f" % data['Mass'])
        f.write("  zp=%3d\n" % data['z'])
    else:
        f.write("  namet='%8s'" % data['Name'][:8])
        f.write("  masst=%8.4f" % data['Mass'])
        f.write("  zt=%3d\n" % data['z'])

class Input:
    def __init__(self, filename: str, comment: str):
        if '/' in filename:
            print("'/' in filename!")
            exit(0)
        if '\n' in comment:
            print("newline in comment!")
            exit(0)
        self.filename = filename
        self.output = self.filename[:self.filename.find('.') if '.' in self.filename else None] + '.out'
        self.comment = comment + '\n'
        self.parameters = {}
        self.projectile = []
        self.target = []
        self.states = []
        self.potentials = []

    def set_parameters(self, par: dict):
        for i in par:
            self.parameters[i] = par[i]

    def set_projectile(self, name: str):
        input_num = ''.join([i for i in name if i.isdigit()])
        input_alp = ''.join([i for i in name if i.isalpha()]).lower().title()
        parsed_str = input_num + input_alp
        if name == 'n' or (input_alp == 'N' and input_num == '1'):
            parsed_str = '1 n'
        if df[df['Name'] == parsed_str].empty:
            print('Cannot find: ' + name)
        else:
            self.projectile.append(df[df['Name'] == parsed_str].index[0])

    def set_target(self, name: str):
        input_num = ''.join([i for i in name if i.isdigit()])
        input_alp = ''.join([i for i in name if i.isalpha()]).lower().title()
        parsed_str = input_num + input_alp
        if name == 'n' or (input_alp == 'N' and input_num == '1'):
            parsed_str = '1 n'
        if df[df['Name'] == parsed_str].empty:
            print('Cannot find: ' + name)
        else:
            self.target.append(df[df['Name'] == parsed_str].index[0])

    def set_state(self, state: dict):
        self.states.append(state)

    def set_pot(self, pot: dict, step=None):
        self.potentials.append(pot)
        if step:
            self.potentials[-1]['step'] = step

    def write_input(self):
        with open(self.filename, 'w') as f:
            f.write(self.comment)
            f.write('NAMELIST\n')
            f.write(' &FRESCO\n')
            for i in self.parameters:
                if type(self.parameters[i]) is list:
                    f.write('  %s=' % i)
                    for j in self.parameters[i]:
                        f.write('%g ' % j)
                    f.write('\n')
                else:
                    f.write(''.join(['  ', i, '=', str(self.parameters[i]), '\n']))
            f.write('  /\n\n')

            f.write(' &PARTITION\n')
            for i in self.projectile:
                write_pt(f, df.iloc[i], 'p')
            f.write('  nex=2  pwf=T\n')
            for i in self.target:
                write_pt(f, df.iloc[i], 't')
            f.write('  qval=0.0000 /\n')

            for i in self.states:
                f.write('   &STATES ')
                for j in i:
                    f.write(' %s=%4g' % (j, i[j]))
                f.write(' /\n')
            f.write(' &partition /\n\n')

            for i in self.potentials:
                f.write(' &pot')
                for j in i:
                    if j != 'step':
                        if type(i[j]) is not list:
                            f.write(' %s=%3s' % (j, i[j]))
                        else:
                            f.write(' %s(1:%d)=' % (j, len(i[j])))
                            for k in i[j]:
                                f.write(' %s' % k)
                f.write(' /\n')
                if 'step' in i:
                    for j in i['step']:
                        f.write('  &step')
                        for k in j:
                            f.write(' %s=%3g' % (k, j[k]))
                        f.write(' /\n')
            f.write(' &pot /\n\n')

            f.write(' &overlap /\n\n')
            f.write(' &coupling /')

    def run(self):
        new_dir = self.filename[:self.filename.find('.') if '.' in self.filename else None]
        try:
            os.mkdir(new_dir)
        except FileExistsError:
            for i in os.listdir(new_dir):
                os.remove('%s/%s' % (new_dir, i))
        self.output = new_dir + '/' + self.output + '.out'
        subprocess.call("fresco < %s > %s" % (self.filename, self.output), shell=True)
        subprocess.call("mv fort.* %s" % new_dir, shell=True)
        subprocess.call("mv %s %s" % (self.filename, self.output))



def parsing_spin(j: str):
    j = ''.join(i for i in j if i not in '()#-+*T=<>.' and not i.islower())
    j = j[:j.find(' ')] if j.find(' ') != -1 else j
    spins = []
    if ':' in j:
        start = float(Fraction(j[:j.find(':')]))
        end = float(Fraction(j[j.find(':') + 1:]))
        while start <= end:
            spins.append(start)
            start += 1
    elif j != '':
        spins = list(map(lambda x: float(Fraction(x) if '/' in x else x), j.split(',')))
    return spins

def get_custom():
    custom = {'Name': input('Name: ')}
    while not df[df['Name'] == custom['Name']].empty:
        print('Duplicated Name!')
        custom['Name'] = input('Name: ')
    for i in ['a', 'z', 'i', 'Ex energy', 'Mass', 'PI', 'J']:
        custom[i] = input(i + ': ')
    custom['Mass'] = float(custom['Mass'])
    if '-' in custom['PI']:
        custom['PI'] = -1
    else:
        custom['PI'] = 1
    if '/' in custom['J']:
        custom['J'] = float(Fraction(custom['J']))
    else:
        custom['J'] = float(custom['J'])
    df.append(custom, ignore_index=True)
    return df[df['Name'] == custom['Name']]

def get_particle(input_name):
    input_num = ''.join([i for i in input_name if i.isdigit()])
    input_alp = ''.join([i for i in input_name if i.isalpha()]).lower().title()
    parsed_str = input_num + input_alp
    if input_name == 'n' or (input_alp == 'N' and input_num == '1'):
        parsed_str = '1 n'
    return df[df['Name'] == parsed_str]

def get_particle_set(comment: str):
    particle_set = []
    while True:
        input_name = input('(input -1 to set custom particle or nothing to finish)\n' +
                           comment + str(len(particle_set) + 1) + ': ')
        if input_name == '':
            break
        elif input_name == '-1':
            particle = get_custom()
        else:
            particle = get_particle(input_name)
        if particle.empty:
            print('There is no such particle in database:', input_name)
            continue
        particle_set.append(particle.index[0])
        print(particle)
        print(comment, list(map(lambda x: df.iloc[x]['Name'], particle_set)))
    return particle_set


# -------------------------------------------------  Data Clearing  -------------------------------------------------- #
# a,z,i,_0,Name,Mass excess,Mass excess err,Ex energy,Ex energy err,flag,Half-life,Unit,_1,JPI,_2,_3,_4,_5
name_list   = [
    'a', 'z', 'i', '_0', 'Name', 'Mass excess', 'Mass excess err',
    'Ex energy', 'Ex energy err', 'flag', 'Half-life', 'Unit',
    '_1', 'JPI', '_2', '_3', '_4', '_5'
]
type_list   = {i: str for i in name_list}
df = pd.read_csv('nubase2016.csv', names=name_list, dtype=type_list, keep_default_na=False)
pd.set_option('display.max_columns', None)
df['a']     = df['a'].apply(int)
df['z']     = df['z'].apply(int)
df['Mass']  = df['a'] + df['Mass excess'].apply(lambda x: float(x.strip('#')) / 931494.013 if x != '' else 0)
df['PI']    = df['JPI'].apply(lambda x: -1 if '-' in x else 1)
df['J']     = df['JPI'].apply(parsing_spin)
df.drop(['_0', '_1', '_2', '_3', '_4', '_5', 'JPI', 'Mass excess', 'Mass excess err', 'flag', 'Half-life', 'Unit'],
        axis=1, inplace=True)
# a, z, i, Name, Ex energy, Ex energy err, flag, Half-life, Unit, Mass, PI, J: list



frin = Input('sample.frin', '11Be + 197Au')
frin.set_parameters({
    'hcm'       : 0.040,
    'rmatch'    : 50.000,
    'rintp'     : 0.24,
    'rasym'     : 340.00,
    'accrcy'    : 0.0010,
    'jtmin'     : 0.0,
    'jtmax'     : 1500.0,
    'absend'    : -0.0010,
    'jump(1:6)' : [10, 20, 0, 0, 0, 0],
    'jbord'     : [300.0, 1300.0, 0, 0, 0, 0],
    'thmin'     : 1.00,
    'thmax'     : 180.00,
    'thinc'     : 1.00,
    'ips'       : 0.0500,
    'iblock'    : 2,
    'chans'     : 1,
    'smats'     : 2,
    'xstabl'    : 1,
    'nlpl'      : 0,
    'elab(1)'   : 42,
})

frin.set_projectile("11Be")
frin.set_target("197Au")

frin.set_state({
    'jp'        : 0.5,
    'ptyp'      : 1,
    'ep'        : 0.0000,
    'cpot'      : 1,
    'jt'        : 0.0,
    'ptyt'      : 1,
    'et'        : 0.0000,
})
frin.set_state({
    'jp'        : 0.5,
    'ptyp'      : -1,
    'ep'        : 0.3200,
    'cpot'      : 1,
    'copyt'     : 1,
})

frin.set_pot({
    'kp'        : 1,
    'type'      : 0,
    'p'         : [197.000, 11.0000, 0.0010],
})
frin.set_pot({
    'kp'        : 1,
    'type'      : 12,
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

frin.set_pot({
    'kp'        : 1,
    'type'      : 1,
    'p'         : [40.000, 1.2290, 0.6120, 15.0000, 1.2290, 0.6120, 0.0000],
})
frin.set_pot({
    'kp'        : 1,
    'type'      : 12,
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

frin.write_input()
frin.run()
