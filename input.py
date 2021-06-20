import os
import subprocess
import pandas as pd
from fractions import Fraction

COULOMB_POTENTIAL = 0
CENTRAL_POTENTIAL_VOLUME = 1
CENTRAL_POTENTIAL_DERIVATIVE = 2
SPIN_ORBIT_PROJECTILE = 3
SPIN_ORBIT_TARGET = 4
TR_TENSOR_FORCE_PROJECTILE = 5
TR_TENSOR_FORCE_TARGET = 6
TENSOR_FORCE_COMBINED = 7
SPIN_FORCE = 8
DEFORMED_PROJECTILE = 10
DEFORMED_TARGET = 11
PROJECTILE_COUPLED = 12
TARGET_COUPLED = 13
L_L1_CENTRAL_POTENTIAL = 30

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
        self.new_dir = self.filename[:self.filename.rfind('.') if '.' in self.filename else None]
        self.output = self.new_dir + '.out'
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
                    f.write('  %s(1:%d)=' % (i, len(self.parameters[i])))
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
        try:
            os.mkdir(self.new_dir)
        except FileExistsError:
            for i in os.listdir(self.new_dir):
                os.remove('%s/%s' % (self.new_dir, i))
        self.output = self.new_dir + '/' + self.output
        subprocess.call("fresco < %s > %s" % (self.filename, self.output), shell=True)
        subprocess.call("mv fort.* %s" % self.new_dir, shell=True)
        subprocess.call("mv %s %s" % (self.filename, self.new_dir), shell=True)
