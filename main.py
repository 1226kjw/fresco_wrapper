import pandas as pd
from fractions import Fraction

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

def get_particle(input_name):
    input_num = ''.join([i for i in input_name if i.isdigit()])
    input_alp = ''.join([i for i in input_name if i.isalpha()]).lower().title()
    parsed_str = input_num + input_alp
    if input_name == 'n' or (input_alp == 'N' and input_num == '1'):
        parsed_str = '1 n'
    idx = df[df['Name'] == parsed_str]
    return idx
def get_particle_set(comment: str):
    particle_set = []
    while True:
        input_name = input(comment + str(len(particle_set) + 1) + ': ')
        if input_name == '':
            break
        particle = get_particle(input_name)
        if particle.empty:
            print('There is no such particle in database:', input_name)
            continue
        particle_set.append(particle.index[0])
        print(particle)
        print(comment, list(map(lambda x: df.iloc[x]['Name'], particle_set)))
    return particle_set


# a,z,i,_0,Name,Mass excess,Mass excess err,Ex energy,Ex energy err,flag,Half-life,Unit,_1,JPI,_2,_3,_4,_5
name_list = ['a', 'z', 'i', '_0', 'Name', 'Mass excess', 'Mass excess err',
             'Ex energy', 'Ex energy err', 'flag', 'Half-life', 'Unit',
             '_1', 'JPI', '_2', '_3', '_4', '_5']
type_list = {i: str for i in name_list}
df = pd.read_csv('nubase2016.csv', names=name_list, dtype=type_list, keep_default_na=False)
pd.set_option('display.max_columns', None)
df['Mass'] = df['a'].apply(int) + df['Mass excess'].apply(lambda x: float(x.strip('#')) / 931494.013 if x != '' else 0)
df['PI'] = df['JPI'].apply(lambda x: -1 if '-' in x else 1)
df['J'] = df['JPI'].apply(parsing_spin)
df.drop(['_0', '_1', '_2', '_3', '_4', '_5', 'JPI', 'Mass excess', 'Mass excess err', 'flag', 'Half-life', 'Unit'],
        axis=1, inplace=True)
# a, z, i, Name, Ex energy, Ex energy err, flag, Half-life, Unit, Mass, PI, J: list

f = open('input.in', 'w')

param_list = ['hcm', 'rmatch', 'rintp', 'rasym', 'accrcy', 'jtmin', 'jtmax', 'absend', 'jump', 'jbord',
              'thmin', 'thmax', 'thinc', 'ips', 'iblock', 'chans', 'smats', 'xstabl', 'nlpl', 'elab']
params = {}
for i in param_list:
    param = input(i)
    params[i] = param
for i in params:
    f.write(i + '=' + params[i] + '\n')

projectile = get_particle_set('Projectile ')
target = get_particle_set('Target ')

