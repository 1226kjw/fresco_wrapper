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



indent = ' '
depth = 1

input_file_name = input('input file name: ')
f = open(input_file_name, 'w')
f.write(input('Title: ') + '\n')

f.write(indent * depth + 'NAMELIST\n')
depth += 1
param_list = ['hcm', 'rmatch', 'rintp', 'rasym', 'accrcy', 'jtmin', 'jtmax', 'absend', 'jump(1:6)', 'jbord',
              'thmin', 'thmax', 'thinc', 'ips', 'iblock', 'chans', 'smats', 'xstabl', 'nlpl', 'elab(1)']
params = {}
f.write(indent * depth + '&FRESCO\n')
depth += 1
for i in param_list:
    param = input(i + ': ')
    if param == '':
        continue
    params[i] = param
for i in params:
    f.write(indent * depth + i + '=' + params[i] + '\n')
f.write(indent * depth + '/\n')
depth -= 1


projectile = get_particle_set('Projectile ')
target = get_particle_set('Target ')

f.write(indent * depth + '&PARTITION\n')
depth += 1
f.write(indent * depth + 'namep=' + "'" + df.iloc[projectile[0]]['Name'][:8] + "'  ")
f.write(indent * depth + 'massp=' + str(df.iloc[projectile[0]]['Mass']) + "   ")


