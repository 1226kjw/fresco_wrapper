import pandas as pd
from fractions import Fraction as frac

def parsing_spin(j: str):
    j = ''.join(i for i in j if i not in '()#-+*T=<>.' and not i.islower())
    j = j[:j.find(' ')] if j.find(' ') != -1 else j
    spins = []
    if ':' in j:
        start = float(frac(j[:j.find(':')]))
        end = float(frac(j[j.find(':') + 1:]))
        while start <= end:
            spins.append(start)
            start += 1
    elif j != '':
        spins = list(map(lambda x: float(frac(x) if not x.isdigit() else x), j.split(',')))

    return spins


name_list = ['a', 'z', 'i', '_0', 'Name', 'Mass excess', 'Mass excess err',
             'Ex energy', 'Ex energy err', 'flag', 'Half-life', 'Unit',
             '_1', 'JPI', '_2', '_3', '_4', '_5']

type_list = {i: str for i in name_list}

# a,z,i,_0,Name,Mass excess,Mass excess err,Ex energy,Ex energy err,flag,Half-life,Unit,_1,JPI,_2,_3,_4,_5
df = pd.read_csv('nubase2016.csv', names=name_list, dtype=type_list, keep_default_na=False)
df['PI'] = df['JPI'].apply(lambda x: 1 if '+' in x else -1)
df['J'] = df['JPI'].apply(parsing_spin)
df.drop(['JPI'], axis=1, inplace=True)
# a,z,i,_0,Name,Mass excess,Mass excess err,Ex energy,Ex energy err,flag,Half-life,Unit,_1,_2,_3,_4,_5,PI,J

input_name = input()
input_num = ''.join([i for i in input_name if i.isdigit()])
input_alp = ''.join([i for i in input_name if i.isalpha()]).lower().title()
parsed_str = input_num + input_alp
if input_name == 'n' or (input_alp == 'N' and input_num == '1'):
    parsed_str = '1 n'
idx = df[df['Name'] == parsed_str]
pd.set_option('display.max_rows', None)


print(df.iloc[341])
