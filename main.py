import pandas as pd

namelist = ['a', 'z', 'i', '_0', 'Name', 'Mass excess1', 'Mass excess2',
            'Excitation energy1', 'Ex energy2', 'flag', 'Half-life', 'Unit',
            '_1', 'JPI', '_2', '_3', '_4', '_5']

df = pd.read_csv('nubase2016.csv', names=namelist, keep_default_na=False)
input_name = input()
input_num = ''.join([i for i in input_name if i.isdigit()])
input_alp = ''.join([i for i in input_name if i.isalpha()]).lower().title()
parsed_str = input_num + input_alp
if input_name == 'n' or (input_alp == 'N' and input_num == '1'):
    parsed_str = '1 n'
idx = df[df['Name'] == parsed_str]
print(idx['Name'])
