import sys
import os
import shlex
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import re
from fractions import Fraction

"""
FRESCO must be executable and added in $PATH to execute by command 'fresco'
"""

if __name__ == '__main__': # ex) "python3 wfresco.py inputs/be11.frin" will make inputs/be11.py
    if len(sys.argv) <= 1:
        print('input file required')
        exit(1)
    for arg in sys.argv[1:]:
        with open(arg, 'r') as f:
            dot = arg.rfind('.')
            fpy = open(arg[:dot if dot != -1 else None] + '.py', 'w')
            fpy.write('from wfresco import *\n\n')
            comment = f.readline().strip()
            if comment[-1] == '\n':
                comment = comment[:-1]
            while True:
                line = f.readline()
                if not line:
                    break
                if line.strip().upper() == 'NAMELIST':
                    fpy.write("f = Wfresco('%s', '%s')\n" % (arg[arg.rfind('/') + 1:], comment))
                else:
                    line = line.strip()
                    if line and '&' == line[0]:
                        line = line.replace('&', 'f.')
                        while '/' not in line:
                            line += ' ' + f.readline().strip()
                        line = re.sub(r"('[^']*')| +", lambda x: x.group() if '\'' in x.group() else ' ', line) #공백여러개를 한개로 줄이기, 따옴표 내부 공백 무시
                        line = re.sub(' ?= ?', '=', line)
                        if ' ' not in line:
                            continue
                        line = line.replace(' ', '(', 1)
                        line = line[:2] + line[2:line.find('(')].upper() + line[line.find('('):]
                        line = re.sub(' */', ')', line)					
                        line = line.replace('!', '#')
                        i = line.find('(') + 1
                        j = line.rfind(')', 0, None if '#' not in line else line.find('#'))
                        parsing = line[i:j]
                        parsing = re.sub(r'\([^)]*\)', '', parsing)
                        parsing_list = re.findall("(?:\'.*?\'|\S)+", parsing) #공백기준 스플릿, 따옴표내부 공백은 무시
                        for pi in range(len(parsing_list) - 1, 0, -1):      #jump(1:6)=0 0 0 0 0 0 -> jump=[0,0,0,0,0,0]
                            if not parsing_list[pi][0].isalpha():
                                parsing_list[pi - 1] += ' ' + parsing_list[pi]
                                del(parsing_list[pi])
                            elif '\'' not in parsing_list[pi] and ' ' in parsing_list[pi]:
                                parsing_list[pi] = parsing_list[pi].replace('=', '=[') + ']'
                        line = line[:i] + re.sub(r"('[^']*')|([a-zA-Z][a-zA-Z0-9]*=)", lambda x: x.group().upper() if '\'' not in x.group() else x.group(), ' '.join(parsing_list)) + line[j:]
                        #등호 이전 키워드 대문자로 변경
                        line = re.sub(r"('[^']*')| +", lambda x: x.group() if x.group()[0] in '"\'' else ', ', line[:line.find('#') if '#' in line else None]) + (line[line.find('#'):] if '#' in line else '')
                        line = line.replace(', #', ' #', 1)
                        line = re.sub(r"=([TF])", lambda x: '=\'' + x.group(1) + '\'', line)
                        if '()' not in line[:line.find('#') if '#' in line else None]:
                            fpy.write(line + '\n')
                    elif line and '!' == line[0]:
                        line = line.replace('!', '#')
                        fpy.write(line + '\n')
                    elif line:
                        fpy.write('#' + line + '\n')
                    else:
                        fpy.write(line + '\n')
            fpy.close()

else:
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

    class Wfresco:
        def __init__(self, filename: str, comment: str):
            pd.set_option('display.max_rows', None)
            pd.set_option('display.width', None)
            if '/' in filename:
                print("'/' in filename!")
                return
            if '\n' in comment:
                print("newline in comment!")
                return
            self.filename = filename
            self.new_dir = self.filename[:self.filename.rfind('.') if '.' in self.filename else None]
            try:
                os.mkdir(self.new_dir)
            except FileExistsError:
                print(self.new_dir+'/', 'is already exists. Overwrite? (y/n)')
                choice = input()
                if choice != 'y' and choice != 'Y':
                    return
                for i in os.listdir(self.new_dir):
                    if i[i.rfind('.') + 1:].isdigit():
                        os.remove('%s/%s' % (self.new_dir, i))
            self.output = self.new_dir + '.out'
            self.comment = comment[:119] + '\n'
            self.partition = []
            self.parameters = {}
            self.projectile = []
            self.target = []
            self.potentials = []
            self.overlap = []
            self.coupling = []
            self.nubase = None
            self.fileinfo = None
            self.gui = True
            self.arginfo = pd.read_excel('nml_fresco.xlsx', sheet_name=None)
            for i in self.arginfo:
                self.arginfo[i].fillna('', inplace=True)
            self.iswritten = False
            self.isexecuted = False

        def parse_nubase(self):
            if self.nubase is not None:
                return
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
            name_list   = [
                'a', 'z', 'i', '_0', 'Name', 'Mass excess', 'Mass excess err',
                'Ex energy', 'Ex energy err', 'flag', 'Half-life', 'Unit',
                '_1', 'JPI', '_2', '_3', '_4', '_5'
            ]
            type_list   = {i: str for i in name_list}
            self.nubase = pd.read_csv('nubase2016.csv', names=name_list, dtype=type_list, keep_default_na=False)
            pd.set_option('display.max_columns', None)
            self.nubase['a']     = self.nubase['a'].apply(int)
            self.nubase['z']     = self.nubase['z'].apply(int)
            self.nubase['Mass']  = self.nubase['a'] + self.nubase['Mass excess'].apply(lambda x: float(x.strip('#')) / 931494.013 if x != '' else 0)
            self.nubase['PI']    = self.nubase['JPI'].apply(lambda x: -1 if '-' in x else 1)
            self.nubase['J']     = self.nubase['JPI'].apply(parsing_spin)
            self.nubase.drop(['_0', '_1', '_2', '_3', '_4', '_5', 'JPI', 'Mass excess', 'Mass excess err', 'flag', 'Half-life', 'Unit'],
                    axis=1, inplace=True)

        def set_parameters(self, **kargs):
            self.iswritten = False
            self.isexecuted = False
            typedef = {'i':[int], 'f':[float, int], 'b':[bool, str], 's':[str], 'l':[list, int, float]}

            if 'theta_range' in kargs:
                theta_range = kargs['theta_range']
                if len(theta_range) == 2:
                    kargs['thmin'] = min(theta_range)
                    kargs['thmax'] = max(theta_range)
                    kargs['thinc'] = 1
                elif len(theta_range) == 3 and theta_range[2] > 0:
                    kargs['thmin'] = theta_range[0]
                    kargs['thmax'] = theta_range[1]
                    kargs['thinc'] = theta_range[2]
                else:
                    print('theta_range invalid')
                    return
                del(kargs['theta_range'])
            if 'jt_range' in kargs:
                kargs['jtmin'] = min(kargs['jt_range'])
                kargs['jtmax'] = max(kargs['jt_range'])
                del(kargs['jt_range'])

            for i in kargs:
                search = self.arginfo['FRESCO'][self.arginfo['FRESCO']['short_name'] == i.lower()]
                if search.empty:
                    print('Unknown arg:', i)
                    return
                if type(kargs[i]) not in typedef[search['type'].iloc[0]]:
                    print('Invalid type:', i)
                    print('Expected:', typedef[search['type'].iloc[0]])
                    print('Given:', type(kargs[i]))
                    return
                if search['type'].iloc[0] == 'b':
                    if kargs[i] is True:
                        kargs[i] = 'T'
                    elif kargs[i] is False:
                        kargs[i] = 'F'
                self.parameters[i] = kargs[i]

        def set_partition(self, qval=None, nex=None, pwf=None, proj=None, massp=None, zp=None, target=None, masst=None, zt=None, states=None):
            self.iswritten = False
            self.isexecuted = False
            partition = {}
            if proj is None or target is None:
                print('proj and target are necessary!')
                return
            if qval is not None:
                partition['qval'] = qval
            if nex is not None:
                partition['nex'] = nex
            if pwf is not None:
                if pwf is True:
                    pwf = 'T'
                elif pwf is False:
                    pwf = 'F'
                partition['pwf'] = pwf
            if massp is None or zp is None:
                self.parse_nubase()
                input_num = ''.join([i for i in proj if i.isdigit()])
                input_alp = ''.join([i for i in proj if i.isalpha()]).lower().title()
                parsed_str = input_num + input_alp
                if proj.lower() == 'n' or (input_alp == 'N' and input_num == '1'):
                    parsed_str = '1 n'
                search = self.nubase[self.nubase['Name'] == parsed_str]
                if search.empty:
                    print('Cannot find:', proj)
                    return
                else:
                    tmp = self.nubase.iloc[search.index[0]].copy()
                    if massp is None:
                        massp = tmp['Mass']
                    if zp is None:
                        zp = tmp['z']
            partition['namep'] = proj
            partition['massp'] = massp
            partition['zp'] = zp

            if masst is None or zt is None:
                self.parse_nubase()
                input_num = ''.join([i for i in target if i.isdigit()])
                input_alp = ''.join([i for i in target if i.isalpha()]).lower().title()
                parsed_str = input_num + input_alp
                if target.lower() == 'n' or (input_alp == 'N' and input_num == '1'):
                    parsed_str = '1 n'
                search = self.nubase[self.nubase['Name'] == parsed_str]
                if search.empty:
                    print('Cannot find:', target)
                    return
                else:
                    tmp = self.nubase.iloc[search.index[0]].copy()
                    if masst is None:
                        masst = tmp['Mass']
                    if zt is None:
                        zt = tmp['z']
            partition['namet'] = target
            partition['masst'] = masst
            partition['zt'] = zt

            if states is not None:
                if type(states) is not list:
                    print('state must be in list')
                    return
                partition['states'] = []
                for i in states:
                    st = {}
                    if (type(i['proj']) == list and len(i['proj']) != 2) or (type(i['target']) == list and len(i['target']) != 2):
                        print('States format invalid')
                        return
                    if type(i['proj']) == list:
                        st['ptyp'] = 1 if '+' in i['proj'][0] else -1
                        st['jp'] = float(''.join([i for i in i['proj'][0] if i not in '-+']))
                        st['ep'] = i['proj'][1]
                    else:
                        st['copyp'] = int(i['proj'])
                    if type(i['target']) == list:
                        st['ptyt'] = 1 if '+' in i['target'][0] else -1
                        st['jt'] = float(''.join([i for i in i['target'][0] if i not in '-+']))
                        st['et'] = i['target'][1]
                    else:
                        st['copyt'] = int(i['target'])
                    st['cpot'] = i['cpot']
                    partition['states'].append(st)
            self.partition.append(partition)

        def set_pot(self, kp, ap=None, at=None, rc=None, pots=None):
            self.iswritten = False
            self.isexecuted = False
            if type(kp) != int:
                print('kp in set_pot is invalid')
                return
            tmp = {}
            if ap is not None:
                tmp['ap'] = ap
            if at is not None:
                tmp['at'] = at
            if rc is not None:
                tmp['rc'] = rc
            if tmp:
                tmp['kp'] = kp
                self.potentials.append(tmp)
            for i in pots:
                if 'type' in i and (not (0 <= i['type'] <= 13 or i['type'] == 30) or i['type'] == 9):
                    print('potential type invalid')
                    return
                i['kp'] = kp
                for j in i:
                    if i[j] is True:
                        i[j] = 'T'
                    elif i[j] is False:
                        i[j] == 'F'
                if 'step' in i and i['step'][-1]['ib'] > 0:
                    i['step'][-1]['ib'] *= -1
                self.potentials.append(i)

        def set_overlap(self, kn1=None, kn2=None, ic1=None, ic2=None, _in=None, kind=None, ch1=None, nn=None, l=None, lmax=None, sn=None, 
                        ia=None, j=None, ib=None, kbpot=None, krpot=None, be=None, isc=None, ipc=None, nfl=None, nam=None, 
                        ampl=None, keep=None, dm=None, nk=None, er=None, e=None):
            self.iswritten = False
            self.isexecuted = False
            over = {'kn1':kn1, 'kn2':kn2, 'ic1':ic1, 'ic2':ic2, 'in':_in, 'kind':kind, 'ch1':ch1, 'nn':nn, 'l':l, 'lmax':lmax, 'sn':sn,
                    'ia':ia, 'j':j, 'ib':ib, 'kbpot':kbpot, 'krpot':krpot, 'be':be, 'isc':isc, 'ipc':ipc, 'nfl':nfl, 'nam':nam,
                    'ampl':ampl, 'keep':keep, 'dm':dm, 'nk':nk, 'er':er, 'e':e}
            self.overlap.append(over)

        def set_coupling(self, icto=None, icfrom=None, kind=None, ip1=None, ip2=None, ip3=None, p1=None, p2=None, jmax=None, rmax=None,
                        kfrag=None, kcore=None, cfp=None):
            self.iswritten = False
            self.isexecuted = False
            coup = {'icto':icto, 'icfrom':icfrom, 'kind':kind, 'ip1':ip1, 'ip2':ip2, 'ip3':ip3, 'p1':p1, 'p2':p2, 'jmax':jmax, 'rmax':rmax,
                    'kfrag':kfrag, 'kcore':kcore, 'cfp':cfp}
            self.coupling.append(coup)

        def __write_parameters(self, f):
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

        def __write_partition(self, f):
            def write_pt(partition):
                for i in partition:
                    if i == 'namep':
                        f.write("  namep='%8s'" % partition[i][:8])
                    elif i == 'massp':
                        f.write("  massp=%8.4f" % partition[i])
                    elif i == 'zp':
                        f.write("  zp=%3d\n" % partition[i])
                    elif i == 'namet':
                        f.write("  namet='%8s'" % partition[i][:8])
                    elif i == 'masst':
                        f.write("  masst=%8.4f" % partition[i])
                    elif i == 'zt':
                        f.write("  zt=%3d\n" % partition[i])
                    
            for i in self.partition:
                f.write(' &PARTITION\n')
                part = '  '
                write_pt(i)
                for j in i:
                    if j in ['states','namep','massp','zp','namet','masst','zt']:
                        continue
                    part += j + '=' + str(i[j])
                    part += ' '
                part += '/\n'
                f.write(part)
                if 'states' in i:
                    for j in i['states']:
                        f.write('   &STATES ')
                        for k in j:
                            f.write(' %s=%4g' % (k, j[k]))
                        f.write(' /\n')
            f.write(' &partition /\n\n')

        def __write_potential(self, f):
            for i in self.potentials:
                f.write(' &POT')
                for j in i:
                    if i[j] is None or j == 'step':
                        continue
                    if type(i[j]) is not list:
                        f.write(' %s=%3s' % (j, i[j]))
                    else:
                        f.write(' %s(1:%d)=' % (j, len(i[j])))
                        for k in i[j]:
                            f.write(' %s' % k)
                f.write(' /\n')
                if 'step' in i and i['step'] is not None:
                    for j in i['step']:
                        f.write('  &step')
                        for k in j:
                            f.write(' %s=%3g' % (k, j[k]))
                        f.write(' /\n')
            f.write(' &pot /\n\n')

        def __write_overlap(self, f):
            for i in self.overlap:
                f.write(' &OVERLAP')
                for j in i:
                    if i[j] is None:
                        continue
                    if type(i[j]) is not list:
                        f.write(' %s=%3s' % (j, i[j]))
                    else:
                        f.write(' %s(1:%d)=' % (j, len(i[j])))
                        for k in i[j]:
                            f.write(' %s' % k)
                f.write(' /\n')
            f.write(' &overlap /\n\n')

        def __write_coupling(self, f):
            for i in self.coupling:
                f.write(' &COUPLING')
                for j in i:
                    if i[j] is None or j == 'cfp':
                        continue
                    if type(i[j]) is not list:
                        f.write(' %s=%3s' % (j, i[j]))
                    else:
                        f.write(' %s(1:%d)=' % (j, len(i[j])))
                        for k in i[j]:
                            f.write(' %s' % k)
                f.write(' /\n')
                if 'cfp' in i and i['cfp'] is not None:
                    for j in i['cfp']:
                        f.write('  &cfp')
                        for k in j:
                            f.write(' %s=%3g' % (k, j[k]))
                        f.write(' /\n')
            f.write(' &coupling /\n\n')

        def write_input(self):
            if self.iswritten:
                return
            self.isexecuted = False
            print("Write input(" + self.filename + ")")
            backup = os.getcwd()
            os.chdir(self.new_dir)
            with open(self.filename, 'w') as f:
                f.write(self.comment)
                f.write('NAMELIST\n')
                self.__write_parameters(f)
                self.__write_partition(f)
                self.__write_potential(f)
                self.__write_overlap(f)
                self.__write_coupling(f)
            self.iswritten = True
            print("Done!")
            os.chdir(backup)

        def run(self):
            if not self.iswritten:
                self.write_input()
                print()
            if self.isexecuted:
                return
            print("Run fresco(fresco < %s > %s)" % (self.filename, self.output))
            backup = os.getcwd()
            os.chdir(self.new_dir)
            ret = subprocess.call("fresco < %s > %s" % (self.filename, self.output), shell=True)
            self.isexecuted = True
            print("Done!")
            os.chdir(backup)
            return ret

        def ls(self, *filenum):
            if not self.isexecuted:
                self.run()
            if self.fileinfo is None:
                self.fileinfo = {}
                with open('file-allo2.txt', 'r') as f:
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        line = [i.strip() for i in line.split('@')]
                        if line[0].isdigit():
                            if line[1] == 'F':
                                line[1] = 'Fix'
                            if line[1] == 'V':
                                line[1] = 'Var'
                            if line[3] == 'S':
                                line[3] = 'Seq'
                            if line[3] == 'R':
                                line[3] = 'Ran'
                            self.fileinfo[int(line[0])] = line[1:]
                        else:
                            self.fileinfo[line[0]] = line[1:]
            if len(filenum) == 0:
                outfile = sorted(os.listdir(self.new_dir), key=lambda x: int(x[x.find('.')+1:]) if 'fort.' in x else 0)
                for i in outfile:
                    if 'fort.' not in i:
                        continue
                    else:
                        try:
                            n = int(i[i.find('.')+1:])
                            print("fort.%-3d  %s" % (n, self.fileinfo[n][4]))
                        except:
                            continue
                print()
            else:
                for i in filenum:
                    info = self.fileinfo[i]
                    column = self.fileinfo['File']
                    print('%10s : %10s' %('File', 'fort.'+str(i)))
                    for j in range(5):
                        print('%10s : %10s' % (column[j], info[j]))
                    print()

        def plot(self, *fileno: tuple):
            if not self.isexecuted:
                self.run()
            graph_len = len(fileno)
            g_row = 1
            g_col = 1
            while g_row * g_col < graph_len:
                if g_col == g_row:
                    g_row += 1
                else:
                    g_col += 1
            plt.rcParams["figure.figsize"] = (6 * g_col, 5 * g_row)
            for num, i in enumerate(fileno):
                if type(i) == str and 'fort.' in i:
                    i = int(i[i.find('.') + 1:])
                elif type(i) == int:
                    graphinput = self.new_dir + '/fort.' + str(i)
                else:
                    print('Invalid characters in wfresco.plot()')
                    return
                print('Graph', num+1, ':', graphinput)
                plt.subplot(g_row, g_col, num + 1)
                with open(graphinput, 'r') as f:
                    d_x = []
                    d_y = []
                    l_string = ''
                    l_istrue = False
                    while True:
                        line = shlex.split(f.readline())
                        if not line:
                            break
                        if line[0][0] == '#':
                            continue
                        if line[0][0] == '@':
                            cmd = line[0][1:].lower()
                            if cmd == 'subtitle':
                                plt.title(line[1])
                            elif cmd == 'legend':
                                if line[1].lower() == 'on':
                                    l_istrue = True
                                elif line[1].lower() == 'string':
                                    l_string = line[3]
                            elif cmd == 'g0':
                                if line[1].lower() == 'type':
                                    if line[2].lower() == 'logy':
                                        plt.yscale('log')
                            elif cmd == 'xaxis':
                                if line[1].lower() == 'label':
                                    plt.xlabel(line[2])
                            elif cmd == 'yaxis':
                                if line[1].lower() == 'label':
                                    plt.ylabel(line[2])
                        elif line[0].lower() == 'end' or line[0][0] == '&':
                            plt.plot(d_x, d_y, marker='.', markersize=3, label=l_string)
                            plt.legend(loc=0)
                            if l_istrue:
                                plt.legend()
                            d_x = []
                            d_y = []
                        elif not line[0].isalpha():
                            d_x.append(float(line[0]))
                            d_y.append(float(line[1]))
                    if d_x and d_y:
                        plt.legend(loc=0)
                        plt.plot(d_x, d_y, marker='.', markersize=3, label=l_string)
                        #if l_istrue:
                        #    plt.legend()
            plt.tight_layout()
            if not self.gui:
                pngname = self.new_dir + '/' + ':'.join(map(str, fileno))
                plt.savefig(pngname)
                subprocess.call("wslview %s" % (pngname + '.png'), shell=True)
            else:
                plt.show()

        def get_table(self, n, df=False):
            if not self.isexecuted:
                self.run()
            if type(n) == str and 'fort.' in n:
                n = int(n[n.find('.') + 1:])
            if type(n) == int:
                fileinput = self.new_dir + '/fort.' + str(n)
            else:
                print("get_table: Invalid argument")
                return
            with open(fileinput, 'r') as f:
                value = []
                table_list = []
                while True:
                    line = shlex.split(f.readline())
                    if not line:
                        break
                    for i in range(len(line)):
                        try:
                            line[i] = float(line[i])
                        except:
                            continue
                    if type(line[0]) == float or line[0][0].isdigit():
                        value.append(line)
                    else:
                        if len(value) != 0:
                            if df == False:
                                table_list.append(value)
                            else:
                                table_list.append(pd.DataFrame(value))
                            value = []
                if len(value) != 0:
                    if df == False:
                        table_list.append(value)
                    else:
                        table_list.append(pd.DataFrame(value))
            return table_list if len(table_list) != 1 else table_list[0]



        def FRESCO(self, **kargs):
            for i in kargs:
                self.parameters[i.lower()] = kargs[i]

        def PARTITION(self, **kargs):
            partition = {}
            for i in kargs:
                partition[i.lower()] = kargs[i]
            self.partition.append(partition)

        def STATES(self, **kargs):
            state = {}
            for i in kargs:
                state[i.lower()] = kargs[i]
            if 'states' not in self.partition[-1]:
                self.partition[-1]['states'] = [state]
            else:
                self.partition[-1]['states'].append(state)

        def POT(self, **kargs):
            pot = {}
            for i in kargs:
                pot[i.lower()] = kargs[i]
            self.potentials.append(pot)

        def STEP(self, **kargs):
            step = {}
            for i in kargs:
                step[i.lower()] = kargs[i]
            if 'step' not in self.potentials[-1]:
                self.potentials[-1]['step'] = [step]
            else:
                self.potentials[-1]['step'].append(step)

        def OVERLAP(self, **kargs):
            overlap = {}
            for i in kargs:
                overlap[i.lower()] = kargs[i]
            self.overlap.append(overlap)

        def COUPLING(self, **kargs):
            coupling = {}
            for i in kargs:
                coupling[i.lower()] = kargs[i]
            self.coupling.append(coupling)

        def CFP(self, **kargs):
            cfp = {}
            for i in kargs:
                cfp[i.lower()] = kargs[i]
            if 'cfp' not in self.coupling[-1]:
                self.coupling[-1]['cfp'] = [cfp]
            else:
                self.coupling[-1]['cfp'].append(cfp)

        def print(self, fileno, limit=0xffffffff):
            if not self.isexecuted:
                self.run()
            if 'fort.' + str(fileno) not in os.listdir(self.new_dir):
                print('fort.'+str(fileno)+' is not exists')
                return
            with open(self.new_dir + '/fort.' + str(fileno)) as wfile:
                while limit != 0:
                    l = wfile.readline()
                    if not l:
                        break
                    print(l, end='')
                    limit -= 1
                print()

        def help(self, func=None, arg=None):
            if func is None:
                NotImplemented
            elif arg is None:
                print(self.arginfo[func])
            else:
                search = self.arginfo[func][self.arginfo[func]['short_name'] == arg]
                if search.empty:
                    print('Cannot find \'', arg, '\' in \'', func, '\'', sep='')
                    return
                else:
                    print(self.arginfo[func].iloc[search.index[0]])