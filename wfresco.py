import os
import math
import shlex
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
			exit(1)
		if '\n' in comment:
			print("newline in comment!")
			exit(1)
		self.filename = filename
		self.new_dir = self.filename[:self.filename.rfind('.') if '.' in self.filename else None]
		self.output = self.new_dir + '.out'
		self.comment = comment + '\n'
		self.parameters = {}
		self.projectile = []
		self.target = []
		self.states = []
		self.potentials = []

	def set_parameters(self, hcm=None, rmatch=None, rintp=None, rasym=None, accrcy=None, jtmin=None, jtmax=None,
					   absend=None, jump=None, jbord=None, theta_range=None, ips=None,
					   iblock=None, chans=None, smats=None, xstabl=None, nlpl=None, elab=None):
		if hcm:
			self.parameters['hcm'] = hcm
		if rmatch:
			self.parameters['rmatch'] = rmatch
		if rintp:
			self.parameters['rintp'] = rintp
		if rasym:
			self.parameters['rasym'] = rasym
		if accrcy:
			self.parameters['accrcy'] = accrcy
		if jtmin:
			self.parameters['jtmin'] = jtmin
		if jtmax:
			self.parameters['jtmax'] = jtmax
		if absend:
			self.parameters['absend'] = absend
		if jump:
			self.parameters['jump'] = jump
		if jbord:
			self.parameters['jbord'] = jbord
		if theta_range:
			if len(theta_range) == 2:
				self.parameters['thmin'] = min(theta_range)
				self.parameters['thmax'] = max(theta_range)
				self.parameters['thinc'] = 1
			elif len(theta_range) == 3 and theta_range[2] > 0:
				self.parameters['thmin'] = theta_range[0]
				self.parameters['thmax'] = theta_range[1]
				self.parameters['thinc'] = theta_range[2]
			else:
				print('theta_range invalid')
				exit(1)
		if ips:
			self.parameters['ips'] = ips
		if iblock:
			self.parameters['iblock'] = iblock
		if chans:
			self.parameters['chans'] = chans
		if smats:
			self.parameters['smats'] = smats
		if xstabl:
			self.parameters['xstabl'] = xstabl
		if nlpl:
			self.parameters['nlpl'] = nlpl
		if elab:
			self.parameters['elab'] = elab

	def set_projectile(self, name: str):
		input_num = ''.join([i for i in name if i.isdigit()])
		input_alp = ''.join([i for i in name if i.isalpha()]).lower().title()
		parsed_str = input_num + input_alp
		if name == 'n' or (input_alp == 'N' and input_num == '1'):
			parsed_str = '1 n'
		if df[df['Name'] == parsed_str].empty:
			print('Cannot find: ' + name)
			exit(2)
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
			exit(2)
		else:
			self.target.append(df[df['Name'] == parsed_str].index[0])

	def set_state(self, proj=None, target=None, cpot=1):
		if (type(proj) != int and len(proj) != 2) or (type(target) != int and len(target) != 2):
			print('States format invalid')
			exit(3)
		state = {}
		if type(proj) == list:
			state['ptyp'] = 1 if '+' in proj[0] else -1
			state['jp'] = float(''.join([i for i in proj[0] if i not in '-+']))
			state['ep'] = proj[1]
		else:
			state['copyp'] = proj
		if type(target) == list:
			state['ptyt'] = 1 if '+' in target[0] else -1
			state['jt'] = float(''.join([i for i in target[0] if i not in '-+']))
			state['et'] = target[1]
		else:
			state['copyt'] = target
		state['cpot'] = cpot
		self.states.append(state)

	def set_pot(self, kp=None, *pots):
		if not kp or type(kp) != int:
			print('kp in set_pot is invalid')
			exit(1)
		for i in pots:
			if not (0 <= i['type'] <= 13 or i['type'] == 30) or i['type'] == 9:
				print('potential type invalid')
				exit(1)
			i['kp'] = kp
			self.potentials.append(i)

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

	def show(self, *fileno: tuple):
		graph_len = len(fileno)
		graph_len_sqrt = math.ceil(math.sqrt(graph_len))
		g_row = 1
		g_col = 1
		while g_row * g_col < graph_len:
			if g_col == g_row:
				g_row += 1
			else:
				g_col += 1
		plt.rcParams["figure.figsize"] = (6 * g_col, 5 * g_row)
		for num, i in enumerate(fileno):
			graphinput = self.new_dir + '/fort.' + str(i)
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
						plt.plot(d_x, d_y, label=l_string)
						plt.legend(loc=0)
						if l_istrue:
							plt.legend()
						#l_string = ''
						d_x = []
						d_y = []
					else:
						d_x.append(float(line[0]))
						d_y.append(float(line[1]))
				if d_x and d_y:
					plt.plot(d_x, d_y, label=l_string)
					plt.legend(loc=0)
					#if l_istrue:
					#    plt.legend()
		plt.show()
