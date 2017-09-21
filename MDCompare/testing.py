

with open("posed.alignment", 'r') as ropen:
	with open("posted2.alignment", 'w+') as wopen:
		for line in ropen.readlines():
			wopen.write("%s\n" % ','.join(line.split()))