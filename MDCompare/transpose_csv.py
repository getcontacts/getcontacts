with open("./sample_inputs/sample_gd.csv", 'r') as ropen:
	lines = [line.strip() for line in ropen.readlines()]

with open("./sample_inputs/sample_gd_2.csv", "w+") as wopen:
	for original_col in xrange(len(lines[0].split(','))):
		new_line = []
		for line in lines:
			token = line.split(',')[original_col]
			new_line.append(token)
		wopen.write("%s\n" % ','.join(new_line))