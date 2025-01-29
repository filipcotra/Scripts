import sys

N2_ratios = open(sys.argv[1], "r")
strain_ratios = open(sys.argv[2], "r")
output = open(sys.argv[3], "w")
line_list = []

with N2_ratios as file:
        for line in file:
                line = line.rstrip("\n")
                line_split = line.split(":")
                line_list.append(line_split[0])
                line_list.append(line_split[1])
N2_ratios.close()

with strain_ratios as file:
        for line in file:
                line = line.rstrip("\n")
                line_split = line.split(":")
                for baseRange in line_list:
                        if line_split[0] == baseRange:
                                numerator = float(line_split[1])
                                denominator = float(line_list[line_list.index(line_split[0]) + 1])
                                normalized = numerator/denominator
                                output.write(str(normalized) + "\n")
strain_ratios.close()
output.close()
