import sys

strainFile = sys.argv[1]
genome_ave = float(sys.argv[2])
output = open(sys.argv[3], "w")
baseRange = "nothing"

with open(strainFile, "r") as file:
    for line in file:
        if "ave" in line:
            line_list = line.split(" ")
            ave = float(line_list[2])
            ratio = ave/genome_ave
            output.write(baseRange + ":" + str(ratio) + "\n")
        else:
            baseRange = line.rstrip("\n")
            baseRange.rstrip("\n")
output.close()
