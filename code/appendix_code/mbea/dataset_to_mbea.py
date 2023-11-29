
#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys, getopt, os

def AdjToMbeaFmt(infile, outfile):
    reader = open(infile, "r")
    writer = open(outfile, "w")
    left = 0 
    right = 0
    edges = 0
    for line in reader:
        nodelist = line.strip().split()
        if len(nodelist) != 0:
            for node in nodelist:
                right = max(right, int(node)+1)
                edges = edges + 1
            left = left + 1
    writer.write(str(left) + " " + str(right) + " " + str(edges) + "\n")
    i = 0
    reader.seek(0, 0)
    for line in reader:
        nodelist = line.strip().split()
        if len(nodelist) != 0:
            for node in nodelist:
                writer.write("u" + str(i) + "\t" + "v" + node + "\n")
            i = i + 1
    reader.close()
    writer.close()

def BatchProcessing(file_dir):
    dir = "db/"
    if not os.path.isdir(dir):
        os.mkdir(dir)
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if os.path.splitext(file)[1] == '.adj':
                AdjToMbeaFmt(root + "/"+file, dir + os.path.splitext(file)[0] + '.mbea')

def main(argv):
    inputfile = ''
    outputfile = ''
    file_dir = '.'
    helpinfo = 'Usage: python dataset_convert.py [option] [arg]\n\
        -i --infile <file> \t: input file\n\
        -o --outfile <file> \t: output file\n\
        -h \t\t:help information\n\
        -b <path> \t: batch processing mode\n'
    try:
        opts, args = getopt.getopt(argv, "hb:i:o:f:",["infile=", "outfile=", "format="])
    except getopt.GetoptError:
        print(helpinfo)
        sys.exit()    
    for opt, arg in opts:
        if opt == '-h':
            print(helpinfo)
            sys.exit()
        elif opt == '-b':
            if len(arg) != 0:
                file_dir = arg
            BatchProcessing(file_dir)
            sys.exit()
        elif opt in ("-i", "--infile"):
            inputfile = arg
        elif opt in ("-o", "--outfile"):
            outputfile = arg
    AdjToMbeaFmt(inputfile, outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])
