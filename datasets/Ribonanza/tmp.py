import glob

files = sorted(glob.glob("DMS_1.0/*"))


for file in files:
    print(file)
    with open(file) as inp:
        lines = inp.readlines()
        print(' '.join([str(x.split()[1]) for x in lines]))
