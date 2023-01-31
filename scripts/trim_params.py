import subprocess
import sys

# In order to run this, it's assumed that mopac is installed in
# /opt/mopac/mopac and there must be a test.mop file in the working
# directory

if len(sys.argv) != 2:
    print("usage: trim_params.py INITIAL_PARAM_FILE")
    exit(1)

params = []
with open(sys.argv[1], "r") as infile:
    params = infile.read().splitlines()

prev = 0
print("Params deleted:")
while prev != len(params):
    prev = len(params)
    with open("params.dat", "w") as out:
        out.write(("\n").join(params))

    subprocess.run(["/opt/mopac/mopac", "test.mop"], capture_output=True)

    with open("test.out", "r") as inf:
        for line in inf:
            if "unrecognized" in line:
                param = line.split("'")[1].split()
                param = f"{param[0]:5}{param[1]:>6}{param[2]:>19}"
                print(param)
                params.remove(param)
