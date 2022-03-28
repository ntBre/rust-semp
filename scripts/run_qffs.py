#!/usr/bin/env python3

import glob
import os
import shutil


def build_dirs(params_log):
    """params_log is a string filename containing the list of parameters for
    each iteration"""
    dirs = []
    with open(params_log, "r") as inf:
        buf = ""
        for line in inf:
            if "Iter" in line:
                it = "%03d" % int(line.split()[1])
                dirs.append(it)
            elif len(line.split()) == 0:
                try:
                    os.mkdir(it)
                except FileExistsError:
                    None
                with open(it + "/params.dat", "w") as p:
                    p.write(buf)
                # requires all the input files for pbqff to be in the current
                # directory
                for f in glob.glob("*.in"):
                    shutil.copy(f, it)
                buf = ""
            else:
                buf += line
    return dirs


def process_output(out):
    "out is a pbqff output filename"
    funds = []
    with open(out, "r") as inf:
        skip = 0
        in_tab = False
        for line in inf:
            if skip > 0:
                skip -= 1
            elif "ZPT" in line:
                skip = 3
                in_tab = True
            elif in_tab and "+" in line:
                break
            elif in_tab:
                sp = line.split()
                funds.append(float(sp[9]))
    return funds


if __name__ == "__main__":
    dirs = build_dirs("params.log")
    pbqff = "/home/brent/programs/pbqff/pbqff"
    for d in dirs:
        os.system("cd " + d + "; " + pbqff + " -o test.in")
        funds = process_output(d + "/test.out")
        print("%5s" % d, end="")
        for f in funds:
            print("%8.1f" % f, end="")
        print()
