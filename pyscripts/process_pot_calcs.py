#!/usr/bin/env python3
"""Process given directory and its subdirectories after OEP
calculations with MOLPRO"""

import os
import glob
import argparse
from pathlib import Path

def get_dirs(path):
    """Constructs list with subdirectories for given directory."""
    dirs = list()
    for directory in os.listdir(path):
        if os.path.isdir(os.path.join(path, directory)):
            dirs.append(directory)
    return dirs

def process_dir(dir_path):
    """Processes a folder after calculation which store dat-file with
    potentials: removes unuseful files, renames dat-files from last iteration
    with potentials along z-axis.
    At the moment, works only for EXX-OEP calculations, but can be generazed
    to RPA-OEP."""

    for filename in Path(dir_path).glob("*node*"):
       os.unlink(filename)

    for filename in Path(dir_path).glob("slurm*"):
       os.unlink(filename)

    for filename in Path(dir_path).glob("error"):
       os.unlink(filename)

    for filename in Path(dir_path).glob("molpro.xml"):
       os.remove(filename)

    for filename in Path(dir_path).glob("run.sh"):
       os.remove(filename)

    for filename in Path(dir_path).glob("*.x"):
       os.unlink(filename)

    for filename in Path(dir_path).glob("*.y*"):
       os.unlink(filename)

    potfiles = glob.glob(dir_path + "/vx-*.z")
    for potfile in sorted(potfiles)[:-1]:
        os.remove(potfile)
    potfiles = glob.glob(dir_path + "/vx-*.z")
    if potfiles != []:
        os.rename(potfiles[0], dir_path + "/vx.z")

    potfiles = glob.glob(dir_path + "/vxa-*.z")
    for potfile in sorted(potfiles)[:-1]:
        os.remove(potfile)
    potfiles = glob.glob(dir_path + "/vxa-*.z")
    if potfiles != []:
        os.rename(potfiles[0], dir_path + "/vxa.z")

    potfiles = glob.glob(dir_path + "/vxb-*.z")
    for potfile in sorted(potfiles)[:-1]:
        os.remove(potfile)
    potfiles = glob.glob(dir_path + "/vxb-*.z")
    if potfiles != []:
        os.rename(potfiles[0], dir_path + "/vxb.z")

    potfiles = glob.glob(dir_path + "/vc-*.z")
    for potfile in sorted(potfiles)[:-1]:
        os.remove(potfile)
    potfiles = glob.glob(dir_path + "/vc-*.z")
    if potfiles != []: 
        os.rename(potfiles[0], dir_path + "/vc.z")

    potfiles = glob.glob(dir_path + "/vca-*.z")
    for potfile in sorted(potfiles)[:-1]:
        os.remove(potfile)
    potfiles = glob.glob(dir_path + "/vca-*.z")
    if potfiles != []: 
        os.rename(potfiles[0], dir_path + "/vca.z")

    potfiles = glob.glob(dir_path + "/vcb-*.z")
    for potfile in sorted(potfiles)[:-1]:
        os.remove(potfile)
    potfiles = glob.glob(dir_path + "/vcb-*.z")
    if potfiles != []: 
        os.rename(potfiles[0], dir_path + "/vcb.z")


def get_args():
    """Gets the command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dir_path",
        nargs="?",
        type=str,
        default=".",
        help="Path to directory (default: current directory)",
    )
    return parser.parse_args()


def main():
    """Gets directory path, process the directory and subdirectories"""
    args = get_args()
    process_dir(args.dir_path)
    subdirs = get_dirs(args.dir_path)
    for subdir in subdirs:
        process_dir(os.path.join(args.dir_path, subdir))


if __name__ == "__main__":
    main()
