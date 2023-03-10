#!/usr/bin/env python3
"""Executes run.sh file (sbatch run.sh) in all subdirectories
of given directory"""

import os
import argparse

def get_args():
    """Gets the command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dir_path",
        nargs="?",
        type=str,
        default=os.getcwd(),
        help="Path to directory (default: current directory)",
    )
    return parser.parse_args()


def get_dirs(path):
    """Constructs list with subdirectories for given directory."""
    dirs = list()
    for directory in os.listdir(path):
        if os.path.isdir(os.path.join(path, directory)):
            dirs.append(directory)
    return dirs


def run_calcs(path):
    """Makes loop over subdirectories and executes run.sh file"""
    dirs = get_dirs(path)

    for directory in sorted(dirs):
        os.chdir(os.path.join(path, directory))
        os.system("sbatch run.sh")


def main():
    """Gets directory path, runs calculation in each subdirectory"""
    args = get_args()
    run_calcs(args.dir_path)


if __name__ == "__main__":
    main()
