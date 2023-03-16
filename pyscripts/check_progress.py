#!/usr/bin/env python3
"""Determines progress for a bunch of Molplro calculations running in
given directory."""

import os
import argparse

def get_dirs(path):
    """Constructs list with subdirectories for given directory."""
    dirs = list()
    for directory in os.listdir(path):
        if os.path.isdir(os.path.join(path, directory)):
            dirs.append(directory)
    return dirs


def get_progress(path):
    """Determines and prints progress of calculations in directory"""
    dirs = get_dirs(path)

    n_started = 0
    n_finished = 0
    for directory in dirs:
        if os.path.isfile(os.path.join(path, directory, "output")):
            n_started += 1
            with open(os.path.join(path, directory, "output")) as file_obj:
                if "Molpro calculation terminated" in file_obj.read():
                    n_finished += 1

    print("Started:", n_started, "/", len(dirs))
    print("Finished:", n_finished, "/", len(dirs))


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
    """Gets directory path, determines progress and prints the information"""
    args = get_args()
    get_progress(args.dir_path)


if __name__ == "__main__":
    main()
