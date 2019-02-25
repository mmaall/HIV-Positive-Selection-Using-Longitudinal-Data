#!/usr/bin/python

import sys
import Bio
import argparse


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument( '-i', dest='filename',  help='file name', required=False)
    args = p.parse_args()

    if args.filename is None:
        p.print_help()
        
    if args.filename != None:

        with open(args.filename,'r') as f:
            for line in f:
                scan(line)