#!/usr/bin/env python3
import os
import sys
import argparse
from stream import stream

"""
Script to show the indexing statistics of a Stream file. The Stream file has to be given as the -i argument. 

Example
------
show_stream_stats.py -i my_fancy_experiment.stream
"""

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = 'Show the indexing statistics of a Stream file')

    parser.add_argument('-i', '--stream_file', type=str, default='input.stream',help='Input stream file')

    args = parser.parse_args()

    #print help if no arguments provided
    if len(sys.argv) < 2:
           parser.print_help()
           sys.exit(1)
    
    if os.path.isfile(args.stream_file):
        stream_file = args.stream_file
    else:
        print("File not found: {:s}".format(args.stream_file))
        sys.exit(1)

    S = stream.Stream(stream_file)
    print(("----> %s <---- " %(stream_file)))
    S.get_stream_summary()
    print("------------------")
