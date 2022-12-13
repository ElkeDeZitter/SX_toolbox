#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolas Coquelle - nicolas.coquelle@ibs.fr
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr


license information
-------
Copyright (c) 2022 Elke De Zitter, Nicolas Coquelle, Jacques-Philippe Collettier
https://github.com/ElkeDeZitter/SX_toolbox/blob/main/LICENSE

show_stream_stats
-------
Script to show the indexing statistics of a Stream file. The Stream file has to be given as the -i argument. 

Usage and example
-------

To get the help message:
python show_stream_stats.py -h

To get the stats of your my_fancy_experiment.stream file:
python show_stream_stats.py -i my_fancy_experiment.stream

"""
import os
import sys
import argparse
from stream import stream


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
