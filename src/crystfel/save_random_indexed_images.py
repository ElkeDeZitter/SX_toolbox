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

save_random_indexed_crystals
-------
Script to select a random number of indexed images to a new stream file.
It will look into indexed images, independent of the number of indeexed crystals.
If you want to save a number of random indexed crystals, then use the 
save_random_indexed_crystals.py script (stil needs to be written)

Usage and example
-------
To get the help message:
python random_images_streamfile.py -h

To save 50 random images to a streamfile called my_output_50indexed.stream:
python random_images_streamfile.py -i my_input.stream -o my_output -n 50

"""
import os
import sys
import re
import argparse
from stream import stream

def select_indexed_images(stream_file, output_prefix, number):

    S = stream.Stream(stream_file)
    print("----> %s <---- " %(stream_file))
    S.select_indexing_methods(S.indexing_methods)
    _ = S.save_random_indexed_images(output_prefix, number)
    print("------------------")
    
def get_filename(fle, suffix):
    if "/" in fle:
        name = re.search(r"\/(.+?)\.%s" %(suffix), fle).group(1).split("/")[-1]
    else:
        name = re.sub("\.%s"%(suffix),"",fle)
        
    return name


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--stream_file', type=str, default='input.stream',help='Input stream file.')
    parser.add_argument('-o', '--output_prefix', type=str, default = None, help='Name prefix for the output stream file. The number of selected images will be mentioned in the output stream file anyway. If not provided, the prefix of the input file will be taken.')
    parser.add_argument('-n', '--number', type=int, default=0, help='Number of random images to be selected')
    
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

    if not isinstance(args.number, int):
        print("number should be an integer")
        sys.exit(1)
    else:
        number = args.number
        
    output_prefix = args.output_prefix
    if output_prefix == None:
        output_prefix = get_filename(stream_file, 'stream')
    
    select_indexed_images(stream_file, output_prefix, number)
