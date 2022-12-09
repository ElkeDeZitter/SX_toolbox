from __future__ import print_function
import os, sys, re, numpy as np#, matplotlib.pyplot as plt
import random
"""
CLass that efficiently reads in a crystfel stream file and allow calculate statistics from the stream file as well as modifying it.
The initiat Stream.py was written by Nicolas Coquelle and later Modified by Elke De Zitter
This class should be opened with python2.7, it is not working with python3 (yet)!
"""

class Stream(object):
    def __init__(self, streamfile):
        self.streamfile = streamfile
        self.frames = []
        self.header = ''
        self.parse_stream()

    def parse_stream(self):
        append = 0
        crystal = 0
        count_shots = 0
        count_crystals = 0
        frame_stream = []
        header = 0
        indexing_methods = []
        
        event = ''

        stream = open(self.streamfile,'r')  # .readlines()

        for index, line in enumerate(stream):
            ### GEt header
            while (header == 0):
                self.header += line
                break

            ### Get beginning of an image
            if 'Begin chunk' in line:
                count_shots += 1
                append = 1
                header = 1
                # frame_stream.append(l)

            ### If
            if 'Image filename' in line: filename = line.split()[2]
            
            if 'Event: //' in line:
                try:
                    event = int(re.search(r'Event:\ \/\/(.+?)$',line).group(1))
                    #event = int(line.split("//")[-1])
                except AttributeError:
                    event = ''
            #else:
                #event = ''

            if 'indexed_by' in line and 'none' not in line:
                count_crystals += 1
                crystal = 1
                frame = Frame()
                frame.indexing = line.split()[2].strip()
                if frame.indexing not in indexing_methods:
                    indexing_methods.append(frame.indexing)
                frame.filename = filename
                frame.event    = event
                try:
                    f = os.path.split(filename)[1]
                    tag = os.path.splitext(f)[0].split('tag_')[1]
                    frame.timeline = tag
                except:
                    pass
            if 'diffraction_resolution_limit' in line:
                res = float(line.split()[5])
                frame.res = res

            if 'Cell parameters' in line:
                a0, b0, c0 = line.split()[2:5]
                frame.a = float(a0)
                frame.b = float(b0)
                frame.c = float(c0)
                alpha0, beta0, gamma0 = line.split()[6:9]
                frame.alpha = float(alpha0)
                frame.beta = float(beta0)
                frame.gamma = float(gamma0)

            if "End chunk" in line:
                if crystal == 1:
                    frame_stream.append(line)
                    frame.all = frame_stream
                    self.frames.append(frame)
                append = 0
                frame_stream = []
                crystal = 0

            if append == 1: frame_stream.append(line)

            if count_shots % 1000 == 0: 
                print('%7i frames parsed, %7i indexed frames found' % (count_shots, count_crystals), end='\r')
                #sys.stdout.flush()

        #print('%7i frames parsed, %7i indexed frames found\r' % (count_shots, count_crystals))
        #sys.stdout.flush()
        self.images           = count_shots
        self.indexed_images   = count_crystals
        self.indexing_methods = indexing_methods
        #self.indexed_images   = len(self.frames)
        #print 'indexing methods:',indexing_methods
        
        if 'Begin chunk' in self.header:
            self.header = re.sub('----- Begin chunk -----\n','',self.header) #quick and dirty removing the 'begin chunk' line that is added to the header
        
    
    def get_index_rate(self):
        """
        returns indexing rate = total number of images in stream / indexed images in stream
        """
        return float(self.indexed_images)/self.images
    
    def get_indexing_per_method(self):
        """
        returns the number of indexed images per indexing method
        """
        d = {}
        for method in self.indexing_methods:
            i = 0
            for f in self.frames:
                if f.indexing == method:
                    i+=1
            d[method]=i
        return d
    
    def get_total_number_of_cystals(self):
        """
        return the amount of crystals. This can be larger than the number of indexed images when multiple crystals were identified on a single image, but can never be smaller than the number of indexed images.
        """
        i = 0
        for f in self.frames:
            i+= self.get_number_of_indexed_crystals_in_frame(f.all)
        return i
    
    def get_stream_summary(self):
        """
        Print some stats about the indexing
        """
        print(("number of processed images: %d" %(self.images)))
        print(("number of indexed images: %d" %(self.indexed_images)))
        d = self.get_indexing_per_method()
        for method in d:
            print(("   %s: %d" %(method, d[method])))
        print(("Indexing rate: %.4f" %(self.get_index_rate())))
        print(("Number of unindexed images: %d" %(self.images - self.indexed_images)))
        print(("Number of crystals: %d" %(self.get_total_number_of_cystals())))
    
    def get_cell_stats(self):
        """
        get statistics on cell_parameters
        """
        aas = np.array([f.a for f in self.frames])
        bbs = np.array([f.b for f in self.frames])
        ccs = np.array([f.c for f in self.frames])
        
        aas_av = np.average(aas)
        aas_stdev = np.std(aas)
        bbs_av = np.average(bbs)
        bbs_stdev = np.std(bbs)
        ccs_av = np.average(ccs)
        ccs_stdev = np.std(ccs)
    
        return aas_av, aas_stdev, bbs_av, bbs_stdev, ccs_av, ccs_stdev
    
    def get_score(self):
        """
        score = indexrate/product-of-stds-on-axis
        """
        rate = self.get_index_rate()
        
        _, aas_stdev, _, bbs_stdev, _, ccs_stdev = self.get_cell_stats()
        stdev_product = aas_stdev * bbs_stdev * ccs_stdev
        
        return rate / stdev_product

    def get_number_of_indexed_crystals_in_frame(self, frame_all):
        """
        get the number of indexed crystals in a frame. Can differ from number of indexed imahges in case the
        --multi keyword was used during indexing
        """
        #return frame.all.count('--- Begin crystal\n') #depends too much on exact spelling
        return len([line for line in frame_all if 'Begin crystal' in line])
               

class Frame():
    def init(self):
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 90.
        self.beta = 90.
        self.gamma = 90.
        self.res = 5.
        self.index = 'none'
        self.filename = 'example.h5'
        self.event = ''
        self.timeline = 0
        self.indexing = ''



