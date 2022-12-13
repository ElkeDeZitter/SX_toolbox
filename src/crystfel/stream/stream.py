# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import re
import numpy as np
import random

class Stream(object):
    """
    Class that efficiently reads in a crystfel stream file and allow calculate statistics from the stream file as well as modifying it.
    
    Parameters
    ----------
    streamfile (str)
        CrystFEL stream file
        
    Output
    ----------
    stream object consisting out of frames.
    
    
    """

    def __init__(self, streamfile):
        self.streamfile = streamfile
        self.frames = []
        self.header = ''
        self.parse_stream()

    def parse_stream(self):
        """
        Read stream file line by line and dedicate lines and info to indexed frames and crystals
        """
        append_frame = 0
        append_crystal = 0
        indexed_crystal = 0
        count_shots = 0
        count_images = 0
        frame_stream = []
        crystal_stream = []
        header = 0
        indexing_methods = []
        end_chunk_line = []
        
        event = ''
        
        with open(self.streamfile,'r') as s:
            stream = s.readlines()

        for index, line in enumerate(stream):
            # GEt header
            while (header == 0):
                self.header += line
                break

            # Get beginning of an image
            if 'Begin chunk' in line:
                count_shots += 1
                append_frame = 1
                header = 1
                # frame_stream.append(l)

            ### If
            elif 'Image filename' in line: filename = line.split()[2]
            
            elif 'Event:' in line:
                try:
                    #Event should be number that succeeds "//" 
                    event = line.rstrip().lstrip().split("//")[-1]
                    if event:
                        try:
                            event = int(event)
                        except ValueError:
                            event = str(event)
                    else:
                        #Event can also be a tag
                        event = line.rstrip().lstrip().split("Event:")[-1]
                        if event:
                            try:
                                event = int(event)
                            except ValueError:
                                event = str(event)
                except AttributeError:
                    event = ''

            elif 'indexed_by' in line and 'none' not in line:
                count_images += 1
                indexed_crystal = 1
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
                
            elif 'Begin crystal' in line:
                append_frame = 0
                append_crystal = 1
                crystal = Crystal()
                
            elif 'diffraction_resolution_limit' in line:
                res = float(line.split()[5])
                crystal.res = res

            elif 'Cell parameters' in line:
                a0, b0, c0 = line.split()[2:5]
                crystal.a = float(a0)
                crystal.b = float(b0)
                crystal.c = float(c0)
                alpha0, beta0, gamma0 = line.split()[6:9]
                crystal.alpha = float(alpha0)
                crystal.beta = float(beta0)
                crystal.gamma = float(gamma0)
                
            elif 'End crystal' in line:
                if indexed_crystal == 1:
                    #not clear if this will be usefull or not
                    crystal.filename = frame.filename
                    crystal.event = frame.event
                    crystal.timeline = frame.timeline
                    crystal.indexing = frame.indexing

                    #Attribute the lines of the crystal to the crystal
                    crystal_stream = [line for line in crystal_stream if line != "\n"]
                    crystal.reflections = crystal_stream
                    
                    frame.crystals.append(crystal)
                append_crystal = 0
                crystal_stream = []

            elif "End chunk" in line:
                if indexed_crystal == 1:
                    #attribute the lines of the chunk to the frame
                    frame.head = frame_stream
                    self.frames.append(frame)
                frame_stream = []
                indexed_crystal = 0
                if not end_chunk_line:
                    end_chunk_line = line

            if append_frame == 1:
                frame_stream.append(line)
            
            if append_crystal == 1:
                crystal_stream.append(line)

            if count_shots % 1000 == 0: 
                print('%7i frames parsed, %7i indexed frames found' % (count_shots, count_images), end='\r')
                #sys.stdout.flush()

        self.end_chunk_line = end_chunk_line
        self.images = count_shots
        self.indexed_images = count_images
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
        return the amount of crystals. This can be larger than the number of indexed images when
        multiple crystals were identified on a single image, but can never be smaller than the number of indexed images.
        """
        return len([f.crystals for f in self.frames])
    
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
        print(("Number of crystals: %d" %(self.get_total_number_of_cystals())))
        print(("Number of unindexed images: %d" %(self.images - self.indexed_images)))
    
    def get_cell_stats(self):
        """
        get statistics on cell_parameters
        """
        
        aas = []
        bbs = []
        ccs = []
        for frame in self.frames:
            for crystal in frame.crystals:
                aas.append(crystal.a)
                bbs.append(crystal.b)
                ccs.append(crystal.c)
                
        aas = np.array(aas)
        bbs = np.array(bbs)
        ccs = np.array(ccs)
        
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


class Frame():
    """
    Frame object initiate an indexed image from a CrystFEL stream file.
    Attributes:
        #index = Not used
                filename = name of the image file (e.g. .h5 file). can be single-event or nulti-event file 
        event = event number. Will be left blanc for signle-event images
        timeline = tag. can be used in time-resolved experiments. Key-word will be replaced by tag.
        indexing = method that was used to index the frame
        head = all info that is listed before the crystal information (including peaks)
        crystals = list with crystal info.
    """
    
    def __init__(self):

        ##self.index = 'none' #Not used
        self.filename = 'example.h5'
        self.event = ''
        self.timeline = 0
        self.indexing = ''
        self.head = ['']
        self.crystals = []
        
class Crystal(Frame):
    """
    Subclass of frame is a crystal since a frame may contain multiple crystals.
    Attributes:
        a = unit cell a axis (in nm)
        b = unit cell b axis (in nm)
        c = unit cell c axis (in nm)
        alpha = unit cell alpha angle (in degr.)
        beta = unit cell beta angle (in degr.)
        gamma = unit cell gamma angle (in degr.)
        res = resolution (in A)
        reflections: h,k,l,I, sigma(I), peak, background, fs/px ss/px panel information
        
    """
    def __init__(self):
        #Not clear if heritage will be usefull or not
        super().__init__()
        
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 90.
        self.beta = 90.
        self.gamma = 90.
        self.res = 5.
        self.reflections = ''

        



