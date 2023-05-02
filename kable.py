'''
MY USAGE
 - Mini set - 
 cd /Users/hendrixj/Dropbox/subDrop/development
 python kable_v0.py build -m toy_mani.tsv -k 5

 - Morpho set - 
cd /Users/hendrixj/Documents/morph/test_kable/inputs_big
python /Users/hendrixj/Dropbox/subDrop/development/kable_v0.py build -m manifest_all.tsv -k 24 -o out_v13_1

KNOWN PROBLEMS
If k-mer size is too small, search might go in circles

'''
import Bio
import argparse
import os
#import gfapy
#import graphviz #conda install python-graphviz pydot
import logging
import pandas as pd # conda install pandas
import subprocess
import sys

from Bio import Align
from Bio.Align import MultipleSeqAlignment
from datetime import datetime
from Levenshtein import distance as lev


# INITIATE LOGS
LOG = logging.getLogger('log_file')

class Card:
    ''' kmer from chopped sequence '''

    def __init__(self, seq, genome, contig, color, pos):
        self.seq = seq
        self.genome = genome
        self.contig = contig
        self.color = color
        self.start = pos
        self.features = []

    def add_feature(self, feat):
        self.features.append(feat)

class Feature:

    def __init__(self, name, ftype, fid, category, loc, color, start=None, source=None):
        self.name       = name      # ex. mmpL4 (gene name)
        self.type       = ftype     # ex. m vs. h
        self.id         = fid       # ex. mod_01
        self.category   = category  # ex. mod
        self.location   = loc
        self.color      = color
        self.kstart     = start
        #self.mod_seq   =   # ex. AATGmATG
        self.source     = source

class Annotation:

    def __init__(self, fname, fid, metaseq, strand, pos, category, source):
        self.name       = fname
        self.fid        = fid
        self.metaseq    = metaseq
        self.strand     = strand
        self.pos        = pos
        self.category   = category
        self.source     = source

class Dir:
    """ Base class for system directories """

    def __init__(self, path):
        self._path = None
        self.path = path

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, value):
        if not os.path.isabs(value):
            value = os.path.join(os.getcwd(), value)
        if os.path.isdir(value):
            self._path = value
        else:
            raise NotADirectoryError(value)

    @property
    def dirname(self):
        return self.path.strip("/").split("/")[-1]

    @property
    def children(self):
        children = [Dir(os.path.join(self.path, subdir)) 
            for subdir in os.listdir(self.path) 
            if os.path.isdir(os.path.join(self.path, subdir))]
        if len(children) > 0:
            return children
        else:
            return None

    @property
    def files(self):
        files = [File(os.path.join(self.path, file))
            for file in os.listdir(self.path)
            if os.path.isfile(os.path.join(self.path, file))]
        if len(files) > 0:
            return files
        else:
            return None

    def join(self, *args):
        return os.path.join(self.path, *args)

    def make_subdir(self, *args):
        """ Makes recursive subdirectories from 'os.path.join' like arguments """
        subdir = self.join(*args)
        return self.make(subdir)

    @classmethod
    def make(cls, path):
        try:
            os.makedirs(path)
            return cls(path)
        except FileExistsError:
            return cls(path)

    def __repr__(self):
        return self.path

class File:
    """ Base class for all file-types """

    def __init__(self, path, file_type=None):
        self._path = None
        self.path = path
        self.file_type = file_type

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, value):
        if not os.path.isabs(value):
            value = os.path.join(os.getcwd(), value)
        if os.path.isfile(value):
            self._path = value
        else:
            raise FileNotFoundError(value)

    #@property
    #def dir(self):
    #   return Dir(os.path.dirname(self.path))

    @property
    def filename(self):
        return os.path.basename(self.path)

    @property
    def file_prefix(self):
        return self.filename.split(".")[0]

    @property
    def extension(self):
        return self.filename.split(".")[-1]

class Graph:
    ''' Graph Object '''

    def __init__(self, nodes, edges, starts, cards, genomes, colors):
        self.nodes = set(nodes)
        self.edges = edges
        self.starts = set(starts)
        self.cards = cards
        self.genomes = genomes
        self.colors = set()

    def add_node(self, node):
        self.nodes.add(node)

    def add_edge(self, edge):
        self.edges.append(edge)

    def find_starts(self, not_starts):
        self.starts = self.nodes - not_starts

    def set_starts(self, starts):
        self.starts = starts

    def add_card(self, key, card):
        if key not in self.cards:
            self.cards[key] = []
        if card not in self.cards[key]:
            self.cards[key].append(card)

    def add_color(self, color):
        self.colors.add(color)

    def add_colors(self, colors):
        for color in colors:
            self.colors.add(color)

    def attach_cards(self, cards):
        self.info = cards

############################################################
## Helper Functions for build and add
def checkEx(ex, ftype):
    ''' 
        Takes a filename and expected file type.
        Returns true if filename matches file type.
    '''
    if ftype == 'fasta':
        suffixList = ['fasta', 'faa', 'fa', 'fna']

    elif ftype == 'gff':
        suffixList = ['gff']

    elif ftype == 'text':
        suffixList = ['tsv', 'txt', 'csv']

    elif ftype == 'graph':
        suffixList = ['gfa']

    elif ftype == 'vcf':
        suffixList = ['vcf']

    # try ex.lower()
    return ex in suffixList
def read_manifest(args):

    LOG.info('READ MANIFEST FILE...')
    # READ IN MANIFEST
    try:
        f = File(args.manifest)
        f1 = open(args.manifest, 'r')
    except:
        print('ERROR: Could not open manifest file. Exit')
        LOG.error('ERROR: Could not open manifest file. Exit')
        exit()

    # Check manifest file type
    if not checkEx(f.extension, 'text'):
        print('ERROR: Manifest file must be a text file. Exit.')
        LOG.error('ERROR: Manifest file must be a text file. Exit.')
        f1.close()
        exit()

    # LOOP OVER SAMPLE LINES
    samples = []
    for l in f1:
        if l.startswith('#'): # User can ignore a line if needed
            continue
        l = l.strip()   # Remove newline
        line = l.split('\t')
        
        # CHECK THAT FILES ARE VALID
        seqFound = False
        s = []
        strand = 'fwd'
        for file in line:
            if file == '+':
                continue
            elif file == '-':
                strand = 'rvs'
                continue

            if not os.path.isfile(file):
                continue
            try:
                f = File(file)
            except:
                print('File ', file, ' could not be found. Skipping...')
                LOG.error('File ', file, ' could not be found. Skipping...')
                continue
            
            if checkEx(f.extension, 'fasta'):
                seq = f
                seqFound = True
            else:
                s.append(f)
        if seqFound:
            s.insert(0, seq)
            entry = [s, strand]
            samples.append(entry)

    f1.close()

    return samples
def get_sequence(ifile):
    """ Get sequences from file """
    
    try:
        f = File(ifile)
        name = f.file_prefix
        f1 = open(f.path, 'r')
    except:
        print('ERROR. Could not open input file ', ifile)
        LOG.error('ERROR. Could not open input file ', ifile)
        return ''

    s = ''
    seq = {}
    for l in f1:
        line = l.strip()

        # START OF NEW SEQUENCE
        if line.startswith('>'):
            # Check if a sequence should be stored
            if len(s) > 0:
                seq[seqname] = str(s)
            # Reset for next seq
            s = ''
            seqname = line.strip('>').split()[0]
        else:       # Add line to sequence
            s = s + line

    # Handle last element
    seq[seqname] = str(s)

    f1.close()
    return seq, name
def read_gff(ifile):
    """ Get annotations from a gff file """

    try:
        f = File(ifile)
        f1 = open(f.path, 'r')
    except:
        print('ERROR: Could not open .gff file')
        LOG.error('ERROR: Could not open .gff file')
        return []
 
    anno = {}
    #name = ''
    #locus = ''
    count = 0
    for l in f1:
        count = count + 1
        if count % 1000 == 0:
            print(count)
        if l.startswith('>'):
            break
        if l.startswith('#'):
            continue
        l = l.strip()               # Remove newline
        l = '\t'.join(l.split())    # Remove whitespace
        line = l.split('\t')

        # Remove duplicate annotations from Prokka output
        if line[1] == 'Prodigal:2.6':
            continue
        # Parse information line
        m = line[8]
        meta = m.split(';')
        locus = ''
        name = ''
        product = ''
        for t in meta:
            tag = t.split('=')
            if tag[0] == 'Name':        # Gene name
                name = tag[1]
            if tag[0] == 'ID':  # ID tag - NOTE: changed from locus_tag
                locus = tag[1]
            if tag[0] == 'product':
                product = tag[1]

        # Lable hypothetical
        if name == '':
            if product != '':
                name = product
            else:
                name = 'hypothetical protein'

        # Build entry
        contig = str(line[0])
        #entry = [str(line[3]), str(line[4]), name, name, locus, str(line[2])]
        ##      Start, stop, Name, Name, Locus, category, source
        ## in gff, the positon is inclusive (i.e. 3-3 vs. 3-4) https://agat.readthedocs.io/en/latest/gxf.html
        entry = [int(line[3]), int(line[4]), name, name, locus, str(line[2]), f.filename]

        if contig not in anno:
            anno[contig] = []
        anno[contig].append(entry)
        #name = ''
        #locus = ''

    f1.close()
    return anno
def read_vcf(ifile):
    """ Get annotations from a vcf file """

    try:
        f = File(ifile)
        f1 = open(f.path, 'r')
    except:
        print('ERROR: Could not open .vcf file')
        LOG.error('ERROR: Could not open .vcf file')
        return []

    anno = {}
    name = ''
    locus = ''
    count = 0
    for l in f1:
        count = count + 1
        if count % 100 == 0:
            print(count)
        if l.startswith('#'):
            continue
        line = l.split('\t')

        # Parse variant type
        ref = line[3]
        alt = line[4]
        lref = len(ref)
        lalt = len(alt)
        if lref == lalt:
            name = 'snp'
        elif lref > lalt:
            name = 'deletion'
        elif lref < lalt:
            name = 'insertion'
        else:
            name = 'NA'

        # Parse information
        contig = str(line[0])
        fstart = int(line[1])
        fstop = int(int(line[1])+lref-1) # Length of reference seq
        #locus = 'NA'
        category = 'var'

        # Set unique name
        locus = '_'.join(('var', str(count).zfill(5)))


        ##  entry = Start, stop, Name, Name, Locus, category, source
        entry = [fstart, fstop, name, name, locus, category, f.filename]

        if contig not in anno:
            anno[contig] = []
        anno[contig].append(entry)
        name = ''
        locus = ''

    f1.close()
    return anno
def chop(sequence, features, genome, contig, strand='fwd', k=3, cyclic=True):
    '''
        Divide sequence into k-mers
        Add meta data
    '''
    # TODO: Allow for additional of other metadata
    # NOTE: i is indexed at 0
    # NOTE: features in gff are indexed at 1 - stop is inclusive
    color = ':'.join((str(genome), str(contig)))

    kmers = []
    counter = 0
    for i in range(0, len(sequence)):
        kend = i+k-1

        counter = counter + 1
        if counter % 100000 == 0:
            LOG.info(f'... ... ... Position: {str(counter)}')

        kmer = sequence[i:i+k].strip()
        if len(kmer) < k:
            continue

        if len(kmer) < 2:   # k-mer - 1 must be greater than 1 bp
            continue
        feats = []

        # POSITION IN SHORT FORM
        ## location: #off,#on,#off

        # ADD FEATURES
        ## Format: type, id, general name, location
        toRemove = []

        for feat in features:
            fstart   = int(feat[0])-1   # Convert to 0-index
            fstop    = int(feat[1])-1   # Convert to 0-index (inclusive)
            ftype    = feat[2]          # ex. m vs. h
            fname    = feat[3]          # ex. mmpL4
            fid      = feat[4]          # ex. mod_01
            category = feat[5]          # ex. mod vs. gene
            source   = feat[6]          # ex. file.gff

            # If annotations are downstream of kmer - stop looking
            # Feature list is ordered by start position
            if kend < fstart:
                break

            # If the feature is upstream of kmer - remove from search
            # Remove element from array if past
            if fstop < i:   # if annotation ends before sequence window begins
                toRemove.append(feat)
                # TODO: keep early features for cyclic wrap-around

            
            # CHECK IF ANY FEATURES EXIST WITHIN KMER
            frange = (fstart, fstop)
            krange = (i, kend)
            # If feature does not overlap -> move on
            if not overlap(frange, krange):
                continue

            loc = ''

            # Annotation starts in middle of kmer
            if fstart > i:
                l = fstart - i
                loc = '0'*l
                lb = fstart
            else:
                lb = i

            # Annotation ends after kmer
            if kend <= fstop:
                on = '1'*((kend-lb)+1)
                loc = loc + on
            else: # Annotationends in middle of kmer
                on = '1'*((fstop-lb)+1)
                off = '0'*(kend-fstop)
                loc = loc + on + off
            
            location = encode_loc(loc)
            if strand == 'rvs':
                loc = location.split('-')
                newLoc = (loc[2], loc[1], loc[0])
                location = '-'.join(newLoc)
            newFeat = Feature(fname, ftype, fid, category, location, color, i+1, source)
            ## Note: start refers to start of kmer, not start of feature
            ### fstart would be misleading because some annnotations are hundreds of k-mers long
            ### Start (i+1) is indexed at 1
            feats.append(newFeat)   # add name of feature
            

        # Remove features that have already been past
        for feat in toRemove:
            features.remove(feat)

        # HANDLE END OF SEQUENCE - cyclic or not
        length = len(kmer)
        if length != k: # if k-mer cut off by end of sequence
            #features don't wrap around...
            if cyclic:
                kmer += sequence[:(k-length)]

        # CREATE INFORMATION CARD
        #color = ':'.join((str(genome), str(contig)))
        if strand == 'rvs':
            kmer = get_reverse_complement(kmer)
        color = ':'.join((genome, contig))
        kcard = Card(kmer, genome, contig, color, i+1)
        if len(feats) > 0:
            for feat in feats:
                kcard.add_feature(feat)
        kmers.append(kcard)
    return kmers, color
def save_gfa(args, edges, tag=''):
    ''' Saves dBG in gfa format '''
    # Currently does not require gfapy
    # Currently alternates between nodes and edges
    ## TODO: save list of nodes and list of edges then write as chunks
    ### make set of nodes - to avoid duplicate lines
    ### make set of edges -to avoid duplicate lines

    LOG.info('SAVING GRAPH TO GFA...')

    # CREATE OUTPUT
    try:
        outdir = '/'.join((args.output_path, args.output_directory))
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        prefix = '/'.join((outdir, args.savename))
        if len(tag) > 0:
            prefix = '_'.join((prefix, tag))
        outf = "_".join((prefix, 'graph.gfa'))
        f1 = open(outf, 'w')

    except:
        print('ERROR: Could not configure output GFA file. Skipping...')
        LOG.error('ERROR: Could not configure output GFA file. Skipping...')
        return

    # ITERATE OVER NODES
    lookup = {}         # Save node names
    seg = 0             # To number nodes
    node_set = set()    # Remove duplicates
    edge_set = set()    # Remove duplicates
    for edge in edges:
        e1 = edge[0]
        e2 = edge[1]

        # NAME NODES
        ## Format: Record type, Name, Sequence, Optional
        # If node e1 has not been seen before...
        if e1 not in lookup:
            # Assign node name
            name = '_'.join(('node', str(seg).zfill(10)))
            lookup[e1] = name
            seg = seg + 1

            # Save node to file
            entry = ('S', name, e1)
            record = '\t'.join(entry) + '\n'
            node_set.add(record)

            # Record node ID for link (below)
            fromNode = name # (saves a lookup later)
        else:
            fromNode = lookup[e1]


        # If node e2 has not been seen before...
        if e2 not in lookup:
            # Assign node name
            name = '_'.join(('node', str(seg).zfill(10)))
            lookup[e2] = name
            seg = seg + 1

            # Save node to file
            entry = ('S', name, e2)
            record = '\t'.join(entry) + '\n'
            node_set.add(record)

            # Record node ID for link (below)
            toNode = name # (saves a lookup later)
        else:
            toNode = lookup[e2]

        # RECORD EDGES AS LINKS
        ## Format: Record Type, From, FromOrient, To, ToOrient, Overlap in CIGAR format, Optional
        ## NOTE: 1M indicates 1bp overlap
        entry = ('L', fromNode, '+', toNode, '+', '1M')
        record = '\t'.join(entry) + '\n'
        edge_set.add(record)

    # WRITE NODES
    for node in node_set:
        f1.write(node)

    # WRITE EDGES
    for edge in edge_set:
        f1.write(edge)

    f1.close()
def save_features(args, info, tag=''):
    
    LOG.info('SAVING FEATURES TO TEXT')

    # CREATE OUTPUT
    try:
        outdir = '/'.join((args.output_path, args.output_directory))
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        prefix = '/'.join((outdir, args.savename))
        if len(tag) > 0:
            prefix = '_'.join((prefix, tag))

        outf = "_".join((prefix, 'feats.tsv'))

        f1 = open(outf, 'w')

    except:
        print('ERROR: Could not configure output for feature data. Skipping...')
        LOG.error('ERROR: Could not configure output for feature data. Skipping...')
        return

    for key in info:
        for card in info[key]:
            fs = []
            for feat in card.features:
                newFeat = (feat.category, feat.id, feat.type, feat.name, feat.location, feat.source)
                newFeat = ','.join(newFeat)
                fs.append(newFeat)
            newFeats = ';'.join((fs))
            entry = (key, card.genome, card.contig, str(card.start), newFeats)
            record = '\t'.join(entry) + '\n'
            f1.write(record)

    f1.close()
def add_to_graph(G, cards):
    '''
        Takes list of kmer cards
        Adds kmers to graph
        Saves metadata in dictionary (info)
    '''

    not_starts = set()

    # NOTE: kmers is a list, not a set - there are duplicates

    # ITERATE OVER CARD INSTANCES FOR THAT KMER
    for card in cards:
        k1 = card.seq[:-1]  # (k-1)-mer
        k2 = card.seq[1:]   # (k-1)-mer

        # ADD NODES AND EDGES TO GRAPH
        G.add_node(k1)
        G.add_node(k2)
        G.add_edge((k1, k2))
        not_starts.add(k2)

        # FILE FEATURES
        seq = card.seq
        G.add_card(seq, card)

    # Remove nodes that have a node pointing to it
    G.find_starts(not_starts)

    return G
def add_to_graph_revC(G, kmers):
    '''
        Takes list of kmer cards
        Adds kmers to graph
        Saves metadata in dictionary (info)
    '''

    not_starts = set()

    # NOTE: kmers is a list, not a set - there are duplicates

    # ITERATE OVER CARD INSTANCES FOR THAT KMER
    for kmer in kmers:
        k1 = kmer.seq[:-1]  # (k-1)-mer
        k2 = kmer.seq[1:]   # (k-1)-mer

        # ADD NODES AND EDGES TO GRAPH
        G.add_node(k1)
        G.add_node(k2)
        G.add_edge((k2, k1))
        not_starts.add(k2)

        # FILE FEATURES
        seq = kmer.seq
        G.add_card(seq, kmer)

    # Remove nodes that have a node pointing to it
    G.find_starts(not_starts)

    return G

############################################################
## Helper Functions for Graphs
def initialize():
    nodes = set()
    edges = []
    starts = list()
    info = {}
    genomes = set()
    colors = set()

    return Graph(nodes, edges, starts, info, genomes, colors)
def find_starts(nodes, edges):
    not_starts = set()

    # RECORD NODES WITH AN IN-DEGREE
    for edge in edges:
        e2 = edge[1]
        not_starts.add(e2)

    # FIND NODES WITH 0 IN-DEGREES
    starts = list(nodes - not_starts)

    return starts
def read_graph(args):
    
    LOG.info('READING GRAPH...')
    # READ GRAPH
    try:
        ifile = File(args.input_graph)
    except:
        print('ERROR: Could not open graph file. Exit.')
        LOG.error('ERROR: Could not open graph file. Exit.')
        exit()
    if not checkEx(ifile.extension, 'graph'):
        print('ERROR: Graph file must be of type .gfa. Exit.')
        LOG.error('ERROR: Graph file must be of type .gfa. Exit.')
        exit()

    nodes = set()
    edges = []

    f1 = open(ifile.path, 'r')
    lookup = {}
    notfound = []
    counter = 0
    for l in f1:
        counter = counter + 1
        if counter % 1000000 == 0:
            LOG.info(f'\t... line: {str(counter)}')
        line = l.strip().split('\t')

        # HANDLE NODES (segments)
        if line[0] == 'S':
            name = line[1]
            seq = line[2]
            lookup[name] = seq

            # ADD NODE TO GRAPH
            nodes.add(seq)

        elif line[0] == 'L':
            # Currently assuming that everything is positive strand
            # Also assuming that all have 1M overlap
            name1 = line[1]
            name2 = line[3]

            # LOOKUP SEQUENCE FOR NODE WITH NAME
            if name1 in lookup and name2 in lookup:
                seq1 = lookup[name1]
                seq2 = lookup[name2]
                # Create edge
                edge = (seq1, seq2)
                edges.append(edge)
            else:
                pair = [name1, name2]
                notfound.append(pair)

    # HANDLE NOT FOUND
    ## In these edge pairs, the edge was found before the node. Attempt to resolve
    ## Do not save an edge node that does not have a sequence
    for pair in notfound:
        name1 = pair[0]
        name2 = pair[1]
        if name1 in lookup and name2 in lookup:
            seq1 = lookup[name1]
            seq2 = lookup[name2]
            # Create edge
            edge = (seq1, seq2)
            edges.append(edge)
        else:
            # TODO: Report failed reads to log?
            continue


    # IDENTIFY NODES WITH 0 IN-DEGREES
    starts = find_starts(nodes, edges)

    f1.close()

    return nodes, edges, starts
def read_feats(args, include_feats = True):
    
    LOG.info('READING FEATURES...')
    try:
        ifile = File(args.input_features)
    except:
        print('ERROR: Could not open features file. Exit.')
        LOG.error('ERROR: Could not open features file. Exit.')
        exit()

    info = {}
    cDict = {}
    genomes = set()
    colors = set()
    f = open(ifile.path, 'r')
    counter = 0
    #print('READING FEATS')
    for l in f:
        counter = counter + 1
        if counter % 1000000 == 0:
            LOG.info(f'\t... line: {str(counter)}')
        if l.startswith('#'):
            continue
        line = l.strip().split('\t')

        kmer = line[0]
        genome = line[1]
        contig = line[2]
        color = ':'.join((genome, contig))
        #print(color)
        pos = int(line[3])

        # CREATE INFORMATION CARD
        kcard = Card(kmer, genome, contig, color, pos)
        if len(line) > 4 and include_feats:
            feats = line[4].split(';')
            # CREATE FEATURE
            for feat in feats:
                entry = feat.split(',')
                newFeat = Feature(entry[3], entry[2], entry[1], entry[0], entry[4], color, pos, entry[5])
                kcard.add_feature(newFeat)

        # ADD CARD TO DICTIONARY
        if kmer not in info:
            info[kmer] = []
        info[kmer].append(kcard)

        # RECORD COLOR
        genomes.add(genome)
        colors.add(color)

        # COUNT KMER PRESENCE
        if kmer not in cDict:
            cDict[kmer] = {}
        if genome not in cDict[kmer]:
            cDict[kmer][genome] = 0
        cDict[kmer][genome] = cDict[kmer][genome] + 1
    f.close()

    list_genomes = list(genomes)
    list_genomes.sort()
    list_colors = list(colors)
    list_colors.sort()
    return(info, cDict, list_genomes, list_colors)

############################################################
## Helper Functions for Tasks
def encode_loc(code):
    
    left = 0
    feature = 0
    right = 0

    featSeen = False
    for i in code:
        if i == '0':
            if featSeen:
                right += 1
            else:
                left += 1
        else:
            feature += 1
            featSeen = True

    counts = (str(left), str(feature), str(right))
    c = '-'.join(counts)

    return c
def decode_loc(c, desc = None):
    code    = c.split('-')

    off1    = '0'*int(code[0])      # Leading 0s
    if desc is None:
        char = '1'
    else:
        if desc == '5mC':
            char = 'm'
        elif desc == '5hmC':
            char = 'h'
    on      = char*int(code[1])
    off2    = '0'*int(code[2])      # Trailing 0s

    # MAKE EXPANDED CODE
    decoded = (off1, on, off2)
    txt = ''.join(decoded)
    
    return txt
def overlap(a, b):
    ''' Determine if two ranges overlap '''
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]
def sort_dictionary(idict):
    ''' 
        Consolodiates keys and
        Sorts a dictionary by key 
    '''

    midDict = {}
    for elem in idict:
        for key in elem:
            for record in elem[key]:
                if key not in midDict:
                    midDict[key] = []
                midDict[key].append(record)

    sdict = {}
    for key in midDict:
        #print(midDict[key])
        sdict[key] = sorted(midDict[key])

    return sdict
def get_genome_list(list_colors):
    genomes = set()
    for color in list_colors:
        gen = color.split(':')[0]
        genomes.add(gen)
    list_genomes = list(genomes)
    list_genomes.sort()
    return list_genomes

############################################################
## Primary Functions
def build(args, command):

    configure(args, 'build', command)

    LOG.info('ENTER BUILD MODULE')
    # GET INPUT FILES
    samples = read_manifest(args)
    if samples is None:
        print('ERROR: Could not find any files. Exit')
        LOG.error('ERROR: Could not find any files. Exit')
        exit()

    
    # INITIAL K-MER SIZE CHECK - prevent too small or negative values
    if args.kmer_size < 3:
        k_size = 3
    else:
        k_size = args.kmer_size

    LOG.info(f'K-mer size = {k_size}')

    # INITIALIZE GRAPH
    G = initialize()

    # ADD SEQUENCES AND META DATA TO GRAPH
    LOG.info('ADD SEQUENCES AND METADATA TO GRAPH')
    for entry in samples:
        # COLLECT DATA FOR SAMPLES - sequence, strand, annotation
        annotations = []
        samp = entry[0]
        strand = entry[1]

        for f in samp:
            # INPUT FAST FILE
            if checkEx(f.extension, 'fasta'):       # Sequence file
                seq, genome = get_sequence(samp[0].path)
                LOG.info(f'. GENOME FILE: {genome}')
            
            # INPUT GFF ANNOTATION RECORDS
            if checkEx(f.extension, 'gff'):     # Annotation file
                LOG.info(f'... Reading annotation file: {f.filename}')
                anno = read_gff(f.path)
                annotations.append(anno)
                LOG.info('... Done.')

            # INPUT VCF ANNOTATION RECORDS
            if checkEx(f.extension, 'vcf'):     # variant file
                LOG.info(f'... Reading variant file: {f.filename}')
                var = read_vcf(f.path)
                annotations.append(var)
                LOG.info('... Done.')


        # SORT DICTIONARY - order annotations by start position
        sanno = sort_dictionary(annotations)

        # ADD CONTIGS TO GRAPH
        LOG.info('... Adding k-mers to graph')
        list_colors = []
        contigs = []
        for key in seq:
            contigs.append(key) # Add contig to list
            if key in sanno:
                annotations = sanno[key]
            else:
                annotations = []

            # DIVIDE SEQUENCE INTO KMERS - and add annotations
            LOG.info(f'... ... Finding k-mers...')
            kmers, color = chop(seq[key], annotations, genome, key, strand, k_size, cyclic=False)
            #print(len(kmers))
            # ADD KMERS TO GRAPH
            LOG.info(f'... ... ... Adding k-mers to graph...')
            if strand == 'fwd':
                G = add_to_graph(G, kmers)
            else:
                G = add_to_graph_revC(G, kmers)
            G.add_color(color)
            list_colors.append(color)


    LOG.info('DONE BUILDING GRAPH')
    # SAVE GRAPH AND FEATURES
    save_gfa(args, G.edges)         # Save graph
    save_features(args, G.cards)    # Save feature data
def add(args, command):

    # CONFIGURE
    configure(args, 'add', command)
    LOG.info('USING EXISTING OBJECTS:')
    LOG.info(f'\tGraph: {args.input_graph}')
    LOG.info(f'\tFeatures: {args.input_features}')

    # GET GRAPH
    nodes, edges, starts = read_graph(args)                 # Read graph file
    if args.exclude_features:
        info, cDict, list_genomes, list_colors = read_feats(args, False)    # Read features file - exclude features
    else:
        info, cDict, list_genomes, list_colors = read_feats(args)           # Read features file
    G = Graph(nodes, edges, starts, info, list_genomes, list_colors)    # Create graph

    # DEDUCE K-MER SIZE
    k_size = len(G.edges[0][0])+1

    # READ MANIFEST FOR ADD ONS
    samples = read_manifest(args)
    if samples is None:
        print('ERROR: Could not find any files. Exit')
        LOG.error('ERROR: Could not find any files. Exit')
        exit()

    # ADD SEQUENCES AND META DATA TO GRAPH
    LOG.info('ADD SEQUENCES AND METADATA TO GRAPH')
    for entry in samples:
        # COLLECT DATA FOR SAMPLES - sequence, strand, annotation
        annotations = []
        samp = entry[0]
        strand = entry[1]

        for f in samp:
            # INPUT FAST FILE
            if checkEx(f.extension, 'fasta'):       # Sequence file
                seq, genome = get_sequence(samp[0].path)
                LOG.info(f'. GENOME FILE: {genome}')
            
            # INPUT GFF ANNOTATION RECORDS
            if checkEx(f.extension, 'gff'):     # Annotation file
                LOG.info(f'... Reading annotation file: {f.filename}')
                anno = read_gff(f.path)
                annotations.append(anno)
                LOG.info('... Done.')

            # INPUT VCF ANNOTATION RECORDS
            if checkEx(f.extension, 'vcf'):     # variant file
                LOG.info(f'... Reading variant file: {f.filename}')
                var = read_vcf(f.path)
                annotations.append(var)
                LOG.info('... Done.')


        # SORT DICTIONARY - order annotations by start position
        sanno = sort_dictionary(annotations)

        # ADD CONTIGS TO GRAPH
        LOG.info('... Adding k-mers to graph')
        list_colors = []
        contigs = []
        for key in seq:
            contigs.append(key) # Add contig to list
            if key in sanno:
                annotations = sanno[key]
            else:
                annotations = []

            # DIVIDE SEQUENCE INTO KMERS - and add annotations
            LOG.info(f'... ... Finding k-mers...')
            kmers, color = chop(seq[key], annotations, genome, key, strand, k_size, cyclic=False)
            #print(len(kmers))
            # ADD KMERS TO GRAPH
            LOG.info(f'... ... ... Adding k-mers to graph...')
            if strand == 'fwd':
                G = add_to_graph(G, kmers)
            else:
                G = add_to_graph_revC(G, kmers)
            G.add_color(color)
            list_colors.append(color)
    #list_samples = []   # Track samples in graph
    

    LOG.info('DONE BUILDING GRAPH')
    # SAVE GRAPH AND FEATURES
    save_gfa(args, G.edges)         # Save graph
    save_features(args, G.cards)    # Save feature data


def find_vars(args, command):

    # CONFIGURE
    configure(args, 'find_vars', command)
    LOG.info('USING EXISTING OBJECTS:')
    LOG.info(f'\tGraph: {args.input_graph}')
    LOG.info(f'\tFeatures: {args.input_features}')

    # GET GRAPH
    nodes, edges, starts = read_graph(args)     # Read graph file
    if args.exclude_features:
        print('exclude')
        info, cDict, list_genomes, list_colors = read_feats(args, False)    # Read features file - exclude features
    else:
        info, cDict, list_genomes, list_colors = read_feats(args)           # Read features file
    G = Graph(nodes, edges, starts, info, list_genomes, list_colors)# Create graph

    # DEDUCE K-MER SIZE
    k_size = len(G.edges[0][0])+1

    # FIND GRAPH STARTING POINTS
    starts = list(G.starts)
    limit = args.max_depth

    RG = remove_exact(args, G, list_colors, cDict)

    # SAVE INTERMEDIATE GRAPH AND FEATURES
    if args.write_intermediates:
        save_gfa(args, RG.edges, 'unique')          # Save graph
        save_features(args, RG.cards, 'unique')     # Save feature data

    # FIND PATHS THROUGH GRAPH
    paths = get_paths_inexact(RG, args.reverse_complement, args.num_snps)

    # COMPARE SEQUENCES AND FEATURES
    LOG.info(f'COMPARING VARIANT SEQUENCES...')
    alignments = get_alignment(RG, paths, k_size, list_colors, args.include_mods)
    variants = get_vars(args, list_colors, alignments, 'snps')

    # GENERATE EXPORTS
    job = 'snp_vars'
    LOG.info(f'WRITING OUTPUT for {job}....')
    export_vcf(args, list_colors, variants, '', job)            # variants
    export_alignments(args, list_colors, alignments, job)       # pseudo alignment
    
    # limited meaning for var search
    export_sum_table(args, alignments, list_colors, job)

    # Does not have much meaning for var search
    #export_mod_table(args, list_colors, alignments, job)   # table of base modifications






def query(args, command):
    ''' Finds all instances of given sequence
        Exact: No variants allowed
        Bubble: TODO 
    '''
    print('not currently opperational')
    exit()

    # CONFIGURE
    configure(args, 'query', command)
    LOG.info('USING EXISTING OBJECTS:')
    LOG.info(f'\tGraph: {args.input_graph}')
    LOG.info(f'\tFeatures: {args.input_features}')

    # GET GRAPH
    nodes, edges, starts = read_graph(args)     # Read graph file
    info, list_colors = read_feats(args)        # Read features file
    G = Graph(nodes, edges, starts, info, list_colors)      # Create graph

    # DEDUCE K-MER SIZE
    k_size = len(G.edges[0][0])+1

    # DIVIDE QUERY INTO K-MERS
    path_set = []
    query = args.query
    kmers = qchop(query, k_size, args.is_circular)
    path_set.append(kmers)
    if args.reverse:
        rcquery = get_reverse_complement(query)
        rckmers = qchop(rcquery, k_size, args.is_circular)
        path_set.append(rckmers)

        kmers = kmers + rckmers

    # FETCH CARDS RELATED TO KMERS
    LOG.info('FETCH CARDS')
    exact = True
    deck = []
    p = []
    if exact:
        for path in path_set:
            # Copy kmers array
            tmp = []
            for kmer in path:
                tmp.append(kmer)
            for kmer in tmp:
                # If kmer not in graph, remove
                if kmer not in G.cards:
                    path.remove(kmer)
                    continue

                for card in G.cards[kmer]:
                    deck.append(card)
            p.append(path)
        paths = [p]
    else:
        print('TODO')
        # Find all paths with up to n SNPs or m indels

    # Temp - testing
    #colors.add('toy:seq3')
    #paths = [['ATGCT', 'TGCTC'], ['ATGGT', 'TGGTC']]

    # ORGANIZE COLORS
    #list_colors = list(colors)
    #list_colors.sort()
    
    print_alignments(args, G, paths, k_size, list_colors, 'query')

def mod_search(args, command):

    # CONFIGURE
    configure(args, 'mod_search', command)
    LOG.info('USING EXISTING OBJECTS:')
    LOG.info(f'\tGraph: {args.input_graph}')
    LOG.info(f'\tFeatures: {args.input_features}')

    # GET GRAPH
    nodes, edges, starts = read_graph(args)                 # Read graph file
    info, cDict, list_genomes, list_colors = read_feats(args)       # Read features file
    G = Graph(nodes, edges, starts, info, list_genomes, list_colors)    # Create graph

    # DEDUCE K-MER SIZE
    k_size = len(G.edges[0][0])+1   

    # FIND KMERS WITH AT LEAST ONE MODIFICATION
    LOG.info('FIND KMERS WITH MODIFICATIONS....')
    mod_deck = []
    ghost_deck = []
    for kmer in G.cards:
        foundMod = find_mod(G.cards[kmer])
        if foundMod:
            for card in G.cards[kmer]:
                mod_deck.append(card)

            if args.reverse_complement:
                krc = get_reverse_complement(kmer)
                try:
                    for card in G.cards[krc]:
                        ghost_deck.append(card)
                except:
                    continue

    # MAKE SUBGRAPH FOR MODIFICATIONS
    LOG.info('MAKE SUBGRAPH....')
    SG = initialize()
    SG = add_to_graph(SG, mod_deck)
    starts = list(SG.starts)
    if args.reverse_complement:
        SG = add_to_graph(SG, ghost_deck)

    # SAVE INTERMEDIATE GRAPH AND FEATURES
    if args.write_intermediates:
        LOG.info('WRITING INTERMEDIARY FILES....')
        save_gfa(args, SG.edges, 'mod')         # Save graph
        save_features(args, SG.cards, 'mod')    # Save feature data

    # FIND PATHS IN SUBGRAPH
    paths = get_paths(SG, starts, args.reverse_complement)

    
    # COMPARE SEQUENCES AND FEATURES
    LOG.info(f'COMPARING MODIFIED SEQUENCES...')
    alignments = get_alignment(SG, paths, k_size, list_colors, True)
    
    # GENERATE ALIGNMENT EXPORTS
    job = 'mod_vars'
    LOG.info(f'WRITING OUTPUT for {job}....')
    export_alignments(args, list_colors, alignments, job)           # pseudo alignment
    export_mod_table(args, list_colors, alignments, 'base_mods')    # table of base modifications
    export_sum_table(args, alignments, list_colors, job)

    # EVALUATE VARIANTS
    variants = get_vars(args, list_colors, alignments, 'mods')
    export_vcf(args, list_colors, variants, '.m', job)              # variants



def bubble_query(args, command):
    print('not currently opperational')
    exit()

    # CONFIGURE
    configure(args, 'bubble_query', command)
    LOG.info('USING EXISTING OBJECTS:')
    LOG.info(f'\tGraph: {args.input_graph}')
    LOG.info(f'\tFeatures: {args.input_features}')

    # GET GRAPH
    nodes, edges, starts = read_graph(args)     # Read graph file
    info, list_colors = read_feats(args)        # Read features file
    G = Graph(nodes, edges, starts, info, list_colors)      # Create graph

    # DEDUCE K-MER SIZE
    k_size = len(G.edges[0][0])+1

    # DIVIDE QUERY INTO K-MERS
    path_set = []
    query = args.query
    kmers = qchop(query, k_size, False)
    pf = kmers[0]
    pr = kmers[-1]
    starts = [pf[:-1]]
    stops = [pr[1:]]
    path_set.append(kmers)
    if args.reverse:
        rpf = get_reverse_complement(pr)
        starts.append(rpf[:-1])
        rpr = get_reverse_complement(pf)
        stops.append(rpr[1:])

    LOG.info(f'Starts: {starts} \nstops: {stops}')

    # CONVERT GRAPH TO DICTIONARY
    m = map(G.edges)

    # NOTE: G is undirected
    #s = set()          # set of found cycles
    #c = set()          # set of computed contigs

    # GET PATHS FOR START-STOP CONSTRAINTS
    paths = []
    #max_depth = (3*k_size) + 2             # Set max depth to search (arbitrary)
    max_depth = 980
    # Check max depth for given system
    ## import sys 
    ##print(sys.getrecursionlimit())
    LOG.info(f'Max depth: {str(max_depth)}')
    count = 0
    for i in range(0, len(starts)):
        count = count + 1
        LOG.info(f'SEARCH FOR PATTERN {str(count)}')
        vertex = starts[i]
        stop = stops[i]
        path_set = get_subgraph(m, vertex, stop, max_depth)
        paths.append(path_set)
        LOG.info(f'\tPaths so far: {paths}')


    # CONDENSE PATH TO PROPER K-MER LENGTH
    cpaths = []
    for path_set in paths:
        cpath_set = []
        for path in path_set:
            cpath_set.append(condense_path(path))
        cpaths.append(cpath_set)

    LOG.info(f'Condensed paths: {cpaths}')

    # PRINT RESULTS
    print_alignments(args, G, cpaths, k_size, list_colors, 'bubble')

    print('done')



def primer_search(args, command):
    print('not currently opperational')
    exit()

    # CONFIGURE
    configure(args, 'primer_search', command)
    LOG.info('USING EXISTING OBJECTS:')
    LOG.info(f'\tGraph: {args.input_graph}')
    LOG.info(f'\tFeatures: {args.input_features}')

    # GET GRAPH
    nodes, edges, starts = read_graph(args)     # Read graph file
    info, list_colors = read_feats(args)        # Read features file
    G = Graph(nodes, edges, starts, info, list_colors)      # Create graph

    # DEDUCE K-MER SIZE
    k_size = len(G.edges[0][0])+1

    # CONVERT GRAPH TO DICTIONARY
    m = map(G.edges)

    # HANDLE CASE WEHRE SEARCH SEQUENCE IS TOO SHORT
    if len(args.start_sequence) < k_size or len(args.end_sequence) < k_size:
        LOG.error('Search sequences must be at least the length of the k-mers. Exit.')
        exit()

    # DIVIDE QUERIES INTO K-MERS
    s_kmers = qchop(args.start_sequence, k_size-1, False)
    e_kmers = qchop(args.end_sequence, k_size-1, False)

    # Only care about inner sequences.... address later
    starts = []
    starts.append(s_kmers[-1])
    stops = []
    stops.append(e_kmers[0])

    LOG.info(f'Starts: {starts} \nstops: {stops}')
    print(G.nodes)

    # GET PATHS FOR START-STOP CONSTRAINTS
    paths = []
    #max_depth = (3*k_size) + 2             # Set max depth to search (arbitrary)
    max_depth = 980
    # Check max depth for given system
    ## import sys 
    ##print(sys.getrecursionlimit())
    LOG.info(f'Max depth: {str(max_depth)}')
    count = 0
    for i in range(0, len(starts)):
        count = count + 1
        LOG.info(f'SEARCH FOR PATTERN {str(count)}')
        vertex = starts[i]
        stop = stops[i]
        path_set = get_subgraph(m, vertex, stop, max_depth)
        paths.append(path_set)
        LOG.info(f'\tPaths so far: {paths}')


    # CONDENSE PATH TO PROPER K-MER LENGTH
    cpaths = []
    for path_set in paths:
        cpath_set = []
        for path in path_set:
            cpath_set.append(condense_path(path))
        cpaths.append(cpath_set)

    LOG.info(f'Condensed paths: {cpaths}')

    # PRINT RESULTS
    print_alignments(args, G, cpaths, k_size, list_colors, 'primer')

    print('done')


#############################################################

def configure(args, prog, command):
    ''' Configure basic aspects of program '''

    global BASEDIR, LOG

    # SET UP OUTPUT
    topdir = Dir(args.output_path)
    BASEDIR = topdir.make_subdir(args.output_directory)

    # INITIATE LOG FILE
    LOG.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    logname = '_'.join(('grapher', prog))
    logname = '.'.join((logname, 'log'))
    LOG_File = BASEDIR.join(logname)
    file_handler = logging.FileHandler(LOG_File)
    file_handler.setFormatter(formatter)

    LOG.addHandler(file_handler)

    LOG.info(f'RUNNING Kable')
    LOG.info(command)

def main(args):
    command = 'Command: %s' % ' '.join(sys.argv)
    print(f'Running: ', command)
    #configure(args, command)
    
    
    args.func(args, command)



if __name__== "__main__":

    cwd = os.getcwd()

    parser = argparse.ArgumentParser(description='program description')
    subparsers = parser.add_subparsers(dest="cmd", help='available actions')
    #subparsers = parser.add_subparsers(title="build", dest="Build a de bruijn graph", help='available actions')
    #subparsers = parser.add_subparsers(title="stat", dest="Get basic statistics", help='available actions')
    subparsers.required = True

    # PARSER : ROOT
    __version__ = "0.0.0"
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))


    

    # DEFINE SUBPARSERS
    parser_build = subparsers.add_parser('build')
    parser_build.set_defaults(func=build)

    parser_add = subparsers.add_parser('add')
    parser_add.set_defaults(func=add)

    parser_query = subparsers.add_parser('query')
    parser_query.set_defaults(func=query)

    parser_mod_search = subparsers.add_parser('mod_search')
    parser_mod_search.set_defaults(func=mod_search)

    parser_find_vars = subparsers.add_parser('find_vars')
    parser_find_vars.set_defaults(func=find_vars)

    parser_bubble_q = subparsers.add_parser('bubble_query')
    parser_bubble_q.set_defaults(func=bubble_query)

    #parser_stat = subparsers.add_parser('stat')#, parents=[parser])
    #parser_stat.set_defaults(func=stat)

    parser_primer_s = subparsers.add_parser('primer_search')
    parser_primer_s.set_defaults(func=primer_search)


    ''' 
    parser_search = subparsers.add_parser('search')#, parents=[parser])
    parser_search.set_defaults(func=search)
    '''
    #build = subparsers.add_parser('build', help='Build a de bruijn graph', parents=[parent_parser])
    
    # PARSR : BUILD
    parser_build.add_argument('-k', '--kmer_size', help='Length of kmer', default=23, type=int)
    parser_build.add_argument('-m', '--manifest', help='File specifying location of sequences and corresponding metadata')
    parser_build.add_argument('-n', '--savename', help='Name of output file', default='kable', type=str)
    parser_build.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_build.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
    #parser_build.add_argument('-r', '--reverse_complement', default=False, help='Align reverse complement sequences', action='store_true')
    #parser_build.add_argument('-s', '--offeset_mod', default=False, help='Adjust modification boundaries to not include end position', action='store_true')
    #parser_build.add_argument('-w', '--write_var_graph', default=False, help='Save graph and features after filtering', action='store_true')

    # PARSER : ADD
    parser_add.add_argument('-e', '--exclude_features', default=False, help='Annotation data will not be read from feature file', action='store_true')  
    parser_add.add_argument('-f', '--input_features', help='File containing feature data')
    parser_add.add_argument('-g', '--input_graph', help='File containing graph')
    parser_add.add_argument('-m', '--manifest', help='File specifying location of sequences and corresponding metadata')
    parser_add.add_argument('-n', '--savename', help='Name of output file', default='kable_add', type=str)
    parser_add.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_add.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)


    # PARSER : MOD_SEARCH
    #parser_mod_search.add_argument('-ag', '--score_gap', help='Gap score during alignment', default=1, type=int)
    #parser_mod_search.add_argument('-am', '--score_match', help='Match score during alignment', default=10, type=int)
    #parser_mod_search.add_argument('-as', '--score_snp', help='SNP score during alignment', default=5, type=int)
    #parser_mod_search.add_argument('-a', '--alignments', help='Top alignments to assess', default=5, type=int)
    parser_mod_search.add_argument('-d', '--max_depth', help='Maximum depth to search for paths', default=900, type=int)
    parser_mod_search.add_argument('-f', '--input_features', help='File containing feature data')
    parser_mod_search.add_argument('-g', '--input_graph', help='File containing graph')
    parser_mod_search.add_argument('-m', '--mods', help='Modifications to look for', default=['m', 'h'], nargs='+')
    parser_mod_search.add_argument('-n', '--savename', help='Name of output file', default='kable', type=str)
    parser_mod_search.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_mod_search.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
    parser_mod_search.add_argument('-r', '--reverse_complement', default=False, help='Align reverse complement sequences', action='store_true')
    parser_mod_search.add_argument('-w', '--write_intermediates', default=False, help='Save graph and features after filtering', action='store_true')


    # PARSER : MOD_SUM

    # PARSER : FIND VARS
    #parser_find_vars.add_argument('-ag', '--score_gap', help='Gap score during alignment', default=1, type=int)
    #parser_find_vars.add_argument('-am', '--score_match', help='Match score during alignment', default=10, type=int)
    #parser_find_vars.add_argument('-as', '--score_snp', help='SNP score during alignment', default=5, type=int)
    #parser_find_vars.add_argument('-a', '--alignments', help='Top alignments to assess', default=5, type=int)
    parser_find_vars.add_argument('-d', '--max_depth', help='Maximum depth to search for paths', default=900, type=int)
    parser_find_vars.add_argument('-e', '--exclude_features', help='Exlude all features (annotations, base modifications, etc.) in variant calling', default=False, action='store_true')
    parser_find_vars.add_argument('-f', '--input_features', help='File containing feature data')
    parser_find_vars.add_argument('-g', '--input_graph', help='File containing graph')
    parser_find_vars.add_argument('-m', '--include_mods', help='Include base modifications in variant calling', default=False, action='store_true')
    parser_find_vars.add_argument('-n', '--savename', help='Name of output file', default='kable', type=str)
    parser_find_vars.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_find_vars.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
    parser_find_vars.add_argument('-r', '--reverse_complement', default=False, help='Align reverse complement sequences', action='store_true')
    parser_find_vars.add_argument('-s', '--num_snps', default=2, help='Number of SNPs to allow', type=int)
    parser_find_vars.add_argument('-w', '--write_intermediates', default=False, help='Save graph and features after filtering', action='store_true')


    # PARSER : QUERY
    parser_query.add_argument('-c', '--is_circular', help='Query sequence is circular', default=False, action='store_true')
    parser_query.add_argument('-f', '--input_features', help='File containing feature data')
    parser_query.add_argument('-g', '--input_graph', help='File containing graph')
    parser_query.add_argument('-n', '--savename', help='Name of output file', default='kable', type=str)
    parser_query.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_query.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
    parser_query.add_argument('-q', '--query', help='Sequence to search for', required=True)
    parser_query.add_argument('-r', '--reverse', help='Search for reverse compliment', default=False, action='store_true')


    # PARSER : BUBBLE QUERY
    parser_bubble_q.add_argument('-f', '--input_features', help='File containing feature data')
    parser_bubble_q.add_argument('-g', '--input_graph', help='File containing graph')
    parser_bubble_q.add_argument('-n', '--savename', help='Name of output file', default='kable', type=str)
    parser_bubble_q.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_bubble_q.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
    parser_bubble_q.add_argument('-q', '--query', help='Sequence to search for', required=True)
    parser_bubble_q.add_argument('-r', '--reverse', help='Search for reverse compliment', default=False, action='store_true')


    # PARSER : PRIMER SEARCH
    parser_primer_s.add_argument('-f', '--input_features', help='File containing feature data')
    parser_primer_s.add_argument('-e', '--end_sequence', help='Sequence to end search', required=True)
    parser_primer_s.add_argument('-g', '--input_graph', help='File containing graph')
    parser_primer_s.add_argument('-n', '--savename', help='Name of output file', default='kable_primer', type=str)
    parser_primer_s.add_argument('-o', '--output_directory', help='Name of output directory', default='kable_output', type=str)
    parser_primer_s.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
    #parser_primer_s.add_argument('-q', '--query', help='Sequence to search for', required=True)
    parser_primer_s.add_argument('-s', '--start_sequence', help='Sequence to start search', required=True)


    args = parser.parse_args()

    main(args)






