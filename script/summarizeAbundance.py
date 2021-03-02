#! /conda2/envs/humann3/bin/python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright 2018-2021, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to summarize abundance of a group like one pathway based on the abundance of each sub item.

<sub-item-abundance matrix>

ID	Samp1	Samp2	Samp3	Samp4
Probe1	1	1	1	1
Probe2	1	1	1	1
Probe3	1	1	1	1
Probe4	1	1	1	1
Probe5	1	1	1	1
Probe6	1	1	1	1

<groupfile>

1. First column should match first column of <sub-item abundance matrix>

ID	Symbol
Probe1	A
Probe2	B
Probe2	C
Probe4	D
Probe5	F
Probe6	A
Probe7	C
Probe4	F
Probe8	F
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import pandas as pd
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print(desc, file=sys.stderr)
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
         help="Sub-item abundance file with format specified above")
    parser.add_option("-m", "--map-file", dest="map_file",
         help="Map file containing group information")
    parser.add_option("-c", "--grpcolumn", dest="grp_column",
            default='2', 
             help="The column(s) contains group information. Multiple columns should be supplied as <2,3> represents the second and third column (1-based).")
    parser.add_option("-s", "--grpsep", dest="grp_sep",
            default=",", 
             help="Separtor(s) for each group. Default <,> represents each group would be split by comma. If each group using special separtors, they should be joined by one <+>. <,+*> represents group separators are <,> or any character (*) (each alphabet would be treat as one group) ")
    parser.add_option("-k", "--subitem-keep", dest="subitem_keep",
        default="all", 
         help="<unique> means only keep sub-items which map to only one group. <all> (default) means keep sub items map to all group.")
    parser.add_option("-e", "--abundance-keep", dest="abundance_keep",
        default="sum", 
         help="This parameter deals with groups containing multiple sub-items. <median>  means using the median abundance of all sub-items map to one group as final group abundance. <min> means using the minimum abundance of all sub-items map to one group as final group abundance. <sum> (default) means using the sum abundance of all sub-items map to one group as final group abundance.")
    parser.add_option("-n", "--norm-type", dest="norm_type",
        default = "raw", 
         help="Specify the output data type, accept <raw>, <cpm> or both <raw,cpm>.")
    parser.add_option("-o", "--output-prefix", dest="output_prefix",
         help="Output file prefix.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def norm_expr(probe_expr_df, norm_type):
    if norm_type == "raw":
        return probe_expr_df
    elif norm_type == "cpm":
        norm_factor = probe_expr_df.sum()
        return probe_expr_df * 1000000 / norm_factor
#----------------------------------------

def main():
    global debug
    options, args = cmdparameter(sys.argv)
    debug = options.debug
    #-----------------------------------
    probe_expr_file = options.filein
    probe_map_file = options.map_file
    subitem_keep    = options.subitem_keep
    abundance_keep     = options.abundance_keep
    grp_columnL     = [int(i)-1 for i in options.grp_column.split(',')]
    grp_column_len = len(grp_columnL)
    grp_sepL     = [i for i in options.grp_sep.split('+')]
    grp_sep_len = len(grp_sepL)

    norm_typeL     = [i for i in options.norm_type.split(',')]
    if grp_sep_len == 1:
        grp_sepL = grp_sepL * grp_column_len
    else:
        assert grp_sep_len == grp_column_len, "Unequal number of columns (-c) and separtors (-s)" 
    if debug:
        print(f'grp_columnL={grp_columnL}')
        print(f'grp_sepL={grp_sepL}')
    #if abundance_keep != "median":
    #    sys.stderr.write(f"{abundance_keep} currently unsupported for <-e>.\n")
    #    return
    output_prefix = options.output_prefix

    # abundance_keep_funD = {'median': median, 'sum': sum, 'minimum': minimum}
    # abundance_keep_fun = abundance_keep_funD[abundance_keep]

    verbose = options.verbose
    #-----------------------------------
    probe_expr_df = pd.read_csv(probe_expr_file, sep="\t", header=0, index_col=0)
    probe_expr_df.index = probe_expr_df.index.map(str)

    if debug:
        sys.stderr.write(f"Probes example in expression file: {probe_expr_df.index[:10]}.\n")
        

    if debug:
        sys.stderr.write(f"Read in probe expression matrix: {probe_expr_df.shape}.\n")

    # probe_expr_df = probe_expr_df.loc[(probe_expr_df > 0).any(axis=1)]
    
    #if debug:
    #    sys.stderr.write(f"Probe expression matrix after filtering unexpressed probes: {probe_expr_df.shape}.\n")

    abundanceL = []
    for i in norm_typeL:
        abundanceL.append(norm_expr(probe_expr_df, i))
    #--------------------------------------------------

    usecols = [0]
    usecols.extend(grp_columnL)
    #print(f"usecols={usecols}")
    probe_map_df = pd.read_csv(probe_map_file, sep="\t", header=0, index_col=None, dtype=str, usecols=usecols)
    #probe_map_df.index = probe_map_df.index.map(str)

    if debug:
        sys.stderr.write(f"Read in probe-map matrix: {probe_map_df.shape}.\n")
        sys.stderr.write(f"Read in probe-map matrix: {probe_map_df.head()}.\n")
    keyName = list(probe_map_df.columns)[0]
    grp_columnL = list(probe_map_df.columns)[1:]
    #print(f"grp_columnL={grp_columnL}")

    for grp, sep in zip(grp_columnL, grp_sepL):
        grpDF = probe_map_df.loc[probe_map_df[grp]!="", [keyName, grp]]
        for exprDF, exprType in zip(abundanceL, norm_typeL):
            #grpMergeDF = pd.concat([grpDF, exprDF], axis=1, sort=False, join="inner")
            grpMergeDF = grpDF.merge(exprDF, left_on=keyName, right_index =True, how="inner")
            if debug:
                sys.stderr.write(f"Read in probe-map matrix: {grpMergeDF.shape}.\n")
                sys.stderr.write(f"Read in probe-map matrix: {grpMergeDF.head()}.\n")
            
            if sep == "*":
                sep = ""
            grpMergeDF[grp] = grpMergeDF[grp].str.split(sep)
            grpMergeDF_explode = grpMergeDF.explode(grp)
            grpMergeDF_explode = grpMergeDF_explode.loc[grpMergeDF_explode[grp]!="", ]
            final_abundanceDF = grpMergeDF_explode.groupby(grp, as_index=False).agg(abundance_keep)
            filename = f"{output_prefix}.{grp}.{exprType}.txt"
            final_abundanceDF.to_csv(filename, sep="\t", index=False)



if __name__ == '__main__':
    main()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------

