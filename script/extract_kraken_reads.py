#!/usr/bin/env python
######################################################################
#extract_kraken_reads.py takes in a kraken-style output and kraken report
#and a taxonomy level to extract reads matching that level
#Copyright (C) 2019 Jennifer Lu, jlu26@jhmi.edu
#
#This file is part of KrakenTools
#KrakenTools is free software; oyu can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.
#
######################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 06/03/2019
#
#This program extracts reads classified by Kraken as a 
#specified taxonomy ID. Those reads are extracted into a new FASTA file.
#
#Required Parameters:
#   -k, --kraken, --kraken-file X.......kraken output file
#   -s, -s1, -1, -U X...................read file 
#                                       [FASTA/FASTQ - may be gzipped]
#   -s2, -2, X..........................second read file if paired 
#                                       [FASTA/FASTQ - may be gzipped]
#   -o, --output X......................output FASTA file with reads 
#   -t, --taxid, --taxids X.............list of taxonomy IDs to extract 
#                                       [separated by spaces]
#   -r, --report-file X.................kraken report file
#                                       [required only with --include-children/parents]
#Optional Parameters:
#   -h, --help..........................show help message.
#   --max X.............................only save the first X reads found
#   --include-children **...............include reads classified at lower levels 
#   --include-parents **................include reads classified at parent levels 
#                                       of taxids 
#   --append............................append extracted reads to output file if existing
#   --noappend..........................rewrite file if existing [default] 
#   --exclude...........................exclude the taxids specified
# ** by default, only reads classified exactly at taxids provided will be extracted
# ** if either of these are specified, a report file must also be provided 
######################################################################
import os, sys, argparse
import gzip
from time import gmtime
from time import strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#################################################################################
#Tree Class 
#usage: tree node used in constructing taxonomy tree  
#   includes only taxonomy levels and genomes identified in the Kraken report
class Tree(object):
    'Tree node.'
    def __init__(self, taxid, level_num, level_id, children=None, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)
#################################################################################
#process_kraken_output
#usage: parses single line from kraken output and returns taxonomy ID and readID
#input: kraken output file with readid and taxid in the
#   second and third tab-delimited columns
#returns: 
#   - taxonomy ID
#   - read ID
def process_kraken_output(kraken_line):
    l_vals = kraken_line.split('\t')
    if len(l_vals) < 5:
        return [-1, '']
    if "taxid" in l_vals[2]:
        temp = l_vals[2].split("taxid ")[-1]
        tax_id = temp[:-1]
    else:
        tax_id = l_vals[2]

    read_id = l_vals[1]
    if (tax_id == 'A'):
        tax_id = 81077
    else:
        tax_id = int(tax_id)
    return [tax_id, read_id]

#process_kraken_report
#usage: parses single line from report output and returns taxID, levelID
#input: kraken report file with the following tab delimited lines
#   - percent of total reads
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level
#       (U, - (root), - (cellular org), D, P, C, O, F, G, S)
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
#   - spaces + name
#returns:
#   - taxonomy ID 
#   - level number (number of spaces before name)
#   - level_type (type of taxonomy level - U, R, D, P, C, O, F, G, S, etc) 
def process_kraken_report(report_line):
    l_vals = report_line.strip().split('\t')
    try:
        int(l_vals[1])
    except ValueError:
        return []
    #Extract relevant information
    level_type = l_vals[3]
    taxid = int(l_vals[4])
    #Get spaces to determine level num
    spaces = 0
    for char in l_vals[-1]:
        if char == ' ':
            spaces += 1
        else:
            break
    level_num = int(spaces/2)
    return[taxid, level_num, level_type]
#################################################################################
#Main method 
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', dest='kraken_file', required=True,
        help='Kraken output file to parse')
    parser.add_argument('-s','-s1', '-1', '-U', dest='seq_file1', required=True,
        help='FASTA/FASTQ File containing the raw sequence letters.')
    parser.add_argument('-s2', '-2', dest='seq_file2', default= "",
        help='2nd FASTA/FASTQ File containing the raw sequence letters (paired).')
    parser.add_argument('-t', "--taxid",dest='taxid', required=True,
        nargs='+',
        help='Taxonomy ID[s] of reads to extract (space-delimited)')
    parser.add_argument('-o', "--output",dest='output_file', required=True,
        help='Output FASTA/Q file containing the reads and sample IDs')
    parser.add_argument('-o2',"--output2", dest='output_file2', required=False, default='',
        help='Output FASTA/Q file containig the second pair of reads [required for paired input]') 
    parser.add_argument('--append', dest='append', action='store_true',
        help='Append the sequences to the end of the output FASTA file specified.')
    parser.add_argument('--noappend', dest='append', action='store_false',
        help='Create a new FASTA file containing sample sequences and IDs \
              (rewrite if existing) [default].')
    parser.add_argument('--max', dest='max_reads', required=False, 
        default=100000000, type=int,
        help='Maximum number of reads to save [default: 100,000,000]')
    parser.add_argument('-r','--report',dest='report_file', required=False,
        default="",
        help='Kraken report file. [required only if --include-parents/children \
        is specified]')
    parser.add_argument('--include-parents',dest="parents", required=False, 
        action='store_true',default=False,
        help='Include reads classified at parent levels of the specified taxids')
    parser.add_argument('--include-children',dest='children', required=False,
        action='store_true',default=False,
        help='Include reads classified more specifically than the specified taxids')
    parser.add_argument('--exclude', dest='exclude', required=False,
        action='store_true',default=False,
        help='Instead of finding reads matching specified taxids, finds all reads NOT matching specified taxids') 
    parser.add_argument('--fastq-output', dest='fastq_out', required=False,
        action='store_true',default=False,
        help='Print output FASTQ reads [requires input FASTQ, default: output is FASTA]')
    parser.set_defaults(append=False)

    args=parser.parse_args()
    
    #Start Program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')
    
    #Check input 
    if (len(args.output_file2) == 0) and (len(args.seq_file2) > 0):
        sys.stderr.write("Must specify second output file -o2 for paired input\n")
        exit(1)

    #Initialize taxids
    save_taxids = {}
    for tid in args.taxid:
        save_taxids[int(tid)] = 0
    main_lvls = ['R','K','D','P','C','O','F','G','S']

    #STEP 0: READ IN REPORT FILE AND GET ALL TAXIDS 
    if args.parents or args.children:
        #check that report file exists
        if args.report_file == "": 
            sys.stderr.write(">> ERROR: --report not specified.")
            exit(1)
        sys.stdout.write(">> STEP 0: PARSING REPORT FILE %s\n" % args.report_file)
        #create tree and save nodes with taxids in the list 
        base_nodes = {} 
        r_file = open(args.report_file,'r')
        prev_node = -1
        for line in r_file:
            #extract values
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            [taxid, level_num, level_id] = report_vals
            if taxid == 0:
                continue 
            #tree root
            if taxid == 1:
                root_node = Tree(taxid, level_num, level_id)
                prev_node = root_node
                #save if needed
                if taxid in save_taxids:
                    base_nodes[taxid] = root_node
                continue
            #move to correct parent
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent 
            #determine correct level ID 
            if level_id == '-' or len(level_id) > 1:
                if prev_node.level_id in main_lvls:
                    level_id = prev_node.level_id + '1'
                else:
                    num = int(prev_node.level_id[-1]) + 1
                    level_id = prev_node.level_id[:-1] + str(num)
            #make node
            curr_node = Tree(taxid, level_num, level_id, None, prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node
            #save if taxid matches
            if taxid in save_taxids:
                base_nodes[taxid] = curr_node 
        r_file.close()
        #FOR SAVING PARENTS
        if args.parents:
            #For each node saved, traverse up the tree and save each taxid 
            for tid in base_nodes:
                curr_node = base_nodes[tid]
                while curr_node.parent != None:
                    curr_node = curr_node.parent
                    save_taxids[curr_node.taxid] = 0
        #FOR SAVING CHILDREN 
        if args.children:
            for tid in base_nodes:
                curr_nodes = base_nodes[tid].children
                while len(curr_nodes) > 0:
                    #For this node
                    curr_n = curr_nodes.pop()
                    if curr_n.taxid not in save_taxids:
                        save_taxids[curr_n.taxid] = 0
                    #Add all children
                    if curr_n.children != None:
                        for child in curr_n.children:
                            curr_nodes.append(child)
                    
    ##############################################################################
    sys.stdout.write("\t%i taxonomy IDs to parse\n" % len(save_taxids))
    sys.stdout.write(">> STEP 1: PARSING KRAKEN FILE FOR READIDS %s\n" % args.kraken_file)
    #Initialize values
    count_kraken = 0
    read_line = -1
    exclude_taxids = {} 
    if args.exclude:
        exclude_taxids = save_taxids 
        save_taxids = {} 
    #PROCESS KRAKEN FILE FOR CLASSIFIED READ IDS
    k_file = open(args.kraken_file, 'r')
    sys.stdout.write('\t0 reads processed')
    sys.stdout.flush()
    #Evaluate each sample in the kraken file
    save_readids = {}
    save_readids2 = {} 
    for line in k_file:
        count_kraken += 1
        if (count_kraken % 10000 == 0):
            sys.stdout.write('\r\t%0.2f million reads processed' % float(count_kraken/1000000.))
            sys.stdout.flush()
        #Parse line for results
        [tax_id, read_id] = process_kraken_output(line)
        if tax_id == -1:
            continue
        #Skip if reads are human/artificial/synthetic
        if (tax_id in save_taxids) and not args.exclude:
            save_taxids[tax_id] += 1
            save_readids2[read_id] = 0
            save_readids[read_id] = 0 
        elif (tax_id not in exclude_taxids) and args.exclude:
            if tax_id not in save_taxids:
                save_taxids[tax_id] = 1
            else:
                save_taxids[tax_id] += 1
            save_readids2[read_id] = 0
            save_readids[read_id] = 0 
        if len(save_readids) >= args.max_reads:
            break 
    #Update user
    k_file.close()
    sys.stdout.write('\r\t%0.2f million reads processed\n' % float(count_kraken/1000000.))
    sys.stdout.write('\t%i read IDs saved\n' % len(save_readids))
    ##############################################################################
    #Sequence files
    seq_file1 = args.seq_file1
    seq_file2 = args.seq_file2
    ####TEST IF INPUT IS FASTA OR FASTQ
    if(seq_file1[-3:] == '.gz'):
        s_file1 = gzip.open(seq_file1,'rt')
    else:
        s_file1 = open(seq_file1,'rt')
    first = s_file1.readline()
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        exit(1)
    s_file1.close()
    if filetype != 'fastq' and args.fastq_out:
        sys.stderr.write('ERROR: for FASTQ output, input file must be FASTQ\n')
        exit(1)
    ####ACTUALLY OPEN FILE
    if(seq_file1[-3:] == '.gz'):
        #Zipped Sequence Files
        s_file1 = gzip.open(seq_file1,'rt')
        if len(seq_file2) > 0:
            s_file2 = gzip.open(seq_file2,'rt')
    else:
        s_file1 = open(seq_file1, 'r')
        if len(seq_file2) > 0:
            s_file2 = open(seq_file2, 'r')
    #PROCESS INPUT FILE AND SAVE FASTA FILE
    sys.stdout.write(">> STEP 2: READING SEQUENCE FILES AND WRITING READS\n")
    sys.stdout.write('\t0 read IDs found (0 mill reads processed)')
    sys.stdout.flush()
    #Open output file
    if (args.append):
        o_file = open(args.output_file, 'a')
        if args.output_file2 != '':
            o_file2 = open(args.output_file2, 'a')
    else:
        o_file = open(args.output_file, 'w')
        if args.output_file2 != '':
            o_file2 = open(args.output_file2, 'w')
    #Process SEQUENCE 1 file 
    count_seqs = 0
    count_output = 0
    for record in SeqIO.parse(s_file1,filetype):
        count_seqs += 1
        #Print update
        if (count_seqs % 1000 == 0):
            sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
            sys.stdout.flush()
        #Check ID 
        test_id = str(record.id)
        if ("/1" in test_id) or ("/2" in test_id):
            test_id = test_id[:-2]
        #Sequence found
        if test_id in save_readids:
            count_output += 1
            #Print update
            sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
            sys.stdout.flush()
            #Save to file
            if args.fastq_out:
                SeqIO.write(record, o_file, "fastq")
            else:
                SeqIO.write(record, o_file, "fasta")
        #If no more reads to find 
        if len(save_readids) == count_output:
            break
    #Close files
    s_file1.close()
    o_file.close()
    sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)\n' % (count_output, float(count_seqs/1000000.)))
    sys.stdout.flush()
    count_output = 0
    count_seqs = 0
    if len(seq_file2) > 0:
        sys.stdout.write('\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
        sys.stdout.flush()
        for record in SeqIO.parse(s_file2, filetype):
            count_seqs += 1
            #Print update
            if (count_seqs % 1000 == 0):
                sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
                sys.stdout.flush()
            test_id = str(record.id)
            if ("/1" in test_id) or ("/2" in test_id):
                test_id = test_id[:-2]
            #Sequence found
            if test_id in save_readids:
                count_output += 1
                sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
                sys.stdout.flush()
                #Save to file
                if args.fastq_out:
                    SeqIO.write(record, o_file2, "fastq")
                else:
                    SeqIO.write(record, o_file2, "fasta")
            #If no more reads to find 
            if len(save_readids) == count_output:
                break
        s_file2.close()
        o_file2.close()
        #End Program
        sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)\n' % (count_output, float(count_seqs/1000000.)))
    
    #End Program
    sys.stdout.write('\t' + str(count_output) + ' reads printed to file\n')
    sys.stdout.write('\tGenerated file: %s\n' % args.output_file)
    if args.output_file2 != '':
        sys.stdout.write('\tGenerated file: %s\n' % args.output_file2)
    
    #End of program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    exit(0)

#################################################################################

if __name__ == "__main__":
    main()

#################################################################################
#################################END OF PROGRAM##################################
#################################################################################
