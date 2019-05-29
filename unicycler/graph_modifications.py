#!/usr/bin/env python3

import argparse
import os
import sys
import random
import tempfile
import subprocess
from .assembly_graph import AssemblyGraph
from .misc import bold, MyHelpFormatter
from . import log


def main():
    """
    Script execution starts here.
    """

    # Fix the random seed so the program produces the same output every time it's run.
    random.seed(0)

    args = get_arguments()

    if args.input_assembly:
        graph = AssemblyGraph(args.input_assembly, args.overlap)
        trim_blunt_ends(graph, graph.overlap)

        output_filename = args.out
        graph.save_to_gfa(output_filename, save_copy_depth_info=False,
                          newline=True, include_insert_size=False)


def get_arguments():
    """
    Parse the command line arguments.
    """

    if '--helpall' in sys.argv or '--allhelp' in sys.argv or '--all_help' in sys.argv:
        sys.argv.append('--help_all')
    show_all_args = '--help_all' in sys.argv

    parser = argparse.ArgumentParser(description="modify gfa graphs", formatter_class=MyHelpFormatter,
                                     add_help=False)

    # Help options
    help_group = parser.add_argument_group('Help')
    help_group.add_argument('-h', '--help', action='help',
                            help='Show this help message and exit')
    help_group.add_argument('--help_all', action='help',
                            help='Show a help message with all program options')

    # input options
    input_group = parser.add_argument_group('Input')
    input_group.add_argument('--input_assembly', required=True,
                             help='gfa from an assembler (required)')
    input_group.add_argument('--overlap', required=False,
                             help='overlap size, default: None - guess from gfa link)', type=int, default=None)

    input_group.add_argument('--makeblastdb_path', type=str, default='makeblastdb',
                                help='Path to the makeblastdb executable'
                                if show_all_args else argparse.SUPPRESS)
    input_group.add_argument('--blastn_path', type=str, default='blastn',
                                help='Path to the blastn executable'
                                if show_all_args else argparse.SUPPRESS)

    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-o', '--out', required=True,
                              help='Output file (required)')
    output_group.add_argument('--verbosity', type=int, required=False, default=2,
                              help='R|Level of stdout and log file information (default: 1)\n  '
                                   '0 = no stdout, 1 = basic progress indicators, '
                                   '2 = extra info, 3 = debugging info')

    # If no arguments were used, print the entire help (argparse default is to just give an error
    # like '--out is required').
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Create an initial logger which doesn't have an output file.
    log.logger = log.Log(None, args.verbosity)

    return args


def make_blast_db(edges, segments, makeblastdb_path, take=None):

    fa_file = tempfile.mktemp(prefix="segments_", suffix=".fa")
    db_file = fa_file + '.db'

    with open(fa_file, 'w') as f:
        for seg_num in edges:
            seq = get_segment(seg_num, segments)
            if take:
                seq = seq[:take]
            f.write('>%s\n%s\n' % (seg_num,seq))

    command = [makeblastdb_path, '-dbtype', 'nucl', '-in', fa_file, '-out', db_file]
    log.log('  ' + ' '.join(command), 2)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = process.communicate()
    if err:
        log.log('\nmakeblastdb encountered an error:\n' + err.decode())
        raise err

    return db_file


class BlastHit(object):

    def __init__(self, blast_line):
        self.qseqid = ''
        self.pident, self.length, self.bitscore  = 0, 0, 0
        self.flip = False

        parts = blast_line.strip().split('\t')
        if len(parts) > 7:
            self.qseqid = parts[0]
            self.pident = float(parts[3])
            self.qstart = int(parts[6]) - 1
            self.bitscore = float(parts[7])
            self.length = float(parts[4])

    def __repr__(self):
        return 'BLAST hit: query=' + self.qseqid + ', ID=' + \
               str(self.pident) + ', length=' + str(self.length) + ', bitscore=' + \
               str(self.bitscore)


def search_blastn(blastdb, query_fa, blastn_path, num_alignments=1):

    command = [blastn_path, '-db', blastdb, '-query', query_fa, '-num_alignments', str(num_alignments), '-outfmt',
               '6 qseqid sstart send pident length qseq qstart bitscore', '-num_threads', '1']
    log.log('  ' + ' '.join(command), 2)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    blast_out, blast_err = process.communicate()
    process.wait()
    if blast_err:
        log.log('\nBLAST encountered an error:\n' + blast_err.decode())

    hits = []

    for line in blast_out.decode().splitlines():
        hit = BlastHit(line)
        hits.append(hit)
    return hits


def get_segment(seg_num, segments):
    return segments[seg_num].forward_sequence if seg_num > 0 else segments[-seg_num].reverse_sequence


def trim_blunt_ends(self, to_trim, makeblastdb_path='makeblastdb', blastn_path='blastn'):
    """
    This function trims overlaps at blunt ends which would normally not be removed. It assumes that
    all overlaps in the graph are the same size.
    """

    # take note of all edges which start or end at a segment
    segment_edges = set()
    for start, ends in self.forward_links.items():
        segment_edges.add(-start)
        for end in ends:
            segment_edges.add(end)

    # make the non-dead-end edge database
    db = make_blast_db(segment_edges, self.segments, makeblastdb_path, take=self.overlap)

    # prepare to search each dead end candidate against the db.
    query_fa = tempfile.mktemp(".fa", "query_de_candidates_")
    with open(query_fa, 'w') as f:
        for seg_num, segment in self.segments.items():
            # dead ends don't have any edges going in or out of them
            if seg_num not in segment_edges:
                seq = get_segment(seg_num, self.segments)
                f.write('>%s\n%s\n' % (seg_num, seq[:self.overlap]))
            if -seg_num not in segment_edges:
                seq = get_segment(-seg_num, self.segments)
                f.write('>%s\n%s\n' % (-seg_num, seq[:self.overlap]))

    hits = filter(lambda x: x.pident == 100, search_blastn(db, query_fa, blastn_path))

    log.log('\nRemoving graph overlaps\n', 3)
    log.log('             Bases     Bases', 3)
    log.log('           trimmed   trimmed', 3)
    log.log(' Segment      from      from', 3)
    log.log('  number     start       end', 3)
    for h in hits:
        # print(h)
        seg_num = int(h.qseqid)
        segment = self.segments[seg_num] if seg_num > 0 else self.segments[-seg_num]
        start_trim = int(h.length) if seg_num > 0 else 0
        end_trim = int(h.length) if seg_num < 0 else 0
        # trim
        if segment.get_length() - start_trim > self.overlap:
            segment.trim_from_start(start_trim)
        else:
            start_trim = 0
        if segment.get_length() - end_trim > self.overlap:
            segment.trim_from_end(end_trim)
        else:
            end_trim = 0
        log.log(str(seg_num).rjust(8) + str(start_trim).rjust(10) + str(end_trim).rjust(10), 3)

    log.log('Graph dead ends trimmed')
