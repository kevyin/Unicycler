#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import random
import itertools
import multiprocessing
from .assembly_graph import AssemblyGraph
from .assembly_graph_copy_depth import determine_copy_depth
from .bridge_long_read_simple import create_simple_long_read_bridges
from .miniasm_assembly import make_miniasm_string_graph
from .bridge_miniasm import create_miniasm_bridges
from .bridge_long_read import create_long_read_bridges
from .bridge_spades_contig import create_spades_contig_bridges
from .bridge_loop_unroll import create_loop_unrolling_bridges
from .misc import int_to_str, float_to_str, quit_with_error, get_percentile, bold, \
    check_input_files, MyHelpFormatter, print_table, get_ascii_art, \
    get_default_thread_count, spades_path_and_version, makeblastdb_path_and_version, \
    tblastn_path_and_version, bowtie2_build_path_and_version, bowtie2_path_and_version, \
    samtools_path_and_version, java_path_and_version, pilon_path_and_version, \
    racon_path_and_version, bcftools_path_and_version, gfa_path, red
from .spades_func import get_best_spades_graph
from .blast_func import find_start_gene, CannotFindStart
from .read_ref import get_read_nickname_dict
from .pilon_func import polish_with_pilon_multiple_rounds, CannotPolish
from .vcf_func import make_vcf
from . import log
from . import settings
from .version import __version__


def main():
    """
    Script execution starts here.
    """

    # Fix the random seed so the program produces the same output every time it's run.
    random.seed(0)

    full_command = ' '.join(('"' + x + '"' if ' ' in x else x) for x in sys.argv)
    args = get_arguments()

    if args.input_assembly:
        graph = AssemblyGraph(args.input_assembly, args.overlap)
        trim_blunt_ends(graph, graph.overlap)

        output_filename = args.out
        graph.save_to_gfa(output_filename, save_copy_depth_info=False,
                          newline=True, include_insert_size=True)

def clean_up_spades_graph(graph):
    log.log_section_header('Cleaning graph')
    log.log_explanation('Unicycler now performs various cleaning procedures on the graph to '
                        'remove overlaps and simplify the graph structure. The end result is a '
                        'graph ready for bridging.', verbosity=1)
    graph.remove_all_overlaps()

    while True:
        graph.repair_multi_way_junctions()
        graph.remove_unnecessary_links()
        graph.expand_repeats()
        if not graph.remove_zero_length_segs():
            break
    while True:
        if not graph.merge_small_segments(5):
            break

    graph.normalise_read_depths()
    graph.renumber_segments()
    graph.sort_link_order()
    

def get_arguments():
    """
    Parse the command line arguments.
    """
    description = bold('graph modifications')
    this_script_dir = os.path.dirname(os.path.realpath(__file__))

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
    help_group.add_argument('--version', action='version', version='Unicycler v' + __version__,
                            help="Show Unicycler's version number")

    # input options
    input_group = parser.add_argument_group('Input')
    input_group.add_argument('--input_assembly', required=True,
                             help='gfa from an assembler (required)')
    input_group.add_argument('--overlap', required=False,
                             help='overlap size, default: None - guess from gfa link)', type=int, default=None)

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


def trim_blunt_ends(self, to_trim):
    """
    This function trims overlaps at blunt ends which would normally not be removed. It assumes that
    all overlaps in the graph are the same size.
    """
    if self.overlap == 0:
        log.log('Graph has no overlaps - overlap removal not needed')
        return

    # take note of all edges which start or end at a segment
    segment_edges = set()
    for start, ends in self.forward_links.items():
        segment_edges.add(start)
        for end in ends:
            segment_edges.add(-end)


    # Now we finally do the segment trimming!
    log.log('\nRemoving graph overlaps\n', 3)
    log.log('             Bases     Bases', 3)
    log.log('           trimmed   trimmed', 3)
    log.log(' Segment      from      from', 3)
    log.log('  number     start       end', 3)
    for seg_num, segment in self.segments.items():
        start_trim = 0
        end_trim = 0
        if seg_num not in segment_edges:
            end_trim = to_trim
        if -seg_num not in segment_edges:
            start_trim = to_trim

        segment.trim_from_start(start_trim)
        segment.trim_from_end(end_trim)
        log.log(str(seg_num).rjust(8) + str(start_trim).rjust(10) + str(end_trim).rjust(10), 3)

    log.log('Graph dead ends trimmed')
