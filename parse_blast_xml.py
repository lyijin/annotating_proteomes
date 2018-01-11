#!/usr/bin/env python3

"""
> filter_blast_xml.py <

Python script parses through BLAST's XML output, and filters out hits that
are below a cutoff e-value (thus removing the need to re-run searches with
varying e-values).

The output of the script is XML (--xml), table & each hit on one line (--table),
table & each hsp on each table (--vtable).
"""

import argparse
import ast
from collections import OrderedDict
import itertools
import sys
import xml.etree.ElementTree

parser = argparse.ArgumentParser(description="""
Python script parses through BLAST's XML output, and filters out hits that
are below a cutoff e-value (thus removing the need to re-run searches with
varying e-values).

The output of the script is XML (--xml), table & each hit on one line (--table),
table & each hsp on each table (--vtable).""")

parser.add_argument('blast_xml', metavar="xml_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="BLAST XML filename.")
parser.add_argument('-e', '--e_value', metavar='e_value', type=float,
                    help='set maximum E value.')
parser.add_argument('-b', '--bit_score', metavar='bit_score', type=float,
                    help='set minimum total bit score value.')
parser.add_argument('-c', '--coverage', metavar='coverage_pct', type=float,
                    help='set minimum coverage (in %%) for individual hits.')
parser.add_argument('-i', '--identity', metavar='identity_pct', type=float,
                    help='set minimum local identity (in %%).')
parser.add_argument('-t', '--top', metavar='n', type=int,
                    help='select for top n hits.')
parser.add_argument('-C', '--compress', action='store_true',
                    help="""print compressed output, showing all hits for each
                    query on a single line. Can be _very_ verbose.""")
parser.add_argument('-N', '--remove_N', metavar='query_file',
                    type=argparse.FileType('r'),
                    help="""discards Ns before calculating coverage %% (best 
                    for scaffolds). Requires FASTA file of query sequences.""")
fasta_opt = parser.add_mutually_exclusive_group(required=True)
fasta_opt.add_argument('--xml', action='store_const', dest='output_format',
                       const='xml', help='produce XML output.')
fasta_opt.add_argument('--table', action='store_const', dest='output_format',
                       const='table', help='produce tabular output.')
fasta_opt.add_argument('--vtable', action='store_const', dest='output_format',
                       const='vtable',
                       help='produce verbose tabular output (1 line per Hsp).')
parser.add_argument('--noheader', action='store_true',
                    help='disable printing of header.')
args = parser.parse_args()

# abbreviated structure BLAST XML, correct as of 2.2.28+:
# <?xml version="1.0"?>
#<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
#<BlastOutput>
#	<BlastOutput_iterations>
#		<Iteration>                                     <!-- block recurs -->
#			<Iteration_iter-num></Iteration_iter-num>
#			<Iteration_query-ID></Iteration_query-ID>
#			<Iteration_query-def></Iteration_query-def>
#			<Iteration_query-len></Iteration_query-len>
#			<Iteration_hits> 
#				<Hit>                                   <!-- block recurs -->
#					<Hit_num></Hit_num>
#					<Hit_id></Hit_id>
#					<Hit_def></Hit_def>
#					<Hit_accession></Hit_accession>
#					<Hit_len></Hit_len>
#					<Hit_hsps>
#						<Hsp>                           <!-- block recurs -->
#							<Hsp_num></Hsp_num>
#							<Hsp_bit-score></Hsp_bit-score>
#							<Hsp_score></Hsp_score>
#							<Hsp_evalue></Hsp_evalue>
#							<Hsp_query-from></Hsp_query-from>
#							<Hsp_query-to></Hsp_query-to>
#							<Hsp_hit-from></Hsp_hit-from>
#							<Hsp_hit-to></Hsp_hit-to>
#							<Hsp_query-frame></Hsp_query-frame>
#							<Hsp_hit-frame></Hsp_hit-frame>
#							<Hsp_identity></Hsp_identity>
#							<Hsp_positive></Hsp_positive>
#							<Hsp_gaps></Hsp_gaps>
#							<Hsp_align-len></Hsp_align-len>
#							<Hsp_qseq></Hsp_qseq>
#							<Hsp_hseq></Hsp_hseq>
#							<Hsp_midline></Hsp_midline>
#						</Hsp>
#					</Hit_hsps>
#				</Hit>
#			</Iteration_hits>
#			<Iteration_stat>
#				<Statistics>
#					<Statistics_db-num></Statistics_db-num>
#					<Statistics_db-len></Statistics_db-len>
#					<Statistics_hsp-len></Statistics_hsp-len>
#					<Statistics_eff-space></Statistics_eff-space>
#					<Statistics_kappa></Statistics_kappa>
#					<Statistics_lambda></Statistics_lambda>
#					<Statistics_entropy></Statistics_entropy>
#				</Statistics>
#			</Iteration_stat>
#		</Iteration>
#	</BlastOutput_iterations>
#</BlastOutput>

def consolidate_coords(zipped_coord_tuples):
    sorted_coords = sorted(list(zipped_coord_tuples))
    consolidated_coords = []
            
    min_start = sorted_coords[0][0]
    max_end = sorted_coords[0][1]
    for a in sorted_coords:
        # max_end + 1 consolidates coords e.g. (1, 3), (4, 6) --> (1, 6)
        if min_start <= a[0] <= max_end + 1:
            max_end = max(max_end, a[1])
        else:
            consolidated_coords.append((min_start, max_end))
            min_start = a[0]
            max_end = a[1]

    consolidated_coords.append((min_start, max_end))
    
    return consolidated_coords

def calc_coords_length(consolidated_coords):
    total_length = 0
    
    for c in consolidated_coords:
        total_length += abs(c[0] - c[1]) + 1
    
    return total_length

def compress_output(output):
    c_output = ''
    if not args.noheader:
        c_output += '\t'.join(['Query', 'Hit accession', 'Hit description',
                               'Query length', 'Hit length', 
                               'Query (start, end)', 'Hit (start, end)',
                               'Frame', 'Max bit score', 'Total bit score',
                               'Identity', 'Identity %',
                               'Coverage %', 'Expect']) + '\n'
    
    # create a dictionary to store lines in the output that belong to each query
    lines_per_query = OrderedDict()
    
    for line in output.split('\n')[1:]:
        if not line: break
        
        cols = line.split('\t')
        if cols[0] not in lines_per_query:
            lines_per_query[cols[0]] = []
        
        lines_per_query[cols[0]].append(line)
    
    for q in lines_per_query:
        # consolidate all the BLAST hits into one line
        hit_accessions = ', '.join([x.split('\t')[1] for x in lines_per_query[q]])
        hit_descriptions = ', '.join([x.split('\t')[2] for x in lines_per_query[q]])
        query_length = lines_per_query[q][0].split('\t')[3]
        hit_lengths = ', '.join([x.split('\t')[4] for x in lines_per_query[q]])
        query_coords = ', '.join([x.split('\t')[5] for x in lines_per_query[q]])
        query_coords = list(ast.literal_eval('[{}]'.format(query_coords)))
        query_coords = consolidate_coords(query_coords)
        hit_coords = '; '.join(['{}: {}'.format(n + 1, x.split('\t')[6]) for n, x in enumerate(lines_per_query[q])])
        frames = '; '.join(['{}: {}'.format(n + 1, x.split('\t')[7]) for n, x in enumerate(lines_per_query[q])])
        max_bit_scores = max([float(x.split('\t')[8]) for x in lines_per_query[q]])
        total_bit_scores = sum([float(x.split('\t')[9]) for x in lines_per_query[q]])
        id_numerator = sum([int(x.split('\t')[10].split('/')[0]) for x in lines_per_query[q]])
        id_divisor = sum([int(x.split('\t')[10].split('/')[1]) for x in lines_per_query[q]])
        identities = '{}/{}'.format(id_numerator, id_divisor)
        id_pct = '{:.2f}%'.format(id_numerator / id_divisor * 100)
        coverage_pct = '{:.2f}%'.format(calc_coords_length(query_coords) / int(query_length) * 100)
        expects = '{:.2e}'.format(min([float(x.split('\t')[13]) for x in lines_per_query[q]]))
        
        o = [q, hit_accessions, hit_descriptions, query_length, hit_lengths,
             query_coords, hit_coords, frames, max_bit_scores, total_bit_scores,
             identities, id_pct, coverage_pct, expects]
        
        c_output += '\t'.join([str(x) for x in o]) + '\n'
    
    return c_output
    
if args.remove_N:
    import parse_fasta
    query_fasta_seq = parse_fasta.get_all_sequences(args.remove_N, 'fasta')

tree = xml.etree.ElementTree.parse(args.blast_xml)
root = tree.getroot()

# get list of <Iteration></Iteration>
blastoutput_iterations = root.find('BlastOutput_iterations')
iterations = blastoutput_iterations.findall('Iteration')

# remove iterations that contain
#   "<Iteration_message>No hits found</Iteration_message>"
for i in iterations:
    if i.find('Iteration_message') is not None:
        if i.find('Iteration_message').text == 'No hits found':
            blastoutput_iterations.remove(i)

# iterate through remaining iterations, and apply specified filters:
#
# in case there are multiple stretches of sequences in the query
# that matches the subject...
# local extreme metrics:
#   e_value: defined to be the minimum e-value
#   max_score: stretch that scores best
#   max_identity: stretch that has the highest % identity
# total metrics:
#   total_score: all scores added together
#   query_coverage: sum(lengths of query with matching pieces) / query_length

for i in blastoutput_iterations.findall('Iteration'):
    iteration_hits = i.find('Iteration_hits')
    hits = iteration_hits.findall('Hit')
    for h in hits:
        hsps = h.find('Hit_hsps').findall('Hsp')
        bit_scores = [float(x.find('Hsp_bit-score').text) for x in hsps]
        e_values = [float(x.find('Hsp_evalue').text) for x in hsps]
        identities = [int(x.find('Hsp_identity').text) for x in hsps]
        align_lens = [int(x.find('Hsp_align-len').text) for x in hsps]
        identity_pcts = [x / y * 100 for x, y in zip(identities, align_lens)]
        query_froms = [int(x.find('Hsp_query-from').text) for x in hsps]
        query_tos = [int(x.find('Hsp_query-to').text) for x in hsps]
        
        # to avoid coverage_pct > 100% if a certain section of the query has
        # multiple matches to the hit.
        c = consolidate_coords(zip(query_froms, query_tos))
        # change divisor if --remove_N is called
        if args.remove_N:
            query_length = len(query_fasta_seq[i.find('Iteration_query-def').text].upper().replace('N', ''))
        else:
            query_length = int(i.find('Iteration_query-len').text)
        
        coverage_pct = sum([abs(x - y) + 1 for x, y in c]) \
                       / query_length * 100
        
        # filter based on specified cutoffs
        remove_flag = False
        if args.e_value is not None:
            if args.e_value < min(e_values):
                remove_flag = True
        
        if args.bit_score is not None:
            if args.bit_score > sum(bit_scores):
                remove_flag = True
        
        if args.coverage is not None:
            #debug: print (i, h, coverage_pct)
            if args.coverage > coverage_pct:
                remove_flag = True
        
        if args.identity is not None:
            if args.identity > max(identity_pcts):
                remove_flag = True
        
        if remove_flag:
            iteration_hits.remove(h)
        
    # remove iteration if all its hits was removed by the cutoffs
    if iteration_hits.find('Hit') is None:
        blastoutput_iterations.remove(i)

if args.output_format == 'xml':
    # print to standard output - use "unicode" to make sure output is writable
    # to screen!
    tree.write(sys.stdout, encoding='unicode')
elif 'table' in args.output_format:
    output = ''
    if args.output_format == 'vtable' and not args.noheader:
        output += '\t'.join(['Query', 'Hit accession', 'Hit description',
                             'Query length', 'Hit length', 
                             'Query start', 'Query end', 'Hit start', 'Hit end',
                             'Frame', 'Bit score', 'Identity', 'Identity %',
                             'Coverage %', 'Expect']) + '\n'
    elif args.output_format == 'table' and not args.noheader:
        output += '\t'.join(['Query', 'Hit accession', 'Hit description',
                             'Query length', 'Hit length', 
                             'Query (start, end)', 'Hit (start, end)',
                             'Frame', 'Max bit score', 'Total bit score', 
                             'Identity', 'Identity %', 
                             'Coverage %', 'Expect']) + '\n'
    
    for i in blastoutput_iterations.findall('Iteration'):
        query_def = i.find('Iteration_query-def').text
        
        hit_count = 0
        iteration_hits = i.find('Iteration_hits')
        hits = iteration_hits.findall('Hit')
        for h in hits:
            hit_count += 1
            if args.top:
                if hit_count > args.top:
                    break
            
            hit_ID = h.find('Hit_id').text
            hit_def = h.find('Hit_def').text
            hit_length = h.find('Hit_len').text
            
            hsps = h.find('Hit_hsps').findall('Hsp')
            query_froms = [int(x.find('Hsp_query-from').text) for x in hsps]
            query_tos = [int(x.find('Hsp_query-to').text) for x in hsps]
            hit_froms = [int(x.find('Hsp_hit-from').text) for x in hsps]
            hit_tos = [int(x.find('Hsp_hit-to').text) for x in hsps]
            query_frames = [int(x.find('Hsp_query-frame').text) for x in hsps]
            bit_scores = [float(x.find('Hsp_bit-score').text) for x in hsps]
            identities = [int(x.find('Hsp_identity').text) for x in hsps]
            align_lens = [int(x.find('Hsp_align-len').text) for x in hsps]
            identity_pcts = [x / y * 100 for x, y in
                             zip(identities, align_lens)]
            e_values = [float(x.find('Hsp_evalue').text) for x in hsps]
            
            # change divisor if --remove_N is called
            if args.remove_N:
                query_length = len(query_fasta_seq[i.find('Iteration_query-def').text].upper().replace('N', ''))
            else:
                query_length = int(i.find('Iteration_query-len').text)
            
            # again, to avoid coverage > 100%
            coverage_pcts = [
                (abs(x - y) + 1) / query_length * 100
                for x, y in zip(query_froms, query_tos)]
            
            c = consolidate_coords(zip(query_froms, query_tos))
            consolidated_coverage_pcts = [
                (abs(x - y) + 1) / query_length * 100
                for x, y in c]
            
            if args.output_format == 'vtable':
                for n, hsp in enumerate(hsps):
                    o = [query_def, hit_ID, hit_def, query_length, hit_length,
                         query_froms[n], query_tos[n], hit_froms[n], 
                         hit_tos[n], query_frames[n],
                         '{:.2f}'.format(bit_scores[n]),
                         '{}/{}'.format(identities[n], align_lens[n]),
                         '{:.2f}%'.format(identity_pcts[n]),
                         '{:.2f}%'.format(coverage_pcts[n]),
                         '{:.2e}'.format(e_values[n])]
                    output += '\t'.join([str(x) for x in o]) + '\n'
            else:
                o = [query_def, hit_ID, hit_def, query_length, hit_length,
                     ', '.join(['({},{})'.format(x, y) for x, y in
                                sorted(zip(query_froms, query_tos))]),
                     ', '.join(['({},{})'.format(x, y) for x, y in
                                sorted(zip(hit_froms, hit_tos))]),
                     ', '.join([str(x) for x in query_frames]),
                     '{:.2f}'.format(max(bit_scores)),
                     '{:.2f}'.format(sum(bit_scores)),
                     '{}/{}'.format(sum(identities), sum(align_lens)),
                     '{:.2f}%'.format(sum(identities) / sum(align_lens) * 100),
                     '{:.2f}%'.format(sum(consolidated_coverage_pcts)),
                     '{:.2e}'.format(e_values[0])]
                output += '\t'.join([str(x) for x in o]) + '\n'
    
    # strip is required to remove trailing \n
    if args.compress:
        print (compress_output(output).strip())
    else:
        print (output.strip())