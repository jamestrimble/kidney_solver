#!/usr/bin/env bash

"""convert_CMU_input_to_Glasgow.py: takes a single .input file in the CMU format
and converts it to separate .input (pairs) and .ndd (NDDs) format for Glasgow"""

import argparse
import os, sys
import logging
import csv

logging.basicConfig()
_log = logging.getLogger('convert')
_log.setLevel(logging.DEBUG)

def convert_and_write(input_file, output_base):
    """Reads CMU input_file and converts to two files, output_base.ginput and
    output_base.gndds containing only pairs and only NDDs, respectively"""

    # Do one full pass to figure out which IDs are NDDs, which are not,
    # and to create mappings of old IDs to new IDs
    pair_id_map = {}
    ndd_id_map = {}
    seen_vert_ids = set()
    dummy_edge_ct = 0
    real_edge_ct = 0
    real_edge_list = []
    _log.info("Figuring out who is an altruist and who is a pair ...")
    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = True
        for row in reader:

            # Files start with a distinguished row <num-verts> <num-edges>
            if header:
                exp_num_verts, exp_num_edges = int(row[0]), int(row[1])
                _log.info("Expecting {0} vertices (pairs+NDDs).".format(exp_num_verts))
                header=False
                continue

            # Files end with a distinguished row -1 -1 -1
            if int(row[0]) == -1:
                break

            # Meat of the file is rows of <src> <tgt> <weight> <is-dummy> <fail-prob>
            # If any vertex has a dummy edge incoming, it is an NDD
            src_id = int(row[0])
            tgt_id = int(row[1])
            weight = float(row[2])
            is_dummy = bool(int(row[3]))

            if is_dummy:
                # Track who is an NDD
                dummy_edge_ct += 1
                if tgt_id not in list(ndd_id_map.keys()):
                    new_ndd_id = len(ndd_id_map)
                    ndd_id_map[tgt_id] = new_ndd_id
            else:
                # Track non-dummy edges
                real_edge_ct += 1
                real_edge_list.append([src_id, tgt_id, weight])

            seen_vert_ids.add(src_id)
            seen_vert_ids.add(tgt_id)

    assert len(seen_vert_ids) == exp_num_verts
    assert dummy_edge_ct+real_edge_ct == exp_num_edges

    # Finish off the ID remapping for the pairs; NDDs were remapped above
    for vert_id in seen_vert_ids:
        if vert_id not in list(ndd_id_map.keys()):
            new_pair_id = len(pair_id_map)
            pair_id_map[vert_id] = new_pair_id

    assert len(ndd_id_map)+len(pair_id_map) == exp_num_verts
    _log.info("Found {0} altruists and {1} pairs; writing to new files ...".format(len(ndd_id_map), len(pair_id_map)))

    # How many edges are outgoing from the NDD set?
    ndd_outgoing_edge_ct = 0
    pairs_outgoing_edge_ct = 0
    for src, tgt, weight in real_edge_list:
        if src in list(ndd_id_map.keys()):
            ndd_outgoing_edge_ct += 1
        else:
            pairs_outgoing_edge_ct += 1

    assert ndd_outgoing_edge_ct+pairs_outgoing_edge_ct == real_edge_ct

    # Do a second pass through the CMU input file (now edges), this time 
    # writing to either the .gndd or the .ginput output file
    pair_outfile = "{0}.ginput".format(output_base)
    ndds_outfile = "{0}.gndds".format(output_base)
    pairs_outgoing_edges_written = 0
    ndds_outgoing_edges_written = 0
    with open(pair_outfile, 'w') as pf, open(ndds_outfile, 'w') as nf:
        pwriter = csv.writer(pf, delimiter='\t')  # writer for pairs file
        nwriter = csv.writer(nf, delimiter='\t')  # writer for ndds file
        
        # Headers for either file are <num-verts-of-that-type> <num-edges>
        pwriter.writerow([len(pair_id_map), pairs_outgoing_edge_ct])
        nwriter.writerow([len(ndd_id_map), ndd_outgoing_edge_ct])

        # Now write the edges, with translated vertex IDs, to the files
        for src, tgt, weight in real_edge_list:
            if src in list(ndd_id_map.keys()):
                assert tgt not in list(ndd_id_map.keys())  # can't have NDD->NDD
                nwriter.writerow([ndd_id_map[src], pair_id_map[tgt], weight])
                ndds_outgoing_edges_written+=1
            else:
                assert tgt not in list(ndd_id_map.keys())  # can't have pair->NDD
                pwriter.writerow([pair_id_map[src], pair_id_map[tgt], weight])
                pairs_outgoing_edges_written+=1
 
        # Each of the Glasgow files ends with -1 -1 -1
        EOF_arr = [-1, -1, -1]
        pwriter.writerow(EOF_arr)
        nwriter.writerow(EOF_arr)

    assert ndds_outgoing_edges_written == ndd_outgoing_edge_ct
    assert pairs_outgoing_edges_written == pairs_outgoing_edge_ct

    _log.info("Done writing to two files!")
        
def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert CMU input to Glasgow input.')
    parser.add_argument('--input-file', dest='input_file', required=True, 
                        help='full path to CMU .input file')
    parser.add_argument('--output-base', dest='output_base', required=True,
                        help='base path for Glasgow output files')
    args = parser.parse_args()
    
    input_file = args.input_file
    output_base = args.output_base
    _log.info("Taking input from {0} and outputting to {1}.(ginput,gndds)".format(input_file, output_base))
    
    # Sanity check arguments
    if not os.path.isfile(input_file):
        _log.error("Could not find input file {0}; quitting.".format(input_file))
        sys.exit(-1)

    # Perform the conversion and output to files
    convert_and_write(input_file, output_base)
    _log.info("Done!")


if __name__ == '__main__':
    main()
