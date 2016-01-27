__author__ = 'julianzaugg'

from sequence import *
from phylo import *
from annotation import *

from itertools import chain
import argparse
from collections import defaultdict

INPUT_SEQS = None
INPUT_SEQS_DICT = None
INPUT_TREE = None
OUT_LOCATION = None
SEQ_TRANSITIONS = [] # pairs of sequence names [["From Seq", "To Seq"]]
POSITIONS = []
SCALE = dict()

def _load_delimited(filename):
    with open(filename, 'r') as fh:
        data = [line.strip().split("\t") for line in fh.readlines()]
    return data

def _get_string_result(seq_pair, positions):
    out = []
    s1,s2 = seq_pair[0], seq_pair[1]
    t_name = s1.name + "_" + s2.name
    mean_value = 0.0
    cnt = 0
    for position in positions:
        s1_content = s1[position]
        s2_content = s2[position]
        # Don't add any score for gapped positions, ignore them and continue
        if s1_content == "-" or s2_content == "-": continue
        s1_position_score = SCALE[s1[position]]
        s2_position_score = SCALE[s2[position]]
        score_difference = s1_position_score - s2_position_score
        cnt += 1
        mean_value += score_difference
        out.append("%s\t%s\t%i\t%s\t%s\t%f\t%f\t%f\t%s" % (s1.name, # Name of first sequence
                                                     s2.name, # Name of second sequence
                                                     position+1, # position (+1 to adjust back to original number)
                                                     s1_content, # AA content for first sequence
                                                     s2_content, # AA content for second sequence
                                                     s1_position_score, # score for first sequence at the position
                                                     s2_position_score, # score for second sequence at the position
                                                     score_difference, # difference in scores
                                                     t_name)) # Name given to the transition
    return out, mean_value/cnt

def _dir_check(directory):
    if not os.path.exists(directory):
        raise StandardError("The directory %s does not exist or the location was incorrectly specified" % directory)

def _process_arguments(my_parser, my_args):
    global INPUT_SEQS, INPUT_TREE, SEQ_TRANSITIONS,INPUT_SEQS_DICT, POSITIONS, SCALE
    OUT_LOCATION = my_args.output_location
    _dir_check(OUT_LOCATION)

    # Load the sequences
    INPUT_SEQS = read_fasta_file(my_args.input_seqs, Protein_Alphabet)
    INPUT_SEQS_DICT = dict([(s.name, s)for s in INPUT_SEQS])
    
    # Load the tree file
    if my_args.tree_file:
        INPUT_TREE = load_tree(my_args.tree_file)
    
    # Load transition pairs, grab and store relevant sequences.     
    if my_args.seq_names_file:
        data = _load_delimited(my_args.seq_names_file)
        for pair in data:
            pair_seqs = []
            for name in pair:
                if name not in INPUT_SEQS_DICT:
                    raise StandardError("The sequence name %s defined in the transition file does not match any of the input sequences" % name)
                pair_seqs.append(INPUT_SEQS_DICT[name])
            SEQ_TRANSITIONS.append(pair_seqs)
    if not my_args.tree_file and not my_args.seq_names_file:
        print "No Tree or transtion file provided, will compare all sequences in input file"
        for seq in INPUT_SEQS:
            for seq2 in INPUT_SEQS:
                if seq != seq2:
                    SEQ_TRANSITIONS.append([seq, seq2])

    # Load Positions if provided
    if my_args.positions_file:
        with open(my_args.positions_file, 'r') as fh:
            POSITIONS = map(lambda x: x - 1, map(int, [line.strip().split() for line in fh.readlines()][0]))
    else: 
        POSITIONS = xrange(len(INPUT_SEQS[0]))
    
    # Load the scale
    with open(my_args.scale_file, 'r') as fh:
        data = [x.split("\t") for line in fh.readlines() for x in line.strip().split("\r")]
        data = [[x,float(y)] for x,y in data]
        SCALE = dict(data)

    # If a tree and transition names have been provided, we will just use the transition name file,
    # scoring just the transitions between specified sequences
    # Otherwise we use the hierarchy of the tree.
    all_result_strings = []
    all_mean_scores = defaultdict(defaultdict)
    if INPUT_TREE and SEQ_TRANSITIONS or SEQ_TRANSITIONS:
        for seq_pair in SEQ_TRANSITIONS:
            print "Scoring %s - %s transition" % (seq_pair[0].name, seq_pair[1].name) 
            result_string_list, mean_score = _get_string_result(seq_pair, POSITIONS)
            all_result_strings += result_string_list
            all_mean_scores[seq_pair[0].name][seq_pair[1].name] = mean_score
    else:
        for node in INPUT_TREE:
            node_seq = INPUT_SEQS_DICT[node.label]
            for child in node.children:
                child_seq = INPUT_SEQS_DICT[child.label]
                print "Scoring %s - %s transition" % (node_seq.name, child_seq.name)
                result_string_list, mean_score = _get_string_result((node_seq, child_seq), POSITIONS)
                all_result_strings += result_string_list
                all_mean_scores[node_seq.name][child_seq.name] = mean_score
    
    with open(my_args.output_location + "position_results.txt", 'w') as fh:
        header =  "Seq1_name\tSeq2_name\tPosition\tSeq1_content\tSeq2_content\tSeq1_score\tSeq2_score\tScore_difference\tTransitionName"
        print >> fh, header
        for r_string in all_result_strings:
            print >> fh, r_string
            
    with open(my_args.output_location + "mean_results.txt", 'w') as fh:
        header = "Seq1_name\tSeq2_name\tMean_score"
        print >> fh, header
        for k,v in all_mean_scores.iteritems():
            for k2,v2 in v.iteritems():
                print >> fh,"%s\t%s\t%f" % (k, k2, v2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Description to be added')

    parser.add_argument('-i', '--input_seqs', help='Input FASTA file', required=True)
    parser.add_argument('-s', '--scale_file', help='Input scale file, should be tab delimited with first column\
            containing the amino acid and the second column the score', required=True)
    parser.add_argument('-o', '--output_location', help='Output location, where to save the results', required=True)

    parser.add_argument('-t', '--tree_file', help='Input tree file, must be in newick format', required=False)
    parser.add_argument('-p', '--positions_file', help='Input positions file, if provided, should be a single column\
            containing position numbers. Assumes indexing starts at 1. If not provided, all positions will be scanned', required=False)
    parser.add_argument('-seq_names', '--seq_names_file', help='Input sequence transitions of interest', required=False)

    args = parser.parse_args()
    # args = parser.parse_args(["-i","Data/seqs.txt" ,"--scale_file", "Data/scale.txt", "-o" ,"Data/", "-seq_names", "Data/temp.txt", "-p", "Data/positions.txt", "--tree_file", "Data/tree.txt"])
    # args = parser.parse_args(["-i","Data/seqs.txt" ,"--scale_file", "Data/scale.txt", "-o" ,"Data/", "-seq_names", "Data/temp.txt", "-p", "Data/positions.txt"])
    # args = parser.parse_args(["-i","Data/seqs.txt" ,"--scale_file", "Data/scale.txt", "-o" ,"Data/", "-p", "Data/positions.txt","--tree_file", "Data/tree.txt"])

    # args = parser.parse_args(["-i", "Data/CYP3_joint_reconstruction.txt", "--output_location", "Data/", "--scale_file", "Data/TP_all_residues.txt", "--seq_names_file", "Data/Transition.txt"])
    _process_arguments(parser, args)

