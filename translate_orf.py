#import needed modules
import sys
import re
from translate import translate_sequence
from find_orf import parse_sequence_from_path, find_first_orf, vet_nucleotide_sequence, vet_codon


def main():
    import argparse

    #Create a command-line parser object
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    #Define command-line arguments this script can receive for parser
    #arg1: the sequence
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    #arg2: path-for sequence in files
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    #arg2: start codon
    parser.add_argument('-s', '--start-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['AUG'],
            help = ('One or more possible start codons.'))
    #arg3: stop codon
    parser.add_argument('-x', '--stop-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['UAA', 'UAG', 'UGA'],
            help = ('One or more possible stop codons.'))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    #write needed code for main

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:            #allocate sequence
        sequence = args.sequence



    #find first orf
    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codons,
            stop_codons = args.stop_codons)

    #translate the orf found
    coded_protein = translate_sequence(orf)

    #print protein to STDOUT
    sys.stdout.write('{}\n'.format(coded_protein))




if __name__ == '__main__':
    main()
