#! /usr/bin/env python3
#the shebang line above


#Objective: The goal is for you to practice writing your own scripts in Python "from scratch"
#This script only works for STRINGS and sequences in FILES.

#import needed modules
import sys
import re


#function: vet sequence
def vet_nucleotide_sequence(sequence):
    """
    Return None if `sequence` is a valid RNA or DNA sequence, else raise exception. 

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)

    Returns
    -------
    None
        Return nothing (None) if sequence is valid, otherwise raise an
        exception.

    Examples
    --------
    >>> vet_nucleotide_sequence('ACGTACGT') == None
    True

    >>> vet_nucleotide_sequence('not a valid sequence')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: 'not a valid sequence'

    Don't allow mixing of DNA and RNA!
    >>> vet_nucleotide_sequence('AUTGC') == None
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: 'AUTGC'

    Don't allow whitespace (or other characters) before, within, or after!
    >>> vet_nucleotide_sequence(' ACGT ACGT ') == None
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: ' ACGT ACGT '
    """
    ##########################################################################
    # `rna_pattern_str` and `dna_pattern_str` are be regular expressions
    # that will match any string of RNA and DNA bases, respectively (and only
    # strings of RNA and DNA bases).
    # Read the docstring above for additional clues.
    rna_pattern_str = r'^[AUCGaucg]*$'
    dna_pattern_str = r'^[ATCGatcg]*$'
    ##########################################################################

    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    if rna_pattern.match(sequence):
        return
    if dna_pattern.match(sequence):
        return
    else:
        raise Exception("Invalid sequence: {0!r}".format(sequence))


#function: vet codon in sequence
def vet_codon(codon):
    """
    Return None if `codon` is a valid RNA codon, else raise an exception. 

    Parameters
    ----------
    codon : str
        A string representing a codon (upper or lower-case)

    Returns
    -------
    None
        Return nothing (None) if codon is valid, otherwise raise an
        exception.

    Examples
    --------
    >>> vet_codon('AUG') == None
    True

    >>> vet_codon('aug') == None
    True

    >>> vet_codon('ATG')
    Traceback (most recent call last):
        ...
    Exception: Invalid codon: 'ATG'

    >>> vet_codon('AUGG')
    Traceback (most recent call last):
        ...
    Exception: Invalid codon: 'AUGG'
    """
    ##########################################################################
    # `codon_pattern_str` is a regular expression that will match any
    # codon (but only one codon). Only valid codons.
    # Read the docstring above for additional clues.
    codon_pattern_str = r'^[AUGCaugc]{3}$'
    ##########################################################################

    codon_pattern = re.compile(codon_pattern_str)

    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))


#function: find first orf
def find_first_orf(sequence,
        start_codons = ['AUG'],
        stop_codons = ['UAA', 'UAG', 'UGA']):
    """
    Return the first open-reading frame in the DNA or RNA `sequence`.

    An open-reading frame (ORF) is the part of an RNA sequence that is
    translated into a peptide. It must begin with a start codon, followed by
    zero or more codons (triplets of nucleotides), and end with a stop codon.
    If there are no ORFs in the sequence, an empty string is returned.

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)
    start_codons : list of strings
        All possible start codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.
    stop_codons : list of strings
        All possible stop codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.

    Returns
    -------
    str
        An uppercase string of the first ORF found in the `sequence` that
        starts with any one of the `start_codons` and ends with any one of the
        `stop_codons`. If no ORF is found an empty string is returned.

    Examples
    --------
    When the whole RNA sequence is an ORF:
    >>> find_first_orf('AUGGUAUAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When the whole DNA sequence is an ORF:
    >>> find_first_orf('ATGGTATAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When there is no ORF:
    >>> find_first_orf('CUGGUAUAA', ['AUG'], ['UAA'])
    ''

    When there is are bases before and after ORF:
    >>> find_first_orf('CCAUGGUAUAACC', ['AUG'], ['UAA'])
    'AUGGUAUAA'
    """
    
    # Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    # Make sure the codons are valid
    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    # Get copies of everything in uppercase
    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    # Make sure seq is RNA
    seq = seq.replace('T', 'U')


    ##########################################################################
    # `orf_pattern_str` is a regular expression that will match an
    # open reading frame within a string of RNA bases. 
    # Read the docstring above for additional clues.
    #I create variables as "strings" to take in start_codons & stop_codons "lists" so it can use them in the "regex"
    
    strt = '|'.join(start_codons) #takes the start codon provided and saves it as a string named 'strt'
    stp = '|'.join(stop_codons) #takes the stop codon provided and takes it as a string named 'stp'

    #here the `regex`
    orf_pattern_str = r'(' + strt +r')([AUCG]{3})*('+ stp +r')'
        #breakdown
            #the regex is in 3 groups
            #first group is the start codon
            #   literarily it is 'AUG', but it is dynamic, and (' + strt +r') fixes that
            #the second group is ([AUGC]{3})* which captures zero or more sets of codon after the start codon
            #the last group is ((UAA)|(UAG)|(UGA)) dynamically captured as ('+ stp +r')
    ##########################################################################

    # Create the regular expression object
    orf_pattern = re.compile(orf_pattern_str)
    # Search the sequence
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''



#function: file read and save file content as sequence
def parse_sequence_from_path(path):
    # Try to open the path to read from it, and handle exceptions if they
    # arrise
    try:
        file_stream = open(path, 'r')
    except FileNotFoundError as e:
        sys.stderr.write("Sorry, couldn't find path {}".format(path))
        raise e
    except IsADirectoryError as e:
        sys.stderr.write("Sorry, path {} appears to be a directory".format(
                path))
        raise e
    except:
        sys.stderr.write("Sorry, something went wrong when trying to open {}".format(
                path))
        raise
    # If we've reached here, the file is open and ready to read
    sequence = ''
    # A for loop to visit each line in the file
    for line in file_stream:
        # Strip whitespace from the line and concatenate it to the end of the
        # sequence
        sequence += line.strip()
    return sequence






#function: translate first orf found
def translate_sequence(sequence, genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}, start_pos = 0):
 

    """Translates the first ORF sequence of RNA into a sequence of amino acids.
    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.
    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).
    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').
    Returns
    -------
    str
        A string of the translated amino acids.
    """
    #find first orf
    first_orf_seq = find_first_orf(sequence)

    # ensure copies of first_orf is uppercase
    first_orf_sequence = first_orf_seq.upper()

    #translate the first orf sequence
    protein = ""
    for i in range(0, len(first_orf_sequence) - (len(first_orf_sequence) % 3), 3):
        codon = first_orf_sequence[i:i + 3]
        if genetic_code[codon] == "*":
            break
        protein += genetic_code[codon]
    return protein



#print(translate_sequence('CCAUGGUAUAACC'))




#Write a main function
#for command-line interface!

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
    #arg3: start codon
    parser.add_argument('-s', '--start-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['AUG'],
            help = ('One or more possible start codons.'))
    #arg4: stop codon
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
