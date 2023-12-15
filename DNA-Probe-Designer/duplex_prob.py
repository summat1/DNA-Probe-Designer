import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from Bio.SeqUtils import GC

###################################################################################################

def calc_duplex_prob(sam_filename, bed_filename, temp):
    '''
    Calculates probability of all probes in a probeset forming a duplex with target sequence at a given temp.
        Arguments:
            - sam_filename [str] : relative path to .sam file containing alignment scores for probe candidates
            - bed_filename [str] : relative path to .bed file containing sequences of final probeset
            - temp [int] : temp at which to predict duplex probability
        Outputs:
            - probs [np.ndarray] : duplex probabilities of probes in probeset
    '''
    # read in final probeset BED file #
    with open(bed_filename) as file:
        probeset = [line.split('\t')[3] for line in file]

    # read and check SAM file format #
    with open(sam_filename) as file:
        sam = []
        for line in file:
            parts = line.split('\t')
            # should be 19 columns
            if len(parts) < 19:
                raise ValueError("SAM file format is incorrect.")
            # probe sequence should only have GATC
            if not set(parts[9]).issubset(set("GATC")):
                raise ValueError("Probe sequence must contain only G, A, T, C in SAM file.")
            if ':' not in parts[12]:
                raise ValueError("Alignment score in SAM file does not contain expected characters.")
            align_parts = parts[12].split(':')
            if len(align_parts) < 3 or not align_parts[2].isnumeric():
                raise ValueError("Alignment score format is incorrect in SAM file.") 
            sam.append(line)

    
    # set up LDA model to predict duplex probability based on probe length, GC content, and alignment score to target seq #
    # Values from models published in Beliveau, et al. (2018) #
    temps = np.array([32, 37, 42, 47, 52, 57])
    if temp not in temps:
        raise ValueError(f"Invalid temperature value: {temp}. Valid values are {temps}")
    coefs = np.array([[-0.14494789, 0.18791679, 0.02588474],
                    [-0.13364364, 0.22510179, 0.05494031],
                    [-0.09006122, 0.25660706, 0.1078303],
                    [-0.01593182, 0.24498485, 0.15753649],
                    [0.01860365, 0.1750174, 0.17003374],
                    [0.03236755, 0.11624593, 0.24306498]])
    intercepts = np.array([-1.17545204, -5.40436344, -12.45549846,
                        -19.32670233, -20.11992898, -23.98652919])
    classes = np.array([-1, 1])

    # initialize LDA model with temp-specific parameters #
    index = np.where(temps==temp)
    clf = LinearDiscriminantAnalysis()
    clf.coef_ = coefs[index]
    clf.intercept_ = intercepts[index]
    clf.classes_ = classes

    # collate inputs for LDA model as [probe length, alignment score, GC content] #
    clf_inputs = []
    for k in range(len(sam)):
        probe_seq = sam[k].split('\t')[9]
        align_score = sam[k].split('\t')[12].split(':')[2]
        GC_content = GC(probe_seq)

        if probe_seq in probeset:
            clf_inputs.append([len(probe_seq), int(align_score), GC_content])

    # predict probabilities (results are [prob_no_dup, prob_dup]) #
    probs = (clf.predict_proba(clf_inputs))[:, 1]

    return probs

###################################################################################################

def filter_duplex_prob(sam_filename, bed_filename, filter_temp, filter_prob):
    '''
    Filters probes based on user-specified duplex probability at one of six user-specified temperatures.
    Writes filtered probe sequences and duplex probabilities to new .bed file.
        Arguments:
            - sam_filename [str] : relative path to .sam file containing alignment scores for probe candidates
            - bed_filename [str] : relative path to .bed file containing sequences of final probeset
            - filter_temp [int] : one of [32, 37, 42, 47, 52, 57] specifying temperature to set probability filter at
            - filter_prob [float] : duplex probability to filter by

    '''
    # make sure files exist #
    if not os.path.exists(sam_filename):
        raise FileNotFoundError(f"The file {sam_filename} does not exist.")
    if not os.path.exists(bed_filename):
        raise FileNotFoundError(f"The file {bed_filename} does not exist.")

    # make sure files are correct type #
    if not sam_filename.endswith('.sam'):
        raise ValueError(f"The file {sam_filename} is not a SAM file.")
    if not bed_filename.endswith('.bed'):
        raise ValueError(f"The file {bed_filename} is not a BED file.")
    
    # read and check BED file format #
    with open(bed_filename) as file:
        seqs = []
        for line in file:
            parts = line.split('\t')
            if len(parts) < 5:
                raise ValueError("BED file schema is incorrect.")
            
            # check if first part is a string #
            if not isinstance(parts[0], str):
                raise ValueError("Invalid chromosome format in BED file.")

            # check if second and third parts are integers (start and end positions) #
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                raise ValueError("Start and end positions must be integers in BED file.")

            # check if fourth part contains only G, A, T, C letters (probe sequence) #
            sequence = parts[3]
            if not set(sequence).issubset(set("GATC")):
                raise ValueError("Probe sequence must contain only G, A, T, C in BED file.")

            # ensure the alignment score is a float #
            score_str = parts[4].strip() 
            try:
                score = float(score_str)
            except ValueError:
                raise ValueError("Score must be a floating-point number in BED file.")
            
            # store the sequence
            seqs.append(sequence)

    # collate probabilities #
    temps = [32, 37, 42, 47, 52, 57]

    # run the duplex prob calculation #
    all_probs = [calc_duplex_prob(sam_filename, bed_filename, temp) for temp in temps]
    all_probs = np.swapaxes(all_probs, 0, 1)

    # filter out probes which do not meet temp / prob thresholds #
    filter_temp = int(filter_temp)
    temp_to_index = {32: 0, 37: 1, 42:2, 47: 3, 52: 4, 57: 5}
    all_probs = np.array([probs for probs in all_probs if probs[temp_to_index[filter_temp]] > filter_prob])

    # for probes which passed filter, write sequences and probabilities to a .bed file #
    output_filename = bed_filename.split('.')[0] + '_pDup_filtered.bed'
    with open(output_filename, 'w') as file:
        # header line #
        file.write(f'{len(all_probs)} probes passed filtering with thresholds set to T={filter_temp}C and PDup={filter_prob} \n')
        # write as probe_number, probe_sequence, and duplex_probabilities #
        for k, seq, probs in zip(range(len(seqs)), seqs, all_probs):
            file.write(f'{k+1} \t {seq} \t {[round(prob, 8) for prob in probs]} \n')

###################################################################################################

def plot_duplex_prob(filtered_filename, probe_num='all'):
    '''
    Plots probabilities of probe forming a duplex with target sequence at all 6 temps.
        Arguments:
            - filtered_filename [str] : relative path to .txt file containing sequences and probabilities for filtered probes
            - probe_num [int] : choose to plot duplex prob for a single probe (default = 'all')
    '''
    # read in seqs and probs #
    with open(filtered_filename, 'r') as file:
        next(file)
        seqs = [line.split('\t')[1].strip() for line in file]
    with open(filtered_filename, 'r') as file:
        next(file)
        all_probs = [line.split('\t')[2:] for line in file]

    # formatting probs back to floats #
    all_probs = [probs[0].split(',') for probs in all_probs]
    all_probs = [[float(prob.strip(' [').strip('] \n')) for prob in probs] for probs in all_probs]

    # set temps
    temps = [32, 37, 42, 47, 52, 57]

    # plot all probes #
    if probe_num == 'all':
        plt.figure(figsize=(8, 6))
        for k in range(len(all_probs)):
            plt.plot(temps, all_probs[k])
        plt.xticks(temps)
        plt.ylim([-0.075, 1])
        plt.xlabel('Temperature (C)', size=15)
        plt.ylabel('Duplex Probability', size=15)
        plt.title('Probability of Probe forming Duplex with its Target Sequence', size=15)
        plt.hlines(0.2, 32, 57, label='Sufficient Binding Strength', linestyle='--', color='k')
        plt.legend()
        plt.show()

    # single probe #
    elif type(probe_num) == int:
        plt.figure(figsize=(8, 6))
        plt.plot(temps, all_probs[probe_num], label=f'Probe {probe_num}: {seqs[probe_num - 1]}')
        plt.xticks(temps)
        plt.ylim([-0.075, 1])
        plt.xlabel('Temperature (C)', size=15)
        plt.ylabel('Duplex Probability', size=15)
        plt.title('Probability of Probe forming Duplex with its Target Sequence', size=15)
        plt.hlines(0.2, 32, 57, label='Sufficient Binding Strength', linestyle='--', color='k')
        plt.legend()
        plt.show()

    else:
        raise ValueError('Probe_num argument needs to be \'all\' (to plot all probes) \
              or an int (to plot a single probe where probe_num corresponds to a line in the .bed file)')
