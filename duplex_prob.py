import numpy as np
import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from Bio.SeqUtils import GC

###################################################################################################

def calc_duplex_prob(sam_filename, probeset_filename, temp):
    '''
    Calculates probability of all probes in a probeset forming a duplex with target sequence at a given temp.
        Arguments:
            - sam_filename [str] : relative path to .sam file containing alignment scores for probe candidates
            - probeset_filename [str] : relative path to .bed file containing sequences of final probeset
            - temp [int] : temp at which to predict duplex probability
        Outputs:
            - seqs [np.ndarray] : sequences of probes in probeset
            - probs [np.ndarray] : duplex probabilities of probes in probeset
    '''
    # read in final probeset BED file #
    with open(probeset_filename) as file:
        probeset = [line.split('\t')[3] for line in file]

    # read in SAM file #
    with open(sam_filename) as file:
        sam = [line for line in file]

    # set up LDA model to predict duplex probability based on probe length, GC content, and alignment score to target seq #
    # Values from models published in Beliveau, et al. (2018) paper. #
    temps = np.array([32, 37, 42, 47, 52, 57])
    if temp not in temps:
        print(f'hey asshole, temp has to be one of {temps}')
        return None
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

    # collate sequences and probs
    seqs = np.array([seq for seq in probeset])

    return seqs, probs