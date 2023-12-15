# Component Specification

### Software Components
**Probe Design Pipeline (via OligoMiner)**<br>
Responsible for generation of probe sequences. Occurs in two primary phases: generation of initial candidates based on thermodynamic properties followed by filtering of candidate probes based on specificity.

This component takes as input a FASTA file specifying a sequence which the user wants to design probes against and returns as output the a set of probe sequences with associated metadata (melting temperature, specificity, etc.)

*Note: This component utilizes OligoMiner, a set of scripts for designing DNA probes.*


**duplex_prob.py**<br>
This module allows users to filter the probeset output by OligoMiner based on the probability of a probe forming a duplex with its target sequence at 6 different temperatures (32, 37, 42, 47, 52, and 57 C). Users can also visualize the result of the filtering by plotting the duplex probability at each of the 6 temperatures for each probe which passed the filtering step.

*calc_duplex_prob*<br>
Duplex probabilities are calculated using a Linear Discriminant Analysis model which takes probe length, GC content, and alignment score with target sequence as inputs and returns a probability of the probe being in the bound or unbound state.

*filter_duplex_prob*<br>
Once duplex probabilities are calculated at each temperature for each probe, they are filtered based on a user-specified duplex probability for a given temperature. For example, a user can input a duplex probability of 0.2 at a temperature of 42C and all probes which have duplex probabilities below 0.2 at 42C will be filtered out.

*plot_duplex_prob*<br>
Once the filtering step has been applied, users can simply press a button in the GUI to create a plot of duplex probability vs temperature for each of the probes which passed the filter, both as a helpful visualization and as a sanity check to ensure that the filtered probeset has the desired characteristics.


**secondary_structure.py**<br>
This module allows users to further filter their probeset based on the minimum free energy of each probe. Nucleic acid sequences can adopt a number of secondary structures which can interfere with duplex formation with the target sequence. The stability of each of these secondary structures is quantified by its free energy; a probe's minimum free energy (MFE) is the free energy of the most stable secondary structure. Lower MFE values indicate a more stable secondary structure and thus a higher chance of poor duplex formation.

*filter_secondary_structure*<br>
MFE values are calculated for each probe via the seqfold package which implements Zuker's dynamic programming algorithm for MFE calculation. Users may specify a MFE threshold in the GUI; only probes with MFE values above this threshold will be passed.


**Graphical User Interface**<br>
To aid users who are not familiar with Python and/or have no programming experience, we have created a graphical user interface (GUI) which is easy to navigate and requires zero coding to use. Users can upload input files, choose filtering parameters, and plot visualizations simply by clicking buttons and navigating file browsers. The GUI also enables visual monitoring of the pipeline's status and prints helpful tips to the user if inputs are not correct.

The GUI was implemented using PyQt6, a package which uses object-oriented programming to allow for streamlined GUI development.

### Interactions to Accomplish Use Case
For a user wishing to design DNA probes from scratch, the probe design pipeline (OligoMiner) and filtering and visualization suite (duplex_prob.py and secondary_structure.py) interact in series, with the outputs of the pipeline (the probe sequences + metadata) being fed into the filtering and visualization suite which in turn generates the relevant plots. The user interface acts as a wrapper for the other components, running in parallel with the pipeline and housing the filtering and visualization software.