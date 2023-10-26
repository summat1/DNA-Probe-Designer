# Functional Specification

### Background
In situ hybridization is a protocol which enables spatial identification of a target nucleic acid sequence within an individual cell. This is achieved via hybridization with detectable oligonucleotides called probes. To ensure accurate results, it is crucial that a probe is highly specific to its target (i.e., does not bind elsewhere in the genome) and binds strongly to its target. Our **DNA probe designer** package allows researchers who perform in situ hybridization experiments to create DNA probes with high specificity and binding affinity.

### User Profile
Users of our package will primarily be experimentalists who plan to run in situ hybridization protocols to study the spatial context of genomic regions of interest. Users must have basic familiarity with Python including the ability to import/export files, assign variables, and use functions.

### Use Cases
##### 1. De Novo Probe Design
<sup>
Designing probes from the ground up for a given target sequence is useful for researchers studying genomic regions for which probes do not already exist. In this case, the user is seeking a sequence of DNA nucleotides that will act as an effective probe for the specified target. The user would input the sequence of their target region and file(s) containing the full relevant genome and receive sequences corresponding to the top probe candidates and their associated performance metrics (specificity, binding affinity, etc.)
</sup>

##### 2. Quantification of Existing Probe Performance
Researchers may wish to quantify the performance of an existing probe as a sanity check before running an experiment or to determine if poor probe performance is causing off-target binding or low signal in their protocols.