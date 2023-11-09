# Component Specification

### Software Components
**Probe Design Pipeline (via OligoMiner)**<br>
Responsible for generation of probe sequences. Occurs in two primary phases: generation of iniital candidates based on thermodynamic properties followed by filtering of candidate probes based on specificity.

This component takes as input a FASTA file specifying a sequence which the user wants to design probes against and returns as output the a set of probe sequences with associated metadata (melting temperature, specificity, etc.)

*Note: This component utilizes OligoMiner, a set of scripts for designing DNA probes.*

**Visualization Suite**<br>
The probe design pipeline will return a set of probes from which the user can select, each with different properties. This component will summarize these properties in a series of easy to interpret plots which the user can utilize to compare each probe and make a decision on which to use for their experiment.

This component takes as input the set of probes returned by the probe design pipeline and the associated metadata and returns as output several visualizations such as duplex probability vs. temperature, quantification of specificity, etc.

**User Interface**<br>
To ensure accessibility, this component provides the user with a graphical interface which will allow for easy uploading of input files, tuning of parameters, and visual monitoring of the pipeline's status.

This component does not have any inputs/outputs as it is essentially a GUI wrapper for the probe design pipeline.

<br>

### Interactions to Accomplish Use Case
For a user wishing to design DNA probes from scratch, the probe design pipeline and visualization suite interact in series, with the outputs of the pipeline (the probe sequences + metadata) being fed into the visualization suite which in turn generates the relevant plots. The user interface acts as a wrapper for the other two components, running in parallel with the pipeline and housing the visualizations.

<br>

### Preliminary Plan
1. Implement probe design pipeline using existing scripts from OligoMiner
2. Decide which visualizations will be included in the final output
3. Implement visualization suite using test outputs from probe design pipeline
4. Wrap up everything nicely in a GUI using PyQT