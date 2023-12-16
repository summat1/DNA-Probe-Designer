# DNA Probe Designer
The DNA-Probe-Designer package aims to streamline the process of probe design for in situ hybridization experiments by lowering the barrier for entry with an accessible graphical user interface and an intuitive probe design pipeline.

***NOTE: REQUIRED SETUP***<br>
The DNA-Probe-Designer package is primarily GUI-based, but does required one setup step to function properly which must be completed on the command line. This step requires installation of *bowtie2*, a software package meant for alignment of nucleic acid sequences of interest to a reference sequence (in this case, a full reference genome). Genome sequences can be downloaded at https://www.ncbi.nlm.nih.gov/datasets/genome/ (download as FASTA file format).

*bowtie2* can be installed at https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/ (remember to add bowtie2 to PATH)

Once *bowtie2* is sucessfully installed and the user has downloaded the reference genome of interest, one command must be entered into the command line to build indices which *bowtie2* will use for running the probe design pipeline:

`bowtie2-build relative/path/to/reference/genome/sequence.fa name_of_indices`

Once this step has been completed, users are ready to run the DNA-Probe-Designer package. (Note: this step will have to be completed each time the user wishes to design probes against a new reference genome).

Once this setup step is complete, the GUI can be accessed with the following Python code:

```
from DNAProbeDesigner.gui import run_gui
run_gui()
```

For help navigating the GUI, watch the demo video (Accessible in the `examples` directory on GitHub: https://github.com/summat1/DNA-Probe-Designer)

**Questions?**<br>
Any questions or issues regarding this package can be directed to longnic@uw.edu or summat@uw.edu.