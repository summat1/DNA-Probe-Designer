## Software components.

High level description of the software components such as: *data manager*, which provides a simplified interface to your data and provides application specific features (e.g., querying data subsets); and *visualization manager*, which displays data frames as a plot. Describe at least 3 components specifying: what it does, inputs it requires, and outputs it provides.

To begin using our package, researchers will need to prepare a FASTA file with the sequence of the gene that they are interested in attaching a probe to. We will parse this file and use it to create a **targetSequence object** - this will contain metadata about the target sequence that the researchers are interested in. Once the sequence information is captured, we will calculate and store various statistics such as the length and GC content of the sequence. Furthermore, we can run various methods like *findAlignment* or *findOffSiteBinding* to generate and refine potential *probe* DNA sequences.

A potential DNA probe will be encoded as a **probeCandidate object**. This class will store critical information about a probe such as:
* Expected levels of off target binding
* Amount of GC content
* Binding affinity strength
* PDupe

Each **probeCandidate** will have various methods to display both numerical/tabular and graphical data about the probe. We will implement methods such as *showPDupeGraph* and *displayOffTargetBinding* that will provide the relevant context about the candidate. 

We plan to develop a GUI so researchers do not have to call these methods themselves - rather the GUI interface will call them in the backend when the user specifies what information they want. 

## Interactions to accomplish use cases. 

Describe how the above software components interact to accomplish at least one of your use cases.

## Preliminary plan. 

A list of tasks in priority order.
