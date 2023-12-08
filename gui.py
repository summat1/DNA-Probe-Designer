from PyQt6.QtWidgets import QApplication, QMainWindow, QTabWidget, QPushButton, QFileDialog, QVBoxLayout, QWidget, QLabel, QLineEdit, QComboBox, QProgressBar
from PyQt6.QtGui import QDoubleValidator
from duplex_prob import filter_duplex_prob, plot_duplex_prob
from secondary_structure import filter_secondary_structure
import sys
import subprocess
import os

# graphical user interface (gui) for DNA probe design
class DNAProbeDesigner(QMainWindow):
    def __init__(self):
        super().__init__()

        # main window setup
        self.setWindowTitle("DNA Probe Designer")
        self.setGeometry(100, 100, 550, 450)

        # main tab widget
        self.tabs = QTabWidget(self)
        self.setCentralWidget(self.tabs)

        # two tabs
        self.dnaProbeTab = QWidget()
        self.plottingTab = QWidget()
        self.tabs.addTab(self.dnaProbeTab, "DNA Probe Design")
        self.tabs.addTab(self.plottingTab, "Analysis and Plotting")
        self.setupDnaProbeTab()
        self.setupPlottingTab()

    # function to set up the DNA probe design tab
    def setupDnaProbeTab(self):
        layout = QVBoxLayout(self.dnaProbeTab)

        # welcome label
        self.welcomeLabel = QLabel("Welcome to the DNA Probe Designer! To get started, please complete the 4 setup steps below:", self)
        self.welcomeLabel.setWordWrap(True)

        # step 1 - Select FASTA Attachment Label + Button
        self.fastaLabel = QLabel("Step 1: Attach target sequence.", self)
        self.fastaLabel.setWordWrap(True)
        self.attachFastaFileBtn = QPushButton("Attach .fasta", self)
        self.attachFastaFileBtn.clicked.connect(self.openFileNameDialog)

        # step 2 - Specify Bowtie Indices Label + Button
        self.bowtieLabel = QLabel("Step 2: Select path to Bowtie indices.", self)
        self.bowtieLabel.setWordWrap(True)
        self.selectBowtieDir = QPushButton("Specify path to indices", self)
        self.selectBowtieDir.clicked.connect(self.selectBowtieDirectory)

        # step 3 - Specify Bowtie Filename
        self.bowtieName = QLabel("Step 3: Specify Bowtie filename.", self)
        self.bowtieName.setWordWrap(True)
        self.typeBowtieName = QLineEdit(self)
        self.typeBowtieName.setPlaceholderText("Enter Bowtie2 name")
        self.typeBowtieName.textChanged.connect(self.updateBowtieIndexName)

        # step 4 - Select Output Directory Label + Button
        self.outputDirLabel = QLabel("Step 4: Select Output Directory.", self)
        self.outputDirLabel.setWordWrap(True)
        self.outputDirBtn = QPushButton("Select Output Directory", self)
        self.outputDirBtn.clicked.connect(self.selectOutputDirectory)
        
        # progress bar
        self.progressBar = QProgressBar(self)
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.setStyleSheet("""
            QProgressBar {
                border: 2px solid grey;
                border-radius: 5px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #05B8CC;
                width: 10px; /* size of the chunks */
                margin: 0.5px;
            }
        """)
       
        # variables to store the paths so that we can change them and know when they have been changed
        self.fastaFilePath = ""
        self.bowtieDirPath = ""
        self.bowtieIndices = ""
        self.outputDirPath = ""
        self.filteredProbeFile = ""
        self.mfeFilteredProbeFile = ""
        self.samFile = ""
        self.bedFile = ""
        self.filterTemp = ""
        self.filterProb = ""
        self.filterMFE = ""
    
        # run probe design
        self.runBtn = QPushButton("Run Probe Design", self)
        self.runBtn.setStyleSheet("""
            QPushButton {
                font-weight: bold;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #ADD8E6;
            }""")
        self.runBtn.clicked.connect(self.runScript)

        # status to inform user on what the issue is
        self.statusLabel = QLabel("", self)

        # adding widgets to the layout
        layout.addWidget(self.welcomeLabel)
        layout.addWidget(self.fastaLabel)
        layout.addWidget(self.attachFastaFileBtn)
        layout.addWidget(self.bowtieLabel)
        layout.addWidget(self.selectBowtieDir)
        layout.addWidget(self.bowtieName)
        layout.addWidget(self.typeBowtieName)
        layout.addWidget(self.outputDirLabel)
        layout.addWidget(self.outputDirBtn)
        layout.addWidget(self.runBtn)
        layout.addWidget(self.statusLabel)
        layout.addWidget(self.progressBar)

    # function to setup the plotting and analysis tab
    def setupPlottingTab(self):
        layout = QVBoxLayout(self.plottingTab)

        # dropdown for temperature
        self.temperatureLabel = QLabel("Select Desired Temperature Cutoff:", self)
        self.temperatureLabel.setWordWrap(True)
        self.temperatureDropdown = QComboBox(self.plottingTab)
        self.temperatureDropdown.addItems(['Select', '37', '42', '47', '52', '57'])
        self.temperatureDropdown.currentTextChanged.connect(self.updateSelectedTemperature)

        # input for probability
        self.probabilityLabel = QLabel("Select Desired Duplex Probability:", self)
        self.probabilityLabel.setWordWrap(True)
        self.probabilityInput = QLineEdit(self.plottingTab)
        self.probabilityInput.setValidator(QDoubleValidator(0.0, 1.0, 10))
        self.probabilityInput.setPlaceholderText("Enter probability (0 - 1)")
        self.probabilityInput.textChanged.connect(self.updateEnteredProbability)

        # input for MFE 
        self.mfeLabel = QLabel("Enter Minimum Free Energy Threshold (optional):", self)
        self.mfeLabel.setWordWrap(True)
        self.mfeInput = QLineEdit(self.plottingTab)
        self.mfeInput.setValidator(QDoubleValidator())
        self.mfeInput.setPlaceholderText("Enter MFE")
        self.mfeInput.textChanged.connect(self.updateEnteredMFE)

        # filter probes
        self.runButton = QPushButton("Filter Probes", self.plottingTab)
        self.runButton.clicked.connect(self.runFilterProbes)
    
        # plot button
        self.plotButton = QPushButton("Plot Results", self.plottingTab)
        self.plotButton.setStyleSheet("""
            QPushButton {
                font-weight: bold;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #FF7F7F;
            }""")
        self.plotButton.clicked.connect(self.plotResults)
        
        # status label
        self.plotStatusLabel = QLabel("", self)

        # add widgets to layout
        layout.addWidget(self.temperatureLabel)
        layout.addWidget(self.temperatureDropdown)
        layout.addWidget(self.probabilityLabel)
        layout.addWidget(self.probabilityInput)
        layout.addWidget(self.mfeLabel)
        layout.addWidget(self.mfeInput)
        layout.addWidget(self.runButton)
        layout.addWidget(self.plotButton)
        layout.addWidget(self.plotStatusLabel)

    # select fasta file dialog
    def openFileNameDialog(self):
        fastaFileName, _ = QFileDialog.getOpenFileName(self, "Select a FASTA file", "", "FASTA Files (*.fasta);;All Files (*)")
        if fastaFileName:
            self.fastaFilePath = fastaFileName
            self.fastaFileName = os.path.basename(self.fastaFilePath).split('.')[0]
            self.fastaLabel.setText(f"FASTA attachment: {fastaFileName}")

    # select bowtie indices dialog
    def selectBowtieDirectory(self):
        bowtieDir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if bowtieDir:
            self.bowtieDirPath = bowtieDir
            self.bowtieLabel.setText(f"Bowtie Directory: {bowtieDir}")
    
    # store the data when the index name changes
    def updateBowtieIndexName(self, text):
        self.bowtieIndices = text

    # select output dialog
    def selectOutputDirectory(self):
        outputDir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if outputDir:
            self.outputDirPath = outputDir
            self.outputDirLabel.setText(f"Output Directory: {outputDir}")

    # progress bar movement
    def updateProgressBar(self, value):
        self.progressBar.setValue(value)

    # run the initial probe design process
    def runScript(self):
        if self.fastaFilePath and self.bowtieDirPath and self.bowtieIndices and self.outputDirPath:
            
            # path to the fastq file
            fastqOutputFile = os.path.join(self.outputDirPath, self.fastaFileName).replace('\\', '/')
            # path to the sam file
            self.samFile = os.path.join(self.outputDirPath, f'{self.fastaFileName}.sam').replace('\\', '/')
            # path to the folder of indices
            bowtiePathArgument = os.path.join(self.bowtieDirPath, self.bowtieIndices).replace('\\', '/')
            # path to the bed file
            self.bedFile = os.path.join(self.outputDirPath, f'{self.fastaFileName}_probes.bed').replace('\\', '/')

            # fastq file generation
            try:
                command = [
                    "python", "Oligominer/blockParse.py",
                    "-f", self.fastaFilePath,
                    "-o", fastqOutputFile,
                    # other arguments could go here
                ]
                self.updateProgressBar(33)
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error: {e}")

            try:
                command = [
                "bowtie2",
                "-x", bowtiePathArgument,  # indices specified
                "-U", f"{fastqOutputFile}.fastq",  # fastq file
                "--no-hd", "-t", "-k", "2", "--local",
                "-D", "20", "-R", "3", "-N", "1", "-L", "20",
                "-i", "C,4", "--score-min", "G,1,4",
                "-S", self.samFile  # output SAM file
                ]
                self.updateProgressBar(66)
                subprocess.run(command, check=True, shell=True)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error in bowtie2: {e}")

            try:
                command = [
                    "python", "Oligominer/outputClean.py",
                    "-T", "42",
                    "-f", self.samFile
                    # other arguments could go here
                ]
                subprocess.run(command, check=True)
                self.updateProgressBar(100)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error: {e}")

        else:
            if not self.fastaFilePath:
                self.statusLabel.setText("No Sequence Selected!")
            elif not self.bowtieDirPath:
                self.statusLabel.setText("Bowtie Indices Directory Not Selected!")
            elif not self.bowtieIndices:
                self.statusLabel.setText("Bowtie Indices Name Not Specified!")
            elif not self.outputDirPath:
                self.statusLabel.setText("Output Directory Not Selected!")

    # update temp if changed
    def updateSelectedTemperature(self, text):
        self.filterTemp = text

    # update Pdupe if changed
    def updateEnteredProbability(self, text):
        try:
            self.filterProb = float(text) if text else None
        except ValueError:
            self.filterProb = None
            self.probabilityInput.setPlaceholderText("Please Enter A Float Between 0 and 1")
    
    # update MFE if changed
    def updateEnteredMFE(self, text):
        try:
            self.filterMFE = float(text) if text else None
        except ValueError:
            self.filterMFE = None
            self.mfeInput.setPlaceholderText("Please Enter A Float Between 0 and 1")

    # filter the probes, differently if the MFE is specified or not
    def runFilterProbes(self):
        if self.samFile and self.bedFile and self.filterTemp and self.filterProb:
            filter_duplex_prob(self.samFile, self.bedFile, self.filterTemp, self.filterProb)
            self.filteredProbeFile = f'{self.bedFile.split(".")[0]}_pDup_filtered.bed'
            if self.filterMFE:
                filter_secondary_structure(self.filteredProbeFile, self.filterMFE)
                self.mfeFilteredProbeFile = f'{self.bedFile.split(".")[0]}_pDup_MFE_filtered.bed'
            self.plotStatusLabel.setText("Filtered Probes Successfully!")
        elif not self.samFile or not self.bedFile:
            self.plotStatusLabel.setText("Please Design Probes First!")
        elif not self.filterTemp:
            self.plotStatusLabel.setText("Please Select Temp!")
        elif not self.filterProb:
            self.plotStatusLabel.setText("Please Select Duplex Probability!")
    
    # plot the filteres probes
    def plotResults(self):
        if self.filterMFE and self.mfeFilteredProbeFile:
            plot_duplex_prob(self.mfeFilteredProbeFile)
        else:
            plot_duplex_prob(self.filteredProbeFile)    
    
# main loop
app = QApplication(sys.argv)
mainWindow = DNAProbeDesigner()
mainWindow.show()
sys.exit(app.exec())
