from PyQt6.QtWidgets import QApplication, QMainWindow, QTabWidget, QPushButton, QFileDialog, QVBoxLayout, QWidget, QLabel, QLineEdit, QComboBox
from PyQt6.QtGui import QDoubleValidator
import sys
import subprocess
import os

class DNAProbeDesigner(QMainWindow):
    def __init__(self):
        super().__init__()

        # Main Window setup
        self.setWindowTitle("DNA Probe Designer")
        self.setGeometry(100, 100, 550, 300)

        # Main tab widget
        self.tabs = QTabWidget(self)
        self.setCentralWidget(self.tabs)

        # Create tabs
        self.dnaProbeTab = QWidget()
        self.plottingTab = QWidget()

        # Add tabs
        self.tabs.addTab(self.dnaProbeTab, "DNA Probe Design")
        self.tabs.addTab(self.plottingTab, "Analysis and Plotting")

        # Setup each tab
        self.setupDnaProbeTab()
        self.setupPlottingTab()

    def setupDnaProbeTab(self):
        layout = QVBoxLayout(self.dnaProbeTab)

        # Welcome Label
        self.welcomeLabel = QLabel("Welcome to the DNA Probe Designer! To get started, please complete the 4 setup steps below:", self)
        self.welcomeLabel.setWordWrap(True)

        # Step 1 - Select FASTA Attachment Label + Button
        self.fastaLabel = QLabel("Step 1: Attach target sequence.", self)
        self.fastaLabel.setWordWrap(True)
        self.attachFastaFileBtn = QPushButton("Attach .fasta", self)
        self.attachFastaFileBtn.clicked.connect(self.openFileNameDialog)

        # Step 2 - Specify Bowtie Indices Label + Button
        self.bowtieLabel = QLabel("Step 2: Select path to Bowtie indices.", self)
        self.bowtieLabel.setWordWrap(True)
        self.selectBowtieDir = QPushButton("Specify path to indices", self)
        self.selectBowtieDir.clicked.connect(self.selectBowtieDirectory)

        # Step 3 - Specify Bowtie Filename
        self.bowtieName = QLabel("Step 3: Specify Bowtie filename.", self)
        self.bowtieName.setWordWrap(True)
        self.typeBowtieName = QLineEdit(self)
        self.typeBowtieName.setPlaceholderText("Enter Bowtie2 name")
        self.typeBowtieName.textChanged.connect(self.updateBowtieIndexName)

        # Step 4 - Select Output Directory Label + Button
        self.outputDirLabel = QLabel("Step 4: Select Output Directory.", self)
        self.outputDirLabel.setWordWrap(True)
        self.outputDirBtn = QPushButton("Select Output Directory", self)
        self.outputDirBtn.clicked.connect(self.selectOutputDirectory)
        
        # Variable to store the paths
        self.fastaFilePath = ""
        self.bowtieDirPath = ""
        self.bowtieIndices = ""
        self.outputDirPath = ""
    
        # Run Button
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

        # Status Label
        self.statusLabel = QLabel("", self)

        # Adding widgets to the layout
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

    def setupPlottingTab(self):
        layout = QVBoxLayout(self.plottingTab)

        # Dropdown for temperature
        self.temperatureDropdown = QComboBox(self.plottingTab)
        self.temperatureDropdown.addItems(['37', '42', '47', '52', '57'])
        layout.addWidget(self.temperatureDropdown)

        # Input for probability
        self.probabilityInput = QLineEdit(self.plottingTab)
        self.probabilityInput.setValidator(QDoubleValidator(0.0, 1.0, 10))
        self.probabilityInput.setPlaceholderText("Enter probability (0 - 1)")
    
        # Run button
        self.runButton = QPushButton("Filter Probes", self.plottingTab)
        #self.runButton.clicked.connect(self.filterProbes)
    
        # Plot button
        self.plotButton = QPushButton("Plot Results", self.plottingTab)
        #self.plotButton.clicked.connect(self.plotResults)
    
        layout.addWidget(self.probabilityInput)
        layout.addWidget(self.runButton)
        layout.addWidget(self.plotButton)

    def openFileNameDialog(self):
        fastaFileName, _ = QFileDialog.getOpenFileName(self, "Select a FASTA file", "", "FASTA Files (*.fasta);;All Files (*)")
        if fastaFileName:
            self.fastaFilePath = fastaFileName
            self.fastaLabel.setText(f"FASTA attachment: {fastaFileName}")

    def selectBowtieDirectory(self):
        bowtieDir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if bowtieDir:
            self.bowtieDirPath = bowtieDir
            self.bowtieLabel.setText(f"Bowtie Directory: {bowtieDir}")
    
    def updateBowtieIndexName(self, text):
        self.bowtieIndices = text

    def selectOutputDirectory(self):
        outputDir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if outputDir:
            self.outputDirPath = outputDir
            self.outputDirLabel.setText(f"Output Directory: {outputDir}")

    def runScript(self):
        if self.fastaFilePath and self.bowtieDirPath and self.bowtieIndices and self.outputDirPath:

            # path to the fastq file
            fastqOutputFile = os.path.join(self.outputDirPath, 'blockParse_output').replace('\\', '/')

            # path to the sam file
            samFile = os.path.join(self.outputDirPath, 'bowtie_output.sam').replace('\\', '/')

            # path to the folder of indices
            bowtiePathArgument = os.path.join(self.bowtieDirPath, self.bowtieIndices).replace('\\', '/')

            # path to the bed file
            bedFile = os.path.join(self.outputDirPath, 'bowtie_output_probes.bed').replace('\\', '/')

            # fastq file generation
            try:
                command = [
                    "python", "Oligominer/blockParse.py",
                    "-f", self.fastaFilePath,
                    "-o", fastqOutputFile,
                    # other arguments could go here
                ]
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error: {e}")

            try:
                command = [
                "bowtie2",
                "-x", bowtiePathArgument,  # Use the user-specified index
                "-U", f"{fastqOutputFile}.fastq",  # Input FASTQ file
                "--no-hd", "-t", "-k", "2", "--local",
                "-D", "20", "-R", "3", "-N", "1", "-L", "20",
                "-i", "C,4", "--score-min", "G,1,4",
                "-S", samFile  # Output SAM file
                ]
                subprocess.run(command, check=True, shell=True)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error in bowtie2: {e}")

            try:
                command = [
                    "python", "Oligominer/outputClean.py",
                    "-T", "42",
                    "-f", samFile
                    # other arguments could go here
                ]
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error: {e}")

            try:
                command = [
                    "python", "Oligominer/structureCheck.py",
                    "-f", bedFile,
                    "-t", "0.4"
                    # other arguments could go here
                ]
                subprocess.run(command, check=True)
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
            
# Main loop
app = QApplication(sys.argv)
mainWindow = DNAProbeDesigner()
mainWindow.show()
sys.exit(app.exec())
