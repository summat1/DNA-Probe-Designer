from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout, QWidget, QLabel, QLineEdit
import sys
import subprocess
import os

class DNAProbeDesigner(QMainWindow):
    def __init__(self):
        super().__init__()

        # Main Window setup
        self.setWindowTitle("DNA Probe Designer")
        self.setGeometry(100, 100, 550, 300)  # Adjusted for additional text

        # Layout
        layout = QVBoxLayout()
        
        # Welcome Label
        self.welcomeLabel = QLabel("Welcome to the DNA Probe Designer! To get started, please complete the 3 setup steps below:", self)
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
        self.outputDirLabel = QLabel("Output Directory: Not selected", self)
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

        # Set central widget and layout
        centralWidget = QWidget()
        centralWidget.setLayout(layout)
        self.setCentralWidget(centralWidget)

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
        if self.fastaFilePath and self.bowtieDirPath and self.bowtieIndices and self.outputDirPath :
            fastqOutputFile = os.path.join(self.outputDirPath, 'blockParse-output')
            samFile = os.path.join(self.outputDirPath, 'bowtie-output')

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
                "-x", self.bowtieIndices,  # Use the user-specified index
                "-U", f"{fastqOutputFile}.fastq",  # Input FASTQ file
                "--no-hd", "-t", "-k", "2", "--local",
                "-D", "20", "-R", "3", "-N", "1", "-L", "20",
                "-i", "C,4", "--score-min", "G,1,4",
                "-S", samFile  # Output SAM file
                ]
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error in bowtie2: {e}")
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
