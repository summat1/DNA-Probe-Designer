from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout, QWidget, QLabel
import sys
import subprocess
import os

class DNAProbeDesigner(QMainWindow):
    def __init__(self):
        super().__init__()

        # Main Window setup
        self.setWindowTitle("DNA Probe Designer")
        self.setGeometry(100, 100, 400, 300)  # Adjusted for additional text

        # Layout
        layout = QVBoxLayout()
        
        # Welcome Label
        self.welcomeLabel = QLabel("Welcome to the DNA Probe Designer!", self)
        self.welcomeLabel.setWordWrap(True)

        # Description Label
        self.descriptionLabel = QLabel("Attach a FASTA file and then click 'Run' to proceed with Oligominer.", self)
        self.descriptionLabel.setWordWrap(True)

        # Attach File Button
        self.attachFileBtn = QPushButton("Attach File", self)
        self.attachFileBtn.clicked.connect(self.openFileNameDialog)

        # In the __init__ method of DNAProbeDesigner
        self.outputDirBtn = QPushButton("Select Output Directory", self)
        self.outputDirBtn.clicked.connect(self.selectOutputDirectory)
        layout.addWidget(self.outputDirBtn)

        # Add a label to display the selected directory
        self.outputDirLabel = QLabel("Output Directory: Not selected", self)
        layout.addWidget(self.outputDirLabel)

        # Variable to store the output directory path
        self.outputDirPath = ""
    
        # Run Button
        self.runBtn = QPushButton("Run", self)
        self.runBtn.clicked.connect(self.runScript)

        # Status Label
        self.statusLabel = QLabel("", self)

        # Adding widgets to the layout
        layout.addWidget(self.welcomeLabel)
        layout.addWidget(self.descriptionLabel)
        layout.addWidget(self.attachFileBtn)
        layout.addWidget(self.runBtn)
        layout.addWidget(self.statusLabel)

        # Set central widget and layout
        centralWidget = QWidget()
        centralWidget.setLayout(layout)
        self.setCentralWidget(centralWidget)

        # Variable to store the file path
        self.filePath = ""

    def openFileNameDialog(self):
        fileName, _ = QFileDialog.getOpenFileName(self, "Select a FASTA file", "", "FASTA Files (*.fasta);;All Files (*)")
        if fileName:
            self.filePath = fileName

    def selectOutputDirectory(self):
        outputDir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if outputDir:
            self.outputDirPath = outputDir
            self.outputDirLabel.setText(f"Output Directory: {outputDir}")


    def runScript(self):
        if self.filePath and self.outputDirPath:
            outputFile = os.path.join(self.outputDirPath, 'blockParse-output')
            try:
                command = [
                    "python", "Oligominer/blockParse.py",
                    "-f", self.filePath,
                    "-o", outputFile,
                    # Add other necessary arguments
                ]
                subprocess.run(command, check=True)
                self.statusLabel.setText("Run Successful. Please check Output Folder specified for FASTQ file.")
            except subprocess.CalledProcessError as e:
                self.statusLabel.setText(f"Error: {e}")
        else:
            if not self.filePath:
                self.statusLabel.setText("No File Selected!")
            elif not self.outputDirPath:
                self.statusLabel.setText("Output Directory Not Selected!")

        

# Main loop
app = QApplication(sys.argv)
mainWindow = DNAProbeDesigner()
mainWindow.show()
sys.exit(app.exec())
