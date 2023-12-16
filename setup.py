from setuptools import setup, find_packages

def setup_config(dependencies):

    with open('README.md', 'r') as file:
        long_descrip = file.read()

    setup(
        name='DNAProbeDesigner',
        version='0.0.19',
        packages=find_packages('DNAProbeDesigner'),
        install_requires=dependencies,
        description='Design DNA probes for ISH experiments.',
        long_description=long_descrip,
        long_description_content_type='text/markdown',
    )

dependencies = [
                'biopython',
                'matplotlib',
                'numpy',
                'pyqt6',
                'scikit-learn',
                'seqfold'
]

setup_config(dependencies)