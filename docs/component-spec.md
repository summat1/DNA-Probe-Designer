# Component Specification

### First Steps
- start with single chromosome (reduced computation time for easy testing)
- initial implementation plan
    1. define target region (maybe coding region +/- a few BP)
    2. sliding window through target region; for each window...
        - calculate # of times that sequence appears in rest of chromosome with exact sequence, one BP different, two BP different, etc.

- will also be helpful to understand how genome files are structured (i.e., what info does a FASTA file contain)