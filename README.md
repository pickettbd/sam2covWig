# sam2covWig

Create coverage Wiggle tracks from SAM-formatted alignments.

## Installation

This is currently implemented as a simple Python script. As such, "installation" is as simple as adding it to a location in your `$PATH` or providing the path to the script as an argument to Python.

### Dependencies

Python v3.6+ (because of f-strings) is required. This implementation depends only on standard modules, specifically: argparse, collections.deque, math.ceil, re, and sys.

## Running

For a details regarding input/output and options, please run the program with `-h` or `--help`. Generally, use of this program should be proceed as follows:

```bash
samtools view input.bam | python3 ./sam2covWig.py > output.coverage.wig
```

## Citation
If you find this script useful in your work, please cite this repository.

###### Last updated: 2 March 2022
