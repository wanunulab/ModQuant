#!/usr/bin/env python3
#Similar to hypermodDistance
import sys
import os
import csv
import re
import itertools
import argparse
import subprocess
#from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import time

# subprocess.Popen( f"samtools view -h -F 256 -q 10 {batchdir}/output.bam > {batchdir}/output.sorted.bam", shell=True).wait()

# with open('big_fastq.fastq', 'w') as f1:
#         for i in range(0,181):
#                 with open('FAL60812_pass_b3a7a5d6_{0}.fastq'.format(str(i)), 'r') as f2:
#                         for line in f2:
#                                 f1.write(line)

# os.system("sed -n '1~4s/^@/>/p;2~4p' big_fastq.fastq > big_fasta.fa")

target_folder = 'directory of folder' #put file path for folder you want to move files to
f5_ID = 'FAL60812_pass_b3a7a5d6' #put the correct ID file name for the fast5 files
for i in range(151):
	subprocess.Popen( f"mv {f5_ID}_{i}.fast5 {target_folder}", shell=True).wait()


