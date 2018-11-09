#!/usr/bin/env python

import os
import subprocess
import re


def run_sdlog_dump(file_path):
    print file_path
    file_out = re.sub('.bin', '.csv', file_path)
    cmd = 'python sdlog2_dump.py {file_path:s} -f {file_out:s}'.format(**locals())
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()

for root, dirs, files in os.walk(os.curdir):
        for file in files:
            if file.endswith('.bin'):
                file_path = os.path.abspath(os.path.join(root, file))
                run_sdlog_dump(file_path)
