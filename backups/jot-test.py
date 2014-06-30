#!/usr/bin/python

import os
os.system("sed -e '/W        W/d' -e '/ION  /d' adhesion-tests/membrane-start.gro > donkey")
