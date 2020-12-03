#! /usr/bin/python3.4
from dsfilter import umifilter
import sys

umifilter(sys.stdin,sys.argv[1])
