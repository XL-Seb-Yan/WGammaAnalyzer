#!/bin/bash

root -l -q analyzer.C+\(\"samples_analyzer.conf\"\)

rm *.so *.d *.pcm
