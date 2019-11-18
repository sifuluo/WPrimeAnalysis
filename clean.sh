#!/bin/bash

rm slurmlog/*
rm results/hypothesis*
hadd results/hypothesis_BG.root results/massresult/hypothesis_BG*
hadd results/hypothesis_FL.root results/massresult/hypothesis_FL*
hadd results/hypothesis_LL.root results/massresult/hypothesis_LL*
rm results/massresult/*.root
