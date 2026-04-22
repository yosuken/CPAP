#!/bin/bash
set -euo pipefail

micromamba env create -f environment.yaml
micromamba run -n CPAP_v0.3.1 bundle install
