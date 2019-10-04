# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
import os
os.system("python3 start_simulation.py params.ini > out.txt")
with open("out.txt", "r") as out:
    for line in out:
        print(line)
        