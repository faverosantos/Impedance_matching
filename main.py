#!/usr/bin/env python3


import rflib as rf
import numpy as np
from matplotlib import pyplot as pp

if __name__ == '__main__':

    ZSOURCE = complex(20, -20)
    ZLOAD = complex(50, -100)

    f = 1E9

    dBm = 27
    mW = 137

    rf.plot_smith_chart()
    #ZA, movement = rf.add_series_inductor(ZSOURCE, 7.611 * 1E-9, f)
    #rf.plot_constant_reactance_movement(movement)
    #print("eita")

    #rf.lumped_match_impedance(ZSOURCE, ZLOAD, f)
    #ZA = rf.add_shunt_capacitor(ZLOAD, 0.886 * 1E-12, f)
    #ZSOURCE = rf.add_series_inductor(ZA, 7.611 * 1E-9, f)
    #print("ZLOAD: " + str(ZSOURCE))






