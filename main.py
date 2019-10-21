#!/usr/bin/env python3


import rflib as rf
import numpy as np
from matplotlib import pyplot as pp


if __name__ == '__main__':

    ZSOURCE = complex(20, -20)
    ZLOAD = complex(50, -100)
    f = 2.4E9

    dBm = 27
    mW = 137
    rf.lumped_match_impedance(ZSOURCE,ZLOAD,f)
    #print(mW, " mW equals to ", rf.mW_2_dBm(mW), " dBm")
    #print(dBm, " dBm equals to ", rf.dBm_2_mW(dBm), " mW")

    #rf.example()


    #FS = 1*1E6
    #TS = 1/FS
    #BW = 100
    #c = 3*1E8
    #f0 = 10
    #t_atraso = 1*(0.1/1E2)
    #print("t_atraso:", t_atraso)
    #n_atraso = int(round(t_atraso/TS,1))
    #print("n_atraso:", n_atraso)

    #print(str(TS*5))

    #S = 1
    #my_chirp_tx = rf.create_chirp(f0, BW, FS, S)
    #pp.plot(my_chirp_tx)

    #TC = TS * len(my_chirp_tx)
    #print("TC: ", str(TC))
    #S = BW / TC
    #print("S: ", str(S))

    #my_chirp_rx = rf.delay_a_chirp(my_chirp_tx, 10)
    #pp.plot(my_chirp_rx)

    #mixed_signal = rf.mix_these_signals(my_chirp_tx, my_chirp_rx)

    #mixed_signal = mixed_signal[n_atraso:int(2*n_atraso):]
    #mean = np.average(mixed_signal)
    #pp.plot(mixed_signal)


    #mixed_signal = mixed_signal - mean

    #pp.plot(mixed_signal)
    #pp.show(block=False)

    #rf.fFtof(my_chirp_tx, FS)
    #rf.fFtof(my_chirp_rx, FS)
    #fft_my_signal = rf.fFtof(mixed_signal, FS)
    #max_val = max(abs(fft_my_signal))
    #abs_fft_my_signal = abs(fft_my_signal)
    #fft_my_signal_list = (abs_fft_my_signal.tolist())
    #max_idx = fft_my_signal_list.index(max_val)
    #fi = max_idx
    #print("fi: ", str(fi))
    #S = BW/TC
    #print("S: ", str(S))
    #d = c*fi/(2*S)
    #print("d: ", str(d))
    #pp.show()






