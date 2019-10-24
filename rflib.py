# rflib: a radio frequency support library
# A code by Fávero Santos
# 13/03/2019

import numpy as np
from smithplot import SmithAxes
from matplotlib import rcParams, pyplot as pp
rcParams.update({"legend.numpoints": 3})

def fFtof(signal, FS):

    fft = np.fft.fft(signal, 1000000)
    f = np.fft.fftfreq(len(fft))
    half_fft = fft[:len(fft) // 2]
    half_f = FS*f[:len(f) // 2]

    pp.figure()
    pp.plot(half_f, np.abs(half_fft))
    pp.show(block = False)

    return fft

# This algorithm was developed based on the definition of the Fourier Transform. It was inspired by the numeric algorithm that prof. Dartora from UFPR teaches during microwave classes.
def F_ft():

    t0 = -10
    tf = 10
    Nt = 1001
    t = np.linspace(t0, tf, Nt)
    dt = t[2] - t[1]
    tau = 0.5
    f = (t >= (-tau/2))*((t <= (tau/2)))*1
    fs = 1/dt
    Nw = 1001
    w = 2 *np.pi*np.linspace(-fs/20, fs/20, Nw)
    fw = np.zeros(Nw)
    for m in range(0, Nw):
        wm = w[m]
        integrando = np.multiply(f, np.exp(-1j * wm * t))
        fw[m] = np.trapz(integrando, t)

    pp.plot(w, abs(fw))
    pp.show()

    return fw

def mix_these_signals(signal1, signal2):
    return np.multiply(signal1,signal2)
    #return signal1 - signal2

# Creates a chirp from starting from frequency f0 and ending on frequency f0+BW
# The chirps are composed by concatenated sinusoids. The sinusoids are created with a standard 10E3 sampling frequency.
def create_chirp(f0, BW, FS, slope):

    amplitude = 1
    frequency = f0
    sampling_frequency = FS
    cycles = 1
    phase = 0
    ff = f0+BW
    #f_chirp = np.arange(f0, final_freq)

    for f in range(int(f0), int(ff)):
        sine, time = generate_sine_wave(amplitude, sampling_frequency, slope*f, cycles, phase)

        if f == f0:
            chirp = []

        chirp = np.concatenate((chirp[:len(chirp)-1:], sine))
        #pp.plot(chirp, '.')
        #pp.show()
        sine = []

    return chirp

#delays a chirp by n_samples
def delay_a_chirp(chirp, n_samples):

    size_of_chirp = len(chirp)
    delayed_chirp = np.zeros(size_of_chirp)

    #delayed_chirp = np.right_shift(chirp, n_samples)

    for i in range(0, n_samples):
        delayed_chirp[i] = 0

    for i in range(n_samples, size_of_chirp):
        delayed_chirp[i] = chirp[i-n_samples]

    return delayed_chirp


# Generate a sine wave
# A = amplitude
# FS = sampling frequency, in Hertz
# f = frequency of the sine wave, in Hertz
# N_CYCLES = number of cycles of the sine wave
# PHASE =  starting time offset of the sine wave, in degrees
def generate_sine_wave(A, FS, f, N_CYCLES, PHASE):

    TS = 1/FS
    T = 1/f
    n = T/TS
    t = np.round(np.arange(0, N_CYCLES*T, TS),6)
    #t = np.arange(0, 2*np.pi, TS)
    w = 2*np.pi*f
    # w*t is in radians

    y = A*np.sin(w*t + degree_to_rad(PHASE))

    #pp.plot(t, y, ".")
    #pp.show()

    return y, t

def degree_to_rad(degree):
    rad = degree*(np.pi/180)
    return rad

# Append N_ZEROS as header and trailer to a signal
def append_header_and_trailer_to_signal(input_signal, N_ZEROS):
    return np.concatenate((np.zeros(N_ZEROS), input_signal, np.zeros(N_ZEROS)))

# Append N_ZEROS as header to a signal
def append_header_to_signal(input_signal, N_ZEROS):
    return np.concatenate((np.zeros(N_ZEROS), input_signal))

# Append N_ZEROS as trailer to a signal
def append_trailer_to_signal(input_signal, N_ZEROS):
    return np.concatenate((input_signal, np.zeros(N_ZEROS)))

def rad_to_degree(rad):
    degree = rad*(180/np.pi)
    return degree

def polar_to_complex(mag, angle):
    if -2*np.pi <= angle <= 2*np.pi:
        x = mag*np.cos(angle)
        y = mag*np.sin(angle)
    else:
        print("[RF LIB ERR] Angle must be in radains!")
        x = 0
        y = 0
    return np.complex(x, y)

def compute_stability_factor(S):
    deltaS = S[0, 0]*S[1, 1] - S[0, 1]*S[1, 0]
    sqr_mod_dS = np.power(abs(deltaS), 2)

    k = (1 - np.power(abs(S[1, 1]), 2) - np.power(abs(S[0, 0]), 2) + sqr_mod_dS)/(2*abs(S[0, 1]*S[1, 0]))
    return k

def arg_of_complex(c):
    arg = np.arctan(np.imag(c)/np.real(c))
    return arg

def mW_2_dBm(mw):
    return 10*np.log10(mw)

def dBm_2_mW(dBm):
    return 10**(dBm/10)

def compute_optimal_impedances(S):
    deltaS = S[0, 0] * S[1, 1] - S[0, 1] * S[1, 0]
    sqr_mod_dS = np.power(abs(deltaS), 2)

    C2 = S[1,1] - deltaS*np.conjugate(S[0,0])
    B2 = 1 + np.power(abs(S[1, 1]), 2) - np.power(abs(S[0, 0]), 2) - sqr_mod_dS

    a = B2/(2*abs(C2))
    b = np.power(a, 2)
    c = np.sqrt(b - 1)
    d = a - c
    T2OPT = d*(np.cos(arg_of_complex(C2)) - 1j*np.sin(arg_of_complex(C2)))
    ZOUTOPT = 50*(1+T2OPT)/(1-T2OPT)

    e = S[0, 0] + (T2OPT*S[0, 1]*S[1, 0])/(1 - T2OPT*S[1, 1])
    T1OPT = np.conjugate(e)
    ZINOPT = 50*(1+T1OPT)/(1-T1OPT)

    return ZINOPT, ZOUTOPT

#from: http://www.edatop.com/down/hwrf/mprf/MP-RF-20986.pdf
# Outras referências: Rizzi PA (1988) Microwave engineering, passive devices. Prentice Hall, New Jersy (Chap 4)
# Ludwig R, Bretchko P (2000) RF circuit design, theory and applications. Prentice Hall, New Jersy (Chap 8)
# https://link.springer.com/referenceworkentry/10.1007%2F978-981-4560-75-7_133-1
def lumped_match_impedance(ZS, ZL, frequency):
    ZLOAD = ZL
    ZSOURCE = ZS
    f = frequency
    #print("ZLOAD is:", ZLOAD)
    #print("ZSOURCE is:", ZSOURCE)
    #print("frequency is:", f)
    RS = np.real(ZSOURCE)
    XS = np.imag(ZSOURCE)
    RL = np.real(ZLOAD)
    XL = np.imag(ZLOAD)
    w = 2 * np.pi * f

    if abs(ZLOAD) > abs(ZSOURCE):
        #print("ZLOAD is greater than ZSOURCE")
        X = XS + np.sqrt(RS * (RL - RS) + (RS / RL) * XL * XL)
        B = (RS - RL) / (RL * XS + RS * XL - RL * X)

        L = X / w
        C = B / w

        print("From Source to Load:")
        print("ZSOURCE:", ZSOURCE)
        print("Series Inductance: ", round(L * 1E9, 3), "nH")
        print("Parallel Capacitance:" , round(C * 1E12,3), "pF")
        print("ZLOAD:", ZLOAD)

        ZA = 1/(1j*B + 1/ZLOAD)
        #print("ZA:", ZA)
        YA = 1 / ZA
        #print("YA:", YA)
        YL = 1 / ZLOAD
        #print("YL:", YL)
        YL_to_YA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        ZA_to_ZS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        step_size = (np.imag(YA) - np.imag(YL)) / 10
        #print("Step size:", step_size)
        for i in range(0, 11):
            real = np.real(YL)
            imaginario = np.imag(YL) + i*step_size
            #print("CPX[i]", complex(real,imaginario), i)
            YL_to_YA[i] = [complex(real, imaginario)]
        #print(YL_to_YA)

        step_size = (np.imag(ZSOURCE) - np.imag(ZA)) / 10
        #print("Step size:", step_size)
        for i in range(0, 11):
            real = np.real(ZS)
            imaginario = np.imag(ZA) + i*step_size
            ZA_to_ZS[i] = [complex(real, imaginario)]
        #print(ZA_to_ZS)

        pp.figure(figsize=(6, 6))
        ax = pp.subplot(1, 1, 1, projection='smith')
        pp.plot([10, 100], markevery=1)

        #pp.plot(ZLOAD, label="ZLOAD", datatype=SmithAxes.Z_PARAMETER)
        #pp.plot(ZSOURCE, label="ZSOURCE", datatype=SmithAxes.Z_PARAMETER)
        #pp.plot(ZA, label="ZA", datatype=SmithAxes.Z_PARAMETER)
        pp.plot(YL_to_YA, markevery=1, label="YLOAD to YA", equipoints=10, datatype=SmithAxes.Y_PARAMETER)
        pp.plot(ZA_to_ZS, markevery=1, label="ZA to ZSOURCE", equipoints=10, datatype=SmithAxes.Z_PARAMETER)

        leg = pp.legend(loc="lower right", fontsize=10)
        pp.title("Impedance matching for ZLOAD > ZSOURCE", y=-0.01)
        pp.show()

    elif abs(ZLOAD) < abs(ZSOURCE):
        #print("ZLOAD is lower than ZSOURCE")
        X = -XL + np.sqrt(RL * (RS - RL) + (RL / RS) * XS * XS)
        B = (RS - RL) / (RS * XL + RL * XS + RS * X)
        L = X / w
        C = B / w

        print("From Source to Load:")
        print("ZSOURCE:", ZSOURCE)
        print("Parallel Capacitance: ", round(C * 1E12, 3), "pF")
        print("Series Inductance:", round(L*1E9,3), "nH")
        print("ZLOAD:", ZLOAD)

        ZA = (RL + 1j*(XL + X))
        #print("ZA:", ZA)
        YA = 1 / ZA
        #print("YA:", YA)
        YL = 1 / ZLOAD
        #print("YL:", YL)
        #print("ZL:", ZL)
        YS = 1/ZSOURCE
        #print("YS:", YS)
        ZL_to_ZA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        YA_to_YS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        step_size = (np.imag(ZA) - np.imag(ZL)) / 10
        #print("Step size:", step_size)
        for i in range(0, 11):
            real = np.real(ZL)
            imaginario = np.imag(ZL) + i*step_size
            #print("CPX[i]", complex(real,imaginario), i)
            ZL_to_ZA[i] = [complex(real, imaginario)]
        #print(ZL_to_ZA)

        step_size = (np.imag(YS) - np.imag(YA)) / 10
        #print("Step size:", step_size)
        for i in range(0, 11):
            real = np.real(YS)
            imaginario = np.imag(YA) + i*step_size
            YA_to_YS[i] = [complex(real, imaginario)]
        #print(YA_to_YS)

        pp.figure(figsize=(6, 6))
        ax = pp.subplot(1, 1, 1, projection='smith')
        pp.plot([10, 100], markevery=1)

        #pp.plot(ZLOAD, label="ZLOAD", datatype=SmithAxes.Z_PARAMETER)
        #pp.plot(ZSOURCE, label="ZSOURCE", datatype=SmithAxes.Z_PARAMETER)
        #pp.plot(ZA, label="ZA", datatype=SmithAxes.Z_PARAMETER)
        pp.plot(ZL_to_ZA, markevery=1, label="ZLOAD to ZA", equipoints=10, datatype=SmithAxes.Z_PARAMETER)
        pp.plot(YA_to_YS, markevery=1, label="YA to YSOURCE", equipoints=10, datatype=SmithAxes.Y_PARAMETER)

        leg = pp.legend(loc="upper right", fontsize=10)
        pp.title("Impedance matching for ZLOAD < ZSOURCE")

        pp.show()

    else:
        print("[RF LIB WAR] No matching needed!")

def example():


    VDC = 4
    IDC = 250*1E-3
    ROPT = complex(VDC/IDC)
    print(ROPT)

    pp.figure(figsize=(6, 6))
    ax = pp.subplot(1, 1, 1, projection='smith')
    pp.plot([10, 100], markevery=1)

    pp.plot(ROPT, label="ROPT", datatype=SmithAxes.Z_PARAMETER)
    #pp.plot(ZL_to_ZA, markevery=1, label="ZLOAD to ZA", equipoints=10, datatype=SmithAxes.Z_PARAMETER)
    #pp.plot(YA_to_YS, markevery=1, label="YA to YSOURCE", equipoints=10, datatype=SmithAxes.Y_PARAMETER)

    leg = pp.legend(loc="upper right", fontsize=10)
    #pp.title("Impedance matching for ZLOAD < ZSOURCE")

    pp.show()

def add_shunt_capacitor(ZS, C, frequency):
    # Y = G + jB
    # admitância = condutância + j*susceptancia
    # Z = R + JX
    # impedância = resistência +j*reatância

    #0. Define the source admitance YS
    #1. Calculate the capacitive reactance XC
    #2. Define the capactive impedance ZC
    #3. Do the parallel(ZA) between XC and ZS
    #4. Invert ZA into admitance(YA)

    ZSOURCE = ZS
    f = frequency
    w = 2 * np.pi * f

    # 0. Define the source admitance YS
    YS = 1/ZS

    # 1. Calculate the capacitive reactance (XC)
    XC = -1 / (w * C)

    # 2. Define the capactive impedance ZC
    RC = 0
    ZC = RC + 1j*XC

    # 3. Do the parallel(ZA) between XC and ZS
    ZA = ZS * ZC / (ZS + ZC)

    # 4. Invert ZA into admitance(YA)
    YA = 1 / ZA

    # As a susceptance is being added, it should move along the constant susceptance circle
    plot_constant_susceptance(YA, YS)
    return ZA

def add_series_capacitor(ZS, C, frequency):
    # Y = G + jB
    # admitância = condutância + j*susceptancia
    # Z = R + JX
    # impedância = resistência +j*reatância

    #0. Define the source admitance YS
    #1. Calculate the capacitive reactance XC
    #2. Define the capactive impedance ZC
    #3. Do the series (ZA) between XC and ZS

    ZSOURCE = ZS
    f = frequency
    w = 2 * np.pi * f

    # 1. Calculate the capacitive reactance (XC)
    XC = -1 / (w * C)

    # 2. Define the capactive impedance ZC
    RC = 0
    ZC = RC + 1j*XC

    # 3. Do the series(ZA) between XC and ZS
    ZA = ZS + ZC

    # As a reactance is being added, it should move along the constant susceptance circle
    plot_constant_reactance(ZA, ZS)
    return ZA

def add_shunt_inductor(ZS, L, frequency):
    # Y = G + jB
    # admitância = condutância + j*susceptancia
    # Z = R + JX
    # impedância = resistência +j*reatância

    #0. Define the source admitance YS
    #1. Calculate the inductive reactance XL
    #2. Define the inductive impedance ZL
    #3. Do the parallel(ZA) between XL and ZS
    #4. Invert ZA into admitance(YA)

    ZSOURCE = ZS
    f = frequency
    w = 2 * np.pi * f

    # 0. Define the source admitance YS
    YS = 1/ZS

    # 1. Calculate the capacitive reactance (XL)
    XL = w*L

    # 2. Define the inductive impedance ZL
    RL = 0
    ZL = RL + 1j*XL

    # 3. Do the parallel(ZA) between XC and ZS
    ZA = ZS * ZL / (ZS + ZL)

    # 4. Invert ZA into admitance(YA)
    YA = 1 / ZA

    # As a susceptance is being added, it should move along the constant susceptance circle
    plot_constant_susceptance(YA, YS)
    return ZA

def add_series_inductor(ZS, L, frequency):
    # Y = G + jB
    # admitância = condutância + j*susceptancia
    # Z = R + JX
    # impedância = resistência +j*reatância

    #0. Define the source admitance YS
    #1. Calculate the capacitive reactance XC
    #2. Define the capactive impedance ZC
    #3. Do the series (ZA) between XC and ZS

    ZSOURCE = ZS
    f = frequency
    w = 2 * np.pi * f

    # 1. Calculate the inductive reactance (XL)
    XL = w*L

    # 2. Define the capactive impedance ZC
    RL = 0
    ZL = RL + 1j*XL

    # 3. Do the series(ZA) between XC and ZS
    ZA = ZS + ZL

    # As a reactance is being added, it should move along the constant susceptance circle
    plot_constant_reactance(ZA, ZS)
    return ZA

def plot_constant_reactance(ZA, ZS):
    ZL_to_ZA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    step_size = (np.imag(ZA) - np.imag(ZS)) / 10
    # print("Step size:", step_size)
    for i in range(0, 11):
        real = np.real(ZS)
        imaginario = np.imag(ZS) + i * step_size
        ZL_to_ZA[i] = [complex(real, imaginario)]

    pp.figure(figsize=(6, 6))
    ax = pp.subplot(1, 1, 1, projection='smith')
    pp.plot([10, 100], markevery=1)

    pp.plot(ZL_to_ZA, markevery=1, label="ZLOAD to ZA", equipoints=10, datatype=SmithAxes.Z_PARAMETER)

    leg = pp.legend(loc="upper right", fontsize=10)
    pp.title("ZS=" + str(ZS) + "\n" + "ZA=" + str(ZA), y=-0.01)
    pp.show()

def plot_constant_susceptance(YA, YS):

    YS_to_YA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    step_size = (np.imag(YA) - np.imag(YS)) / 10
    for i in range(0, 11):
        real = np.real(YS)
        imaginario = np.imag(YS) + i * step_size
        YS_to_YA[i] = [complex(real, imaginario)]

    pp.figure(figsize=(6, 6))
    ax = pp.subplot(1, 1, 1, projection='smith')
    pp.plot([10, 100], markevery=1)

    pp.plot(YS_to_YA, markevery=1, label="YSOURCE to YA", equipoints=10, datatype=SmithAxes.Y_PARAMETER)

    #leg = pp.legend(loc="lower right", fontsize=10)
    pp.title("ZS="+str(1/YS)+"\n"+"ZA="+str(1/YA), y=-0.01)
    pp.show()


