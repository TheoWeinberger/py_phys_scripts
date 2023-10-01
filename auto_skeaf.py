#!/usr/bin/env python3
import sys
import glob
import subprocess
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import re
import os
import numpy as np
import scienceplots

plt.style.use(['science'])

cp = sns.color_palette("bright")

def insert(originalfile,string):
    with open(originalfile,'r') as f:
        with open('config_temp.in','w') as f2: 
            f2.write(str(string))
            f2.write(f.read())
    f2.close()

def insert2(originalfile,string):
    with open(originalfile,'r') as f:
        with open('config.in','w') as f2: 
            f2.write(str(string) + "\n")
            f2.write(f.read())
    f2.close()

#get base file name
file_name = sys.argv[1]

#get Ef max shift
min_shift = float(sys.argv[2])
max_shift = float(sys.argv[3])

#get Ef shift intervals
shift_intervals = float(sys.argv[4])

files = glob.glob(file_name + ".bxsf.band-*")

if shift_intervals > 0:

    for Ef_shift in np.arange(min_shift, max_shift, shift_intervals):

        Ef = 0

        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

        for file in files:
            #generate config file for ca
            file_io = open(file, "r")

            #get Fermi energy
            for line in file_io:
                if re.search("Fermi Energy:", line):
                    line = line.replace(" ", "")
                    Ef = line.removeprefix("FermiEnergy:")
                    break

            #Screen for Ef crossing

            data_vals = pd.read_csv(file, skiprows=13, skipfooter=2, header=None, delimiter=" ", engine='python')
            
            #get max value 
            max_energy = data_vals.max().max()

            #get min energy 
            min_energy = data_vals.min().min()

            Ef_float = float(Ef)

            Ef_float_shifted = Ef_float + Ef_shift 

            if min_energy <= Ef_float_shifted and max_energy >= Ef_float_shifted:

                band_name = file.split(".bxsf.")[-1]

                #add file name and fermi energy to config
                insert("config_ca.in", str(Ef_float_shifted) + "\n")
                insert2("config_temp.in", file)

                popen = subprocess.Popen([r"skeaf", "-rdcfg", "-nodos"])

                popen.wait()

                data = pd.read_csv("results_freqvsangle.out")

                os.rename("results_freqvsangle.out", file_name + "_" + str(Ef_float_shifted) + "." + band_name + ".c-a.out")

                ax1.scatter(data["Phi(deg)"], data["Freq(kT)"], label = band_name, s = 10)
                ax1.set_ylabel("Frequency (kT)")
                ax1.set_xlabel("c-a angle")
                #ax1.legend()

                insert("config_cb.in", str(Ef_float_shifted) + "\n")
                insert2("config_temp.in", file)

                #generate config file for ca
                popen = subprocess.Popen([r"skeaf", "-rdcfg", "-nodos"])

                data = pd.read_csv("results_freqvsangle.out")

                os.rename("results_freqvsangle.out", file_name + "_" + str(Ef_float_shifted) + "." + band_name + ".c-b.out")

                popen.wait()

                ax2.scatter(data["Phi(deg)"], data["Freq(kT)"], label = band_name, s = 10)
                ax2.set_xlabel("c-b angle")
                ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.savefig(file_name + "_" + str(Ef_float_shifted) + ".png", format='png', dpi=1200)

else:

    Ef = 0

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    for file in files:
        #generate config file for ca
        file_io = open(file, "r")

        #get Fermi energy
        for line in file_io:
            if re.search("Fermi Energy:", line):
                line = line.replace(" ", "")
                Ef = line.removeprefix("FermiEnergy:")
                break

        #Screen for Ef crossing

        data_vals = pd.read_csv(file, skiprows=13, skipfooter=2, header=None, delimiter=" ", engine='python')
        
        #get max value 
        max_energy = data_vals.max().max()

        #get min energy 
        min_energy = data_vals.min().min()

        Ef_float = float(Ef)

        #get band name 

        if min_energy <= Ef_float and max_energy >= Ef_float:

            band_name = file.split(".bxsf.")[-1]

            #add file name and fermi energy to config
            insert("config_ca.in", Ef)
            insert2("config_temp.in", file)

            popen = subprocess.Popen([r"skeaf", "-rdcfg", "-nodos"])

            popen.wait()

            data = pd.read_csv("results_freqvsangle.out")

            os.rename("results_freqvsangle.out", file_name + "_" + Ef.strip() + "." + band_name + ".c-a.out")

            ax1.scatter(data["Phi(deg)"], data["Freq(kT)"], label = band_name, s = 10)
            ax1.set_ylabel("Frequency (kT)")
            ax1.set_xlabel("c-a angle")
            #ax1.legend()

            insert("config_cb.in", Ef)
            insert2("config_temp.in", file)

            #generate config file for ca
            popen = subprocess.Popen([r"skeaf", "-rdcfg", "-nodos"])

            popen.wait()

            data = pd.read_csv("results_freqvsangle.out")

            os.rename("results_freqvsangle.out", file_name + "_" + Ef.strip() + "." + band_name + ".c-b.out")

            ax2.scatter(data["Phi(deg)"], data["Freq(kT)"], label = band_name, s = 10)
            ax2.set_xlabel("c-b angle")
            ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))


    plt.savefig(file_name + "_" + Ef.strip() + ".png", format='png', dpi=1200)
