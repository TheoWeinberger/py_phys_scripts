#!/usr/bin/env python3
import re
import sys
import ast
"""
returns dict of band energy
to irp k
"""
def read_energies(case, band_index):

    filename = case + '.energyup'
    
    f = open(filename, "r")
    lines=f.readlines()
    energy_dict = {}
    irp = 1
    for line in lines:

        if ' '+str(band_index)+'  ' in line:
            line = line.split('  ')
            line = [val for val in line if val != '']
            if len(line) == 3:
                energy_dict[irp] = float(line[1])
                irp += 1

    return energy_dict



"""
get fermi energy
"""
def get_fermi(case):

    filename = case + '.in1'
    
    f = open(filename, "r")
    lines=f.readlines()
    
    """
    this corresponds to code for gettinf the fermi eneryg
    """
    # get fermi energy
    x = lines[0].split(" ")
    x = x[2].split("=")

    # convert ef to a number
    e_f = float(x[-1])

    # get fermi energy

    return e_f

"""
get dict of k to irrep k value
"""
def read_k(case):
      
    filename = case + '.outputkgen'
    f = open(filename, "r")
    lines=f.readlines()

    counter = 1
    for line in lines:
        if "DIVISION OF RECIPROCAL LATTICE VECTORS (INTERVALS)" in line:
            # get dimensions
            x = re.findall(r'\d+', line)
            k_dict = {}
            #convert dimesions to integer values
            dimensions = [int(dim) + 1 for dim in x]

        #find line when k mesh begins
        if "point    coordinates     relation" in line:
            break

        counter += 1

    all_reference = []
    for i in range(counter, counter+dimensions[0]*dimensions[1]*dimensions[2]):
        x = re.findall(r'\d+', lines[i])
        all_reference += [x[4]]
        k_dict[x[0]] = x[4]

   
    reference_to_irp = {}
    ref = 1
    for item in list(dict.fromkeys(all_reference)):
        reference_to_irp[item] = ref
        ref += 1


    return k_dict, reference_to_irp, dimensions



def gen_bxsf(case, bands):

    e_f = get_fermi(case)

    k_dict, reference_to_irp, dimensions = read_k(case)
        
    f2 = open(case +".bxsf.bands", "w")

    #do it again for combined one
    f2.write(" BEGIN_INFO\n")
    f2.write("  Fermi Energy:     " + str(e_f) + "\n")
    f2.write(" END_INFO\n")
    f2.write(" BEGIN_BLOCK_BANDGRID_3D\n")
    f2.write(" band_energies\n")
    f2.write(" BANDGRID_3D_BANDS\n")
    f2.write(" " + str(len(bands)) + "\n")
    f2.write(" " + str(dimensions[0]) + " " + str(dimensions[1]) + " " + str(dimensions[2]) +"\n")

    """
    get some k space details
    """
    filename = case + '.outputkgen'
    k_file = open(filename, "r")
    lines = k_file.readlines()

    for line in lines:
        if "SHIFT" in line:
            shift = re.findall(r'\d+', line)
            shift = [float(val) for val in shift]
            f2.write("      " + str(shift[0]) + "0000000      " + str(shift[1]) + "0000000      " + str(shift[2]) +"0000000\n")
            break
    
    counter = 1
    for line in lines:
        if "G1        G2        G3" in line:
            break
        counter += 1
    

            

    #get basis vectors
    vec1 = lines[counter].split(" ")
    vec1 = [val.strip("\n") for val in vec1]
    vec1 = [val.strip(" ") for val in vec1]
    vec1 = [val for val in vec1 if val != ""]
    #vec1 = vec1[1:]

    vec2 = lines[counter+1].split(" ")
    vec2 = [val.strip("\n") for val in vec2]
    vec2 = [val.strip(" ") for val in vec2]
    vec2 = [val for val in vec2 if val != ""]
    #vec2 = vec2[1:]

    vec3 = lines[counter+2].split(" ")
    vec3 = [val.strip("\n") for val in vec3]
    vec3 = [val.strip(" ") for val in vec3]
    vec3 = [val for val in vec3 if val != ""]
    #vec3 = vec3[1:]

    vec_list = [vec1, vec2, vec3]

    for vec in vec_list:

        counter = 0

        for val in vec:
            if counter == 0:
                
                if float(val) >= 0:
                    f2.write("      " + str(val))
                else:
                    f2.write("     " + str(val))
            
            else:
                if float(val) >= 0:
                    f2.write("00      " + str(val))
                else:
                    f2.write("00     " + str(val))

                
            if counter == 2:

                f2.write("00\n")

            counter += 1


    for band in bands:
        f = open(case +".bxsf.band-" + str(band), "w")

        """
        general .bxsf rubbish
        """
        f.write(" BEGIN_INFO\n")
        f.write("  Fermi Energy:     " + str(e_f) + "\n")
        f.write(" END_INFO\n")
        f.write(" BEGIN_BLOCK_BANDGRID_3D\n")
        f.write(" band_energies\n")
        f.write(" BANDGRID_3D_BANDS\n")
        f.write(" 1\n")
        f.write(" " + str(dimensions[0]) + " " + str(dimensions[1]) + " " + str(dimensions[2]) +"\n")

        """
        get some k space details
        """
        filename = case + '.outputkgen'
        k_file = open(filename, "r")
        lines = k_file.readlines()

        for line in lines:
            if "SHIFT" in line:
                shift = re.findall(r'\d+', line)
                shift = [float(val) for val in shift]
                f.write("      " + str(shift[0]) + "0000000      " + str(shift[1]) + "0000000      " + str(shift[2]) +"0000000\n")
                break

        #get basis vectors
        vec1 = lines[counter].split(" ")
        vec1 = [val.strip("\n") for val in vec1]
        vec1 = [val.strip(" ") for val in vec1]
        vec1 = [val for val in vec1 if val != ""]
        #vec1 = vec1[1:]

        vec2 = lines[counter+1].split(" ")
        vec2 = [val.strip("\n") for val in vec2]
        vec2 = [val.strip(" ") for val in vec2]
        vec2 = [val for val in vec2 if val != ""]
        #vec2 = vec2[1:]

        vec3 = lines[counter+2].split(" ")
        vec3 = [val.strip("\n") for val in vec3]
        vec3 = [val.strip(" ") for val in vec3]
        vec3 = [val for val in vec3 if val != ""]
        #vec3 = vec3[1:]


        for vec in vec_list:

            counter = 0

            for val in vec:
                if counter == 0:
                    
                    if float(val) >= 0:
                        f.write("      " + str(val))
                    else:
                        f.write("     " + str(val))
                
                else:
                    if float(val) >= 0:
                        f.write("00      " + str(val))
                    else:
                        f.write("00     " + str(val))

                    
                if counter == 2:

                    f.write("00\n")

                counter += 1
                

        f.write(" BAND:   " + str(band) + "\n")


        energies_unfold = []

        energies = read_energies(case, band)

        for key in k_dict.keys():

            energies_unfold += [energies.get(int(reference_to_irp.get(k_dict.get(key))))]
            #if energies.get(int(k_dict.get(key))) is None:
                #print(k_dict.get(key))

        for i in range(dimensions[0]*dimensions[1]*dimensions[2]):
           
            f.write(" " + str(energies_unfold[i]))
            if (i+1)%6 == 0 :
                 f.write("\n")
        
        #last new line
        if (i+1)%6 == 0 :
            pass
        else:
            f.write("\n")
        
        f.write(" END_BANDGRID_3D\n")
        f.write(" END_BLOCK_BANDGRID_3D\n")
        f.close()
        
        f2.write(" BAND:   " + str(band) + "\n")


        energies_unfold = []

        energies = read_energies(case, band)

        for key in k_dict.keys():

            energies_unfold += [energies.get(int(reference_to_irp.get(k_dict.get(key))))]
            #if energies.get(int(k_dict.get(key))) is None:
                #print(k_dict.get(key))

        for i in range(dimensions[0]*dimensions[1]*dimensions[2]):
           
            f2.write(" " + str(energies_unfold[i]))
            if (i+1)%6 == 0 :
                 f2.write("\n")

        if (i+1)%6 == 0 :
            pass
        else:
            f2.write("\n")
        
    f2.write(" END_BANDGRID_3D\n")
    f2.write(" END_BLOCK_BANDGRID_3D\n")
    f2.close()


        

            
        
        
        
"""
Input case name and desired bands here

requires case.outputkgen, case.in1, case.energyup

therefore to look at non sp cases rename case.energy to case.energyup
for down electrons rename case.energydn to case.energyup
for so rename case.energyso to case.energyup etc...
"""
"""
Input case name and desired bands here

requires case.outputkgen, case.in1, case.energyup

therefore to look at non sp cases rename case.energy to case.energyup
for down electrons rename case.energydn to case.energyup
for so rename case.energyso to case.energyup etc...
"""
#get base file name
file_name = sys.argv[1]
band_range = ast.literal_eval(sys.argv[2])
#
gen_bxsf(file_name, band_range)





