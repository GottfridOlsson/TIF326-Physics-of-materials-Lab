##====================================================##
##     Project: TIF326 LAB HEAT TREATMENT OF Al-7075
##        File: formatting-raw-CSV.py
##      Author: GOTTFRID OLSSON 
##     Created: 2023-02-14
##     Updated: 2023-02-16
##       About: Formatting raw CSV from Canvas and
##              calculating sigma_UTS, sigma_0.2
##====================================================##

import matplotlib.pyplot as plt
import numpy as np
import util

# Paths
path_raw_folder = "raw-CSV/"
paths_start = ["Group1-130C", "Group2-160C", "Group3-190C", "Group5-250C"] #group 4 is missing from canvas, 2023-02-14
path_raw_end = "for2h.is_tens.raw.csv"

path_formatted_folder = "formatted-CSV/"
path_formatted_end = "-2h_tensile-test.csv"

# Measured cross-section area of samples
area_mm2 = [12.25*1.94, 12.06*1.93, 12.08*1.94, 12.19*1.94] #group 1,2,3,5; measurements in (mm)^2
area_m2 = np.array(area_mm2) * 1e-6

def get_index_of_first_instance_where_element_in_list_is_larger_than(list, value):
    index = 0
    for i,element in enumerate(list):
        if element > value:
            break
        index += 1
    return index

for i, path_start in enumerate(paths_start):

    # Paths
    path_raw       = path_raw_folder + path_start + path_raw_end
    path_formatted = path_formatted_folder + path_start + path_formatted_end

    # Read data
    data = np.genfromtxt(path_raw, delimiter=';', skip_header=1)
    time_second    = data[:,0]
    extension_mm   = data[:,1]
    strain_percent = data[:,2]
    load_newton    = data[:,3]

    # Calculate sigma
    stress_Pascal  = load_newton/area_m2[i]

    # Print formatted data
    util.print_arrays_to_CSV(f"{path_formatted_folder}{paths_start[i]}{path_formatted_end}", 
                                "Strain (percent)", strain_percent, 
                                "Stress (MPa)", stress_Pascal*1e-6, 
                                "Extension (mm)", extension_mm, 
                                "Time (s)", time_second, 
                                print_message=True)

    # Calculate specific values
    #if i == 1:
    sigma_UTS_Pascal = np.max(stress_Pascal)
    index_sigma_point_two_percent = get_index_of_first_instance_where_element_in_list_is_larger_than(strain_percent, 0.2)
    sigma_point_two_percent_Pascal = stress_Pascal[index_sigma_point_two_percent]
    fracture_strain_percent = np.max(strain_percent)
    print(f"Group 2: sigma_UTS={sigma_UTS_Pascal*1e-6:.1f} MPa, sigma_0.2%={sigma_point_two_percent_Pascal*1e-6:.1f} MPa, fracture_strain={fracture_strain_percent:.2f} %")


    plt.plot(strain_percent, stress_Pascal*1e-6, linestyle='-', marker='.')
    plt.xlabel("Strain (%)")
    plt.ylabel("Stress (MPa)")
    plt.show()
