# fastOffsetSearch.pyx
import numpy as np
cimport numpy as np
from libc.math cimport fabs
import bisect
from tqdm import tqdm

cdef class FastOffsetShiftCalculator:

    cdef float threshold
    cdef str count_mode
    cdef float lower_range
    cdef float upper_range

    def __init__(self, threshold, count_mode, lower_range, upper_range):
        self.threshold = threshold
        self.count_mode = count_mode
        self.lower_range = lower_range
        self.upper_range = upper_range

    cpdef count_matches(self, exp_ions, theo_ions):
        cdef int exp_pointer = 0
        cdef int theo_pointer = 0
        cdef int matches = 0
        cdef float exp_mass
        cdef float theo_mass
        cdef float ppm_difference

        while exp_pointer < len(exp_ions) and theo_pointer < len(theo_ions):
            exp_mass = exp_ions[exp_pointer]
            theo_mass = theo_ions[theo_pointer]
            ppm_difference = fabs(exp_mass - theo_mass) / theo_mass * 1e6

            if ppm_difference < self.threshold:
                matches += 1
                exp_pointer += 1
                if self.count_mode == "Theoretical ions":
                    theo_pointer += 1
                elif self.count_mode == "Observed ions":
                    pass  # No change in theo_pointer
            elif exp_mass < theo_mass:
                exp_pointer += 1
            else:
                theo_pointer += 1

        return matches

    cpdef count_matched_ions(self, peak_list, theo_ions):
        cdef int match_count = 0
        cdef int num_matches
        cdef int ion
        cdef int low, high
        cdef float ppm_threshold = self.threshold

        for ion in theo_ions:
            low = bisect.bisect_left(peak_list, ion * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(peak_list, ion * (1 + ppm_threshold / 1e6))
            num_matches = high - low

            if num_matches >= 1:
                if self.count_mode == "Theoretical ions":
                    match_count += 1
                elif self.count_mode == "Observed ions":
                    match_count += num_matches

        return match_count

    cpdef offset_scan(self, peak_list, theo_ions):
        cdef float scan_spacing = self.threshold / 1000 / 1.5 # 2/3 the threshold at m/z 1000
        cdef int num_points = int((self.upper_range - self.lower_range) // scan_spacing)
        cdef np.ndarray offset_list = np.linspace(self.lower_range, self.upper_range, num=num_points)
        cdef list match_count_list = []

        for offset in tqdm(offset_list, total=len(offset_list)):
            temp_theo_ions = [ion_mass + offset for ion_mass in theo_ions]
            temp_theo_ions = [ion_mass for ion_mass in temp_theo_ions if ion_mass >= 0]
            match_count = self.count_matches(peak_list, temp_theo_ions)
            match_count_list.append(match_count)

        return offset_list, match_count_list