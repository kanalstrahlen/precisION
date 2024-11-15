# fastDatabaseSearch.pyx
import numpy as np
cimport numpy as np
from libc.math cimport fabs
from tqdm import tqdm


cdef class FastDbSearchFunctions:

    cdef float threshold

    def __init__(self, threshold):
        self.threshold = threshold

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
                theo_pointer += 1
            elif exp_mass < theo_mass:
                exp_pointer += 1
            else:
                theo_pointer += 1

        return matches


    cpdef count_matches_for_proteins(self, theoretical_ions_for_each_proteoform, observed_ions):
        observed_ions.sort()
        matches_per_proteoform = []

        # iterate through each proteoform and count matches
        for theoretical_ions in tqdm(theoretical_ions_for_each_proteoform, total=len(theoretical_ions_for_each_proteoform)):
            theoretical_ions.sort()
            theoretical_ions = [ion for ion in theoretical_ions if ion > 0]
            match_count = self.count_matches(observed_ions, theoretical_ions)
            matches_per_proteoform.append(match_count)

        return matches_per_proteoform
