# fastOffsetSearch.pyx
# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
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

        # Convert inputs to contiguous C buffers once, instead of rebuilding a Python
        # list of shifted ions on every offset. The experimental peaks are sorted so the
        # two-pointer merge below is correct regardless of input order; theo_ions are
        # already ascending (b-ions, and y-ions are reversed to ascending upstream).
        cdef double[::1] exp = np.sort(np.asarray(peak_list, dtype=np.float64))
        cdef double[::1] base = np.ascontiguousarray(theo_ions, dtype=np.float64)
        cdef np.ndarray[np.int64_t, ndim=1] counts = np.empty(num_points, dtype=np.int64)
        cdef np.int64_t[::1] cv = counts

        cdef int n = exp.shape[0]
        cdef int m = base.shape[0]
        cdef bint theo_mode = (self.count_mode == "Theoretical ions")

        cdef Py_ssize_t k = 0
        cdef int i, j, jstart, matches
        cdef double off
        cdef float exp_mass
        cdef float theo_mass
        cdef float ppm_difference

        for off_val in tqdm(offset_list, total=len(offset_list)):
            off = off_val

            # Equivalent to filtering out ions where (ion_mass + offset) < 0: because
            # base is ascending these are always a leading prefix, so skip past them.
            jstart = 0
            while jstart < m and (base[jstart] + off) < 0.0:
                jstart += 1

            i = 0
            j = jstart
            matches = 0
            while i < n and j < m:
                exp_mass = <float>exp[i]
                theo_mass = <float>(base[j] + off)
                ppm_difference = <float>(fabs(exp_mass - theo_mass) / theo_mass * 1e6)

                if ppm_difference < self.threshold:
                    matches += 1
                    i += 1
                    if theo_mode:
                        j += 1
                elif exp_mass < theo_mass:
                    i += 1
                else:
                    j += 1

            cv[k] = matches
            k += 1

        return offset_list, counts