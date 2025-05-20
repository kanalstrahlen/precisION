import subprocess
import os
import warnings
import numpy as np
from psims.mzml import MzMLWriter
from pyteomics import mzml
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import convolve
from tqdm import tqdm
import plotly.graph_objects as go

warnings.filterwarnings('ignore')

class RLPeakPicking():
    def __init__(self, file_path, directory_path, html_file_save):
        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(self.file_path).replace(".txt", "")
        self.output_profile = os.path.join(self.directory_path, f"{self.basename}_profile.mzML")
        self.output_centroid = os.path.join(self.directory_path, f"{self.basename}_centroid.mzML")
        self.output_plot = os.path.join(self.directory_path, f"{self.basename}.resolutionPlot.pdf")
        self.output_html = os.path.join(self.directory_path, f"{self.basename}.deconvolution.html")

        cwt_peaks = self.cwt_centroid()
        noise_coefficients = self.noise_level_calculation(cwt_peaks) # coefficients for quadratic fit of noise level
        interpolated_mz, interpolated_intensity = self.interpolate_spectrum()
        res_coefficients = self.resolution_estimate(interpolated_mz, interpolated_intensity)
        rl_peaks, convolved_spectrum = self.richardson_lucy_centroid(
            interpolated_mz,
            interpolated_intensity,
            res_coefficients,
            noise_coefficients,
            html_file_save
        )

        mz = [item[0] for item in rl_peaks]
        intensity = [item[1] for item in rl_peaks]
        self.make_mzml(mz, intensity, 'centroid')

        self.save_data(cwt_peaks, rl_peaks, convolved_spectrum)


    def save_data(self, cwt_peaks, rl_peaks, convolved_spectrum):
        output_cwt = os.path.join(self.directory_path, f"{self.basename}.centroid.cwt.txt")
        output_rl = os.path.join(self.directory_path, f"{self.basename}.centroid.rl.txt")
        output_convolved = os.path.join(self.directory_path, f"{self.basename}.profile.rl.txt")

        with open(output_cwt, "w") as file:
            for row in cwt_peaks:
                file.write(f"{row[0]}\t{row[1]}\n")

        with open(output_rl, "w") as file:
            for row in rl_peaks:
                file.write(f"{row[0]}\t{row[1]}\n")

        with open(output_convolved, "w") as file:
            for row in convolved_spectrum:
                file.write(f"{row[0]}\t{row[1]}\n")


    def richardson_lucy_centroid(self, mz, intensity, res, noise, html_file_save):
        num_windows = int((max(mz) - min(mz))/100)
        window_width = int(len(mz) / num_windows)

        peaks_mz = []
        peaks_intensity = []
        total_convolution = []
        total_mz = []

        print("Running RL deconvolution.")
        for i in tqdm(range(num_windows)):
            window_min = i * window_width
            window_max = (i+1) * window_width

            mz_values = mz[window_min:window_max]
            observed_spectrum = intensity[window_min:window_max]
            num_points = len(mz_values)
            nonzero_indices = np.where(
                np.logical_and(
                    np.roll(observed_spectrum, 1) > 1e-10,
                    np.roll(observed_spectrum, -1) > 1e-10
                )
            )[0]
            observed_spectrum[~np.isin(np.arange(len(observed_spectrum)), nonzero_indices)] = 1e-10

            # speed is primarily affected by spacing, i have reduced this to a feasible level now w/ interpolation

            # Initialize parameters
            num_iterations = 100000
            convergence_threshold = 1e-4 # low convergence for intial solution as much harder to converge on true solution

            estimated_signal = np.copy(observed_spectrum)
            prev_estimated_signal = np.copy(estimated_signal)

            fac = 1.6

            psf_sigma = (res[0] * ((np.min(mz_values) - res[2]) ** res[1])) + res[3]


            for iteration in range(num_iterations):
                prev_estimated_signal = np.copy(estimated_signal)
                gaussian_psf = np.exp(-(mz_values - mz_values[num_points // 2])**2 / (2 * fac * psf_sigma **2))
                convolved_estimate = convolve(estimated_signal, gaussian_psf, mode='same')
                ratio = observed_spectrum / convolved_estimate
                estimated_signal *= ratio
                squared_diff_sum = np.sum(np.diff(estimated_signal - prev_estimated_signal)**2)
                totalsum = np.sum(estimated_signal)
                change = squared_diff_sum/totalsum
                if change < convergence_threshold:
                    break

            window_size = 2.355 * psf_sigma
            window_size = int(window_size  / (mz_values[1] - mz_values[0]))
            clustered_estimate = np.zeros_like(estimated_signal)
            clustered_estimate[:] = 1e-10
            sorted_indices = np.argsort(estimated_signal)[::-1]
            clustered_indices = np.zeros_like(estimated_signal, dtype=bool)
            temp_estimate = estimated_signal.copy()
            for sorted_index in sorted_indices:
                if not clustered_indices[sorted_index]:
                    mz_start_idx = max(0, sorted_index - window_size // 2)
                    mz_end_idx = min(len(estimated_signal)-1, sorted_index + window_size // 2 + 1)
                    com_index = sorted_index
                    total_intensity = sum(temp_estimate[mz_start_idx:mz_end_idx])
                    temp_estimate[mz_start_idx:mz_end_idx] = 1e-10
                    clustered_estimate[com_index] = total_intensity
                    clustered_indices[mz_start_idx:mz_end_idx] = True

            peak_indices = np.where(clustered_estimate >= 1e-2)[0]
            estimated_signal = np.copy(clustered_estimate)

            convergence_threshold = 1e-10
            for iteration in range(num_iterations):
                prev_estimated_signal = np.copy(estimated_signal)
                gaussian_psf = np.exp(-(mz_values - mz_values[num_points // 2])**2 / (2 * psf_sigma **2))
                convolved_estimate = convolve(estimated_signal, gaussian_psf, mode='same')
                ratio = observed_spectrum / convolved_estimate
                estimated_signal *= ratio
                estimated_signal[np.isin(np.arange(len(estimated_signal)), peak_indices, invert=True)] = 1e-10
                squared_diff_sum = np.sum(np.diff(estimated_signal - prev_estimated_signal)**2)
                totalsum = np.sum(estimated_signal)
                change = squared_diff_sum/totalsum
                if change < convergence_threshold:
                    break

            window_size = 2.355 * psf_sigma
            window_size = int(window_size  / (mz_values[1] - mz_values[0]))
            clustered_estimate = np.zeros_like(estimated_signal)
            clustered_estimate[:] = 1e-10
            sorted_indices = np.argsort(estimated_signal)[::-1]
            clustered_indices = np.zeros_like(estimated_signal, dtype=bool)
            temp_estimate = estimated_signal.copy()
            for sorted_index in sorted_indices:
                if not clustered_indices[sorted_index]:
                    mz_start_idx = max(0, sorted_index - window_size // 2)
                    mz_end_idx = min(len(estimated_signal)-1, sorted_index + window_size // 2 + 1)
                    com_index = sorted_index
                    total_intensity = sum(temp_estimate[mz_start_idx:mz_end_idx])
                    temp_estimate[mz_start_idx:mz_end_idx] = 1e-10
                    clustered_estimate[com_index] = total_intensity
                    clustered_indices[mz_start_idx:mz_end_idx] = True

            peak_indices = np.where(clustered_estimate >= 1e-2)[0]
            noise_level = noise[0] * mz_values[0] ** 2 + noise[1] * mz_values[0] + noise[2]

            for i in peak_indices:
                if clustered_estimate[i] >= (0.1 * noise_level) :
                    peaks_mz.append(mz_values[i+1])
                    peaks_intensity.append(clustered_estimate[i])

            clustered_estimate = [value if value >= 0.1 * noise_level else 1e-10 for value in clustered_estimate]

            convolved = convolve(clustered_estimate, gaussian_psf, mode='same')
            total_mz.extend(mz_values)
            total_convolution.extend(convolved)

        if html_file_save:
            self.generate_decon_plot(peaks_mz, peaks_intensity, total_mz, total_convolution, mz, intensity)

        return list(zip(peaks_mz, peaks_intensity)), list(zip(total_mz, total_convolution))


    def generate_decon_plot(self, rl_cent_x, rl_cent_y, rl_x, rl_y, mz, intensity):
        peaks = go.Scatter(
            x=rl_cent_x,
            y=rl_cent_y,
            name='Peaks',
            mode='markers',
            visible=True,
            marker=dict(color='blue', line=dict(color='white', width=1))
        )

        line_spectrum = go.Scatter(
            x=mz, y=intensity,
            mode='lines',
            visible=True,
            name='Spectrum',
            line=dict(color='black')
        )
        line_rl_profile = go.Scatter(
            x=rl_x,
            y=rl_y,
            mode='lines',
            visible=True,
            name='Modelled spectrum',
            line=dict(color='red')
        )

        all_traces = [line_spectrum, line_rl_profile, peaks]

        max_y_value = max(rl_y)

        # Create layout
        layout = go.Layout(
            title='RL deconvolution results',
            showlegend=True,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            yaxis=dict(range=[0, 1.1 * max_y_value])
        )

        fig = go.Figure(data=all_traces, layout=layout)

        button_labels = ['Show all', 'Show profile', 'Show raw data']
        button_list = []

        button_list.append(
            dict(
                label='Show all',
                method='update',
                args=[{'visible': [True, True, True]}]
            )
        )

        button_list.append(
            dict(
                label='Show profile',
                method='update',
                args=[{'visible': [True, True, False]}]
            )
        )

        button_list.append(
            dict(
                label='Show peaks',
                method='update',
                args=[{'visible': [True, False, True]}]
            )
        )

        button_list.append(
            dict(
                label='Show raw spectrum',
                method='update',
                args=[{'visible': [True, False, False]}]
            )
        )

        fig.update_layout(
            updatemenus=[
                dict(
                    type='buttons',
                    showactive=True,
                    buttons=button_list,
                    x=1,
                    xanchor='left',
                    y=0.1,
                    yanchor='top'
                ),
            ]
        )

        fig.write_html(self.output_html)


    def interpolate_spectrum(self):
        try:
            mz, intensity = self.read_scan_data(self.file_path)
        except:
            mz, intensity = self.read_scan_data_space(self.file_path)

        interpolation_function = interp1d(mz, intensity)

        #spacing = self.calc_min_spacing(mz)/1.1 #this was too slow in some cases
        #replaced with 0.001 which is a good comprimise in most cases using suggested
        #settings
        # if there is an issue with rl, this is probably it
        spacing = 0.001
        num_mz_points = int((max(mz) - min(mz)) / spacing)

        mz = np.linspace(min(mz), max(mz), num=num_mz_points)
        intensity = interpolation_function(mz)

        return mz, intensity


    def resolution_estimate(self, mz, intensity):
        max_mz = max(mz)
        min_mz = min(mz)
        self.min_mz = min_mz

        spacing = mz[1] - mz[0]

        resolution200 = 240000
        resolution_max = resolution200 * (200 / min_mz) ** 0.5
        min_fwhm = min_mz / resolution_max

        slices = int((max_mz - min_mz) / 40)
        peak_index_list = []
        for i in range(slices):
            peak_index = np.argmax(intensity[int(i * len(mz) / slices):int((i+1) * len(mz) / slices) - 1]) + int(i * len(mz) / slices) - 1
            peak_index_list.append(peak_index)

        peak_param_list = []
        for peak_index in peak_index_list:
            params_peak = self.fit_gaussian(
                mz[peak_index - int(1/spacing):peak_index + int(1/spacing)],
                intensity[peak_index - int(1/spacing):peak_index + int(1/spacing)],
                min_fwhm
            )
            peak_param_list.append(params_peak)

        mz_res_fit = []
        sigma_res_fit = []

        for index, params in enumerate(peak_param_list):
            mz_fit = params[1]
            sigma_fit = params[2]
            if len(sigma_res_fit) >= 1:
                if np.logical_or(
                    sigma_fit > 2 * max(sigma_res_fit),
                    sigma_fit == 0
                    ):
                    pass
                else:
                    mz_res_fit.append(mz_fit)
                    sigma_res_fit.append(sigma_fit)
            else:
                mz_res_fit.append(mz_fit)
                sigma_res_fit.append(sigma_fit)

        plt.scatter(mz_res_fit, sigma_res_fit)
        plt.savefig(self.output_plot)
        plt.clf()
        plt.cla()


        a, b, c, d = self.fit_fwhm_curve(mz_res_fit, sigma_res_fit)


        def power_law(x, a, b, c, d):
            return (a * ((x-c) ** b)) + d

        # generate y values from the fitted curve
        y_fit = power_law(mz_res_fit, a, b, c, d)

        # [lot the original data and the fitted curve
        plt.scatter(mz_res_fit, sigma_res_fit, color='orange')
        plt.plot(mz_res_fit, y_fit, color='black')
        plt.xlabel('m/z')
        plt.ylabel('sigma')
        plt.savefig(self.output_plot)

        return [a, b, c, d]


    def fit_fwhm_curve(self, mz_res_fit, sigma_res_fit):
        def power_law(x, a, b, c, d):
            return (a * ((x-c) ** b)) + d

        initial_a = 1e-6
        initial_b = 1.5
        initial_c = 0
        initial_d = 0
        bounds = ([0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, self.min_mz, np.inf])
        popt, _ = curve_fit(power_law, mz_res_fit, sigma_res_fit, p0=[initial_a, initial_b, initial_c, initial_d], bounds=bounds, maxfev=20000)

        a_fit, b_fit, c_fit, d_fit = popt

        return a_fit, b_fit, c_fit, d_fit


    def calc_min_spacing(self, lst):
        min_diff = float('inf')
        for i in range(len(lst)-1):
            diff = abs(lst[i+1] - lst[i])
            if diff <= min_diff:
                min_diff = diff
        return min_diff


    def gaussian(self, x, amplitude, mean, stddev):
        return amplitude * np.exp(-((x - mean)**2 / (2 * stddev**2)))


    def fit_gaussian(self, x, y, fwhm):
        try:
            # Initial parameter estimates for the Gaussian fit
            amplitude_guess = np.max(y)
            mean_guess = x[np.argmax(y)]
            stddev_guess = fwhm/ (2 * np.sqrt(2 * np.log(2)))  # FWHM to standard deviation conversion

            # Perform curve fitting
            popt, _ = curve_fit(self.gaussian, x, y, p0=[amplitude_guess, mean_guess, stddev_guess])
        except:
            popt = [0,0,0]
        return popt  # Returns the optimized parameters of the Gaussian fit


    def noise_level_calculation(self, peaks):
        num_splits = 50
        split_size = len(peaks) // num_splits
        splits = [peaks[i:i + split_size] for i in range(0, len(peaks), split_size)]

        average_mz_list = []
        modal_intens_list = []

        for split in splits:
            average_mz = np.mean([t[0] for t in split])

            intensities = [t[1] for t in split]
            hist, edges = np.histogram(intensities, bins=500)
            mode_bin_index = np.argmax(hist)
            mode_value = (edges[mode_bin_index] + edges[mode_bin_index + 1]) / 2

            average_mz_list.append(average_mz)
            modal_intens_list.append(mode_value)

        coefficients = np.polyfit(average_mz_list, modal_intens_list, 2)
        a, b, c = coefficients

        return [a, b, c]


    def cwt_centroid(self):
        try:
            mz_values, intensities = self.read_scan_data(self.file_path)
        except:
            mz_values, intensities = self.read_scan_data_space(self.file_path)
        self.make_mzml(mz_values, intensities, 'profile')
        self.make_centroid_mzml()
        mz, intensity = self.mzml_to_xy(self.output_centroid)

        return list(zip(mz, intensity))


    def read_scan_data(self, input_scan):
        mz_values = []
        intensities = []
        with open(input_scan, "r") as file:
            for line in file:
                mz, intensity = line.strip().split("\t")
                mz_values.append(float(mz))
                intensities.append(float(intensity))
        return mz_values, intensities


    def read_scan_data_space(self, input_scan):
        mz_values = []
        intensities = []
        with open(input_scan, "r") as file:
            for line in file:
                mz, intensity = line.strip().split(" ")
                mz_values.append(float(mz))
                intensities.append(float(intensity))
        return mz_values, intensities


    def make_mzml(self, mz_values, intensities, spectrum_type):
        #write a fake mzml file w/o refrencing scan header

        if spectrum_type == "profile":
            writer = MzMLWriter(self.output_profile)
        elif spectrum_type == "centroid":
            writer = MzMLWriter(self.output_centroid)
        with writer:
            writer.controlled_vocabularies()
            writer.file_description(["MSn spectrum"])
            writer.software_list([{"id": "psims-writer", "version": "0.1.2", "params": ["python-psims"]}])
            source = writer.Source(1, ["nanoelectrospray", "nanospray inlet"])
            analyzer = writer.Analyzer(2, ["fourier transform ion cyclotron resonance mass spectrometer"])
            detector = writer.Detector(3, ["inductive detector"])
            config = writer.InstrumentConfiguration(
                id="IC1",
                component_list=[source, analyzer, detector],
                params=["high resolution mass spectrometer"]
            )
            writer.instrument_configuration_list([config])
            methods = []
            methods.append(writer.ProcessingMethod(
                    order=0,
                    software_reference="psims-writer",
                    params=["Conversion to mzML"]
                )
            )
            processing = writer.DataProcessing(methods, id="DP1")
            if spectrum_type == 'profile':
                writer.data_processing_list([processing])
            elif spectrum_type == 'centroid':
                methods = []
                methods.append(writer.ProcessingMethod(
                        order=1,
                        software_reference="psims-writer",
                        params=["MS:1000035", "peak picking"]
                    )
                )
                processing2 = writer.DataProcessing(methods, id="DP2")
                writer.data_processing_list([processing, processing2])

            with writer.run(id=1, instrument_configuration="IC1"):
                with writer.spectrum_list(count=1):
                    if spectrum_type == "profile":
                        writer.write_spectrum(
                            mz_values, intensities, id="controllerType=0 controllerNumber=1 scan=1",
                            centroided=False, scan_start_time=0.01,
                            scan_window_list=[(min(mz_values), max(mz_values))],
                            params=[{"ms level": 2}, {"total ion current": sum(intensities)}]
                        )
                    elif spectrum_type == "centroid":
                        writer.write_spectrum(
                            mz_values, intensities, id="controllerType=0 controllerNumber=1 scan=1",
                            centroided=True, scan_start_time=0.01,
                            scan_window_list=[(min(mz_values), max(mz_values))],
                            params=[{"ms level": 2}, {"total ion current": sum(intensities)}]
                        )
        #need to change ID to make work
        old_id = "1002991"
        new_id = "1000615"

        if spectrum_type == "profile":
            with open(self.output_profile, "r") as file:
                content = file.read()
            modified_content = content.replace(old_id, new_id)
            with open(self.output_profile, "w") as file:
                file.write(modified_content)

        elif spectrum_type == "centroid":
            with open(self.output_centroid, "r") as file:
                content = file.read()
            modified_content = content.replace(old_id, new_id)
            with open(self.output_centroid, "w") as file:
                file.write(modified_content)
                print(f"save centroid file as {self.output_centroid}")


    def make_centroid_mzml(self):
        pwiz_command = [
            "./pwiz/msconvert",
            "--filter",
            "peakPicking cwt snr=0.05 peakSpace=0.005 msLevel=1-",
            self.output_profile,
            "-o",
            self.directory_path,
            "--outfile",
            self.output_centroid
        ]

        try:
            subprocess.run(pwiz_command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error while executing proteowizard: {e}")

        #replace scan description
        old_scan = "single_scan"
        new_scan = "controllerType=0 controllerNumber=1 scan=1"
        with open(self.output_centroid, "r") as file:
            content = file.read()
        modified_content = content.replace(old_scan, new_scan)

        with open(self.output_centroid, "w") as file:
            file.write(modified_content)


    def mzml_to_xy(self, mzml_file):
        reader = mzml.read(mzml_file)
        spectrum = next(reader)
        mzs = np.array(spectrum["m/z array"])
        intensities = np.array(spectrum["intensity array"])

        return mzs, intensities



class CWTPeakPicking():
    def __init__(self, file_path, directory_path, html_file_save):
        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(self.file_path).replace(".txt", "")
        self.output_profile = os.path.join(
            self.directory_path,
            f"{self.basename}_profile.mzML"
        )
        self.output_centroid = os.path.join(
            self.directory_path,
            f"{self.basename}_centroid.mzML"
        )
        self.output_html = os.path.join(
            self.directory_path,
            f"{self.basename}.deconvolution.html"
        )

        cwt_peaks = self.cwt_centroid()


        mz = [item[0] for item in cwt_peaks]
        intensity = [item[1] for item in cwt_peaks]

        self.make_mzml(mz, intensity, 'centroid')

        self.save_data(cwt_peaks)
        if html_file_save:
            self.generate_decon_plot(
                mz,
                intensity,
                self.mz_values,
                self.intensities
            )


    def save_data(self, cwt_peaks):
        output_cwt = os.path.join(self.directory_path, f"{self.basename}.centroid.cwt.txt")

        with open(output_cwt, "w") as file:
            for row in cwt_peaks:
                file.write(f"{row[0]}\t{row[1]}\n")


    def cwt_centroid(self):
        try:
            mz_values, intensities = self.read_scan_data(self.file_path)
        except:
            mz_values, intensities = self.read_scan_data_space(self.file_path)

        self.mz_values = mz_values
        self.intensities = intensities
        self.make_profile_mzml(mz_values, intensities)
        self.make_centroid_mzml()
        mz, intensity = self.mzml_to_xy(self.output_centroid)

        cwt_peaks = list(zip(mz, intensity))
        noise = self.noise_level_calculation(cwt_peaks) # coefficients for quadratic fit of noise level

        filtered_cwt_peaks = []
        for peak in cwt_peaks:
            mz = peak[0]
            noise_level = noise[0] * peak[0] ** 2 + noise[1] * peak[0] + noise[2]

            if peak[1] > 0.1 * noise_level:
                filtered_cwt_peaks.append(peak)

        return filtered_cwt_peaks


    def read_scan_data(self, input_scan):
        mz_values = []
        intensities = []
        with open(input_scan, "r") as file:
            for line in file:
                mz, intensity = line.strip().split("\t")
                mz_values.append(float(mz))
                intensities.append(float(intensity))
        return mz_values, intensities


    def read_scan_data_space(self, input_scan):
        mz_values = []
        intensities = []
        with open(input_scan, "r") as file:
            for line in file:
                mz, intensity = line.strip().split(" ")
                mz_values.append(float(mz))
                intensities.append(float(intensity))
        return mz_values, intensities


    def make_profile_mzml(self, mz_values, intensities):
        #write a fake mzml file w/o refrencing scan header
        writer = MzMLWriter(self.output_profile)
        with writer:
            writer.controlled_vocabularies()
            writer.file_description(["MSn spectrum"])
            writer.software_list([{"id": "psims-writer", "version": "0.1.2", "params": ["python-psims"]}])
            source = writer.Source(1, ["nanoelectrospray", "nanospray inlet"])
            analyzer = writer.Analyzer(2, ["fourier transform ion cyclotron resonance mass spectrometer"])
            detector = writer.Detector(3, ["inductive detector"])
            config = writer.InstrumentConfiguration(
                id="IC1",
                component_list=[source, analyzer, detector],
                params=["high resolution mass spectrometer"]
            )
            writer.instrument_configuration_list([config])
            methods = []
            methods.append(writer.ProcessingMethod(
                    order=1,
                    software_reference="psims-writer",
                    params=["Conversion to mzML"]
                )
            )
            processing = writer.DataProcessing(methods, id="DP1")
            writer.data_processing_list([processing])

            with writer.run(id=1, instrument_configuration="IC1"):
                with writer.spectrum_list(count=1):
                    writer.write_spectrum(
                        mz_values, intensities, id="single_scan",
                        centroided=False, scan_start_time=0.01,
                        scan_window_list=[(min(mz_values), max(mz_values))],
                        params=[{"ms level": 2}, {"total ion current": sum(intensities)}]
                    )

        #need to change ID to make work
        old_id = "1002991"
        new_id = "1000615"
        with open(self.output_profile, "r") as file:
            content = file.read()
        modified_content = content.replace(old_id, new_id)
        with open(self.output_profile, "w") as file:
            file.write(modified_content)


    def make_centroid_mzml(self):
        pwiz_command = [
            "./pwiz/msconvert",
            "--filter",
            "peakPicking cwt snr=0.05 peakSpace=0.005 msLevel=1-",
            self.output_profile,
            "-o",
            self.directory_path,
            "--outfile",
            self.output_centroid
        ]

        try:
            subprocess.run(pwiz_command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error while executing proteowizard: {e}")

        #replace scan description
        old_scan = "single_scan"
        new_scan = "controllerType=0 controllerNumber=1 scan=1"
        with open(self.output_centroid, "r") as file:
            content = file.read()
        modified_content = content.replace(old_scan, new_scan)

        with open(self.output_centroid, "w") as file:
            file.write(modified_content)


    def mzml_to_xy(self, mzml_file):
        reader = mzml.read(mzml_file)
        spectrum = next(reader)
        mzs = np.array(spectrum["m/z array"])
        intensities = np.array(spectrum["intensity array"])

        return mzs, intensities


    def make_mzml(self, mz_values, intensities, spectrum_type):
        #write a fake mzml file w/o refrencing scan header

        if spectrum_type == "profile":
            writer = MzMLWriter(self.output_profile)
        elif spectrum_type == "centroid":
            writer = MzMLWriter(self.output_centroid)
        with writer:
            writer.controlled_vocabularies()
            writer.file_description(["MSn spectrum"])
            writer.software_list([{"id": "psims-writer", "version": "0.1.2", "params": ["python-psims"]}])
            source = writer.Source(1, ["nanoelectrospray", "nanospray inlet"])
            analyzer = writer.Analyzer(2, ["fourier transform ion cyclotron resonance mass spectrometer"])
            detector = writer.Detector(3, ["inductive detector"])
            config = writer.InstrumentConfiguration(
                id="IC1",
                component_list=[source, analyzer, detector],
                params=["high resolution mass spectrometer"]
            )
            writer.instrument_configuration_list([config])
            methods = []
            methods.append(writer.ProcessingMethod(
                    order=0,
                    software_reference="psims-writer",
                    params=["Conversion to mzML"]
                )
            )
            processing = writer.DataProcessing(methods, id="DP1")
            if spectrum_type == 'profile':
                writer.data_processing_list([processing])
            elif spectrum_type == 'centroid':
                methods = []
                methods.append(writer.ProcessingMethod(
                        order=1,
                        software_reference="psims-writer",
                        params=["MS:1000035", "peak picking"]
                    )
                )
                processing2 = writer.DataProcessing(methods, id="DP2")
                writer.data_processing_list([processing, processing2])

            with writer.run(id=1, instrument_configuration="IC1"):
                with writer.spectrum_list(count=1):
                    if spectrum_type == "profile":
                        writer.write_spectrum(
                            mz_values, intensities, id="controllerType=0 controllerNumber=1 scan=1",
                            centroided=False, scan_start_time=0.01,
                            scan_window_list=[(min(mz_values), max(mz_values))],
                            params=[{"ms level": 2}, {"total ion current": sum(intensities)}]
                        )
                    elif spectrum_type == "centroid":
                        writer.write_spectrum(
                            mz_values, intensities, id="controllerType=0 controllerNumber=1 scan=1",
                            centroided=True, scan_start_time=0.01,
                            scan_window_list=[(min(mz_values), max(mz_values))],
                            params=[{"ms level": 2}, {"total ion current": sum(intensities)}]
                        )
        #need to change ID to make work
        old_id = "1002991"
        new_id = "1000615"

        if spectrum_type == "profile":
            with open(self.output_profile, "r") as file:
                content = file.read()
            modified_content = content.replace(old_id, new_id)
            with open(self.output_profile, "w") as file:
                file.write(modified_content)

        elif spectrum_type == "centroid":
            with open(self.output_centroid, "r") as file:
                content = file.read()
            modified_content = content.replace(old_id, new_id)
            with open(self.output_centroid, "w") as file:
                file.write(modified_content)


    def noise_level_calculation(self, peaks):
        num_splits = 50
        split_size = len(peaks) // num_splits
        splits = [peaks[i:i + split_size] for i in range(0, len(peaks), split_size)]

        average_mz_list = []
        modal_intens_list = []

        for split in splits:
            average_mz = np.mean([t[0] for t in split])

            intensities = [t[1] for t in split]
            hist, edges = np.histogram(intensities, bins=500)
            mode_bin_index = np.argmax(hist)
            mode_value = (edges[mode_bin_index] + edges[mode_bin_index + 1]) / 2

            average_mz_list.append(average_mz)
            modal_intens_list.append(mode_value)

        coefficients = np.polyfit(average_mz_list, modal_intens_list, 2)
        a, b, c = coefficients

        return [a, b, c]


    def generate_decon_plot(self, cent_x, cent_y, mz, intensity):
        peaks = go.Scatter(
            x=cent_x,
            y=cent_y,
            name='Peaks',
            mode='markers',
            visible=True,
            marker=dict(color='blue', line=dict(color='white', width=1))
        )

        line_spectrum = go.Scatter(
            x=mz,
            y=intensity,
            mode='lines',
            visible=True,
            name='Spectrum',
            line=dict(color='black')
        )

        all_traces = [line_spectrum, peaks]

        max_y_value = max(intensity)

        # Create layout
        layout = go.Layout(
            title='CWT deconvolution results',
            showlegend=True,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            yaxis=dict(range=[0, 1.1 * max_y_value])
        )

        fig = go.Figure(data=all_traces, layout=layout)

        button_list = []

        button_list.append(
            dict(
                label='Show all',
                method='update',
                args=[{'visible': [True, True]}]
            )
        )
        
        button_list.append(
            dict(
                label='Show raw spectrum',
                method='update',
                args=[{'visible': [True, False]}]
            )
        )


        fig.update_layout(
            updatemenus=[
                dict(
                    type='buttons',
                    showactive=True,
                    buttons=button_list,
                    x=1,
                    xanchor='left',
                    y=0.1,
                    yanchor='top'
                ),
            ]
        )

        fig.write_html(self.output_html)
