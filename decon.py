import os
import threading
import subprocess
import wx
import pandas as pd
from peakProperties import PeakProperties

# to do: https://stackoverflow.com/questions/20295646/python-ascii-plots-in-terminal

class Cluster:
    def __init__(self, file_path, directory_path, cluster_ppm):
        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(self.file_path).replace(".txt", "")
        self.output_cluster_csv = os.path.join(self.directory_path, f"{self.basename}.clusters.csv")

        self.cluster_ppm = float(cluster_ppm)
        
        # list of tuples containing key info about each peak from thrash
        tuple_list = self.make_tuple_list()
        # add topfd info
        tuple_list.extend(self.read_topfd_file())
        # join peaks with same mass (+/- ppm) and charge
        clusters = self.cluster_tuples_by_mass_and_charge(tuple_list, self.cluster_ppm)
        self.get_cluster_summary(clusters)


    def read_thrash_file(self, file_name):
        input_file_path = os.path.join(self.directory_path, f"{self.basename}.{file_name}.csv")
        low_score_file_path = os.path.join(self.directory_path, f"{self.basename}.thrash.0.4.csv")
        if os.path.exists(input_file_path):
            df = pd.read_csv(input_file_path)
        elif os.path.exists(low_score_file_path):
            df = pd.read_csv(low_score_file_path)
            score = float(file_name[-3:])
            df = df[df['fit'] <= score]

        output_list = list(
            zip(
                df["monoisotopic_mw"],
                df["charge"],
                df["abundance"],
                df["fit"],
                df["interference_score"]
            )
        )
        return output_list


    def make_tuple_list(self):
        list_names = ["thrash.0.1", "thrash.0.2", "thrash.0.3", "thrash.0.4"]
        full_data = []

        for list_name in list_names:
            data_list = [
                (mw, charge, abundance, fit, interference, list_name)
                for mw, charge, abundance, fit, interference in self.read_thrash_file(list_name)
            ]
            full_data.extend(data_list)
        return full_data


    def read_topfd_file(self):
        input_file_path = os.path.join(self.directory_path, f"{self.basename}.msdeconv.csv")
        with open(input_file_path, "r") as file:
            data = file.readlines()

        start_index = data.index("PRECURSOR_INTENSITY=0.00\n") + 1
        end_index = data.index("END IONS\n")

        spectral_data = data[start_index:end_index]
        spectral_list = []
        for line in spectral_data:
            values = line.strip().split("\t")
            monoisotopic_mw = float(values[0])
            abundance = float(values[1])
            charge = int(values[2])
            spectral_list.append((
                monoisotopic_mw,
                charge,
                abundance,
                float("NaN"),
                float("NaN"),
                "msdeconv"
            ))

        input_file_path = os.path.join(self.directory_path, f"{self.basename}.envcnn.csv")
        with open(input_file_path, "r") as file:
            data = file.readlines()

        start_index = data.index("PRECURSOR_INTENSITY=0.00\n") + 1
        end_index = data.index("END IONS\n")

        spectral_data = data[start_index:end_index]
        for line in spectral_data:
            values = line.strip().split("\t")
            monoisotopic_mw = float(values[0])
            abundance = float(values[1])
            charge = int(values[2])
            spectral_list.append((
                monoisotopic_mw,
                charge,
                abundance,
                float("NaN"),
                float("NaN"),
                "envcnn"
            ))
        return spectral_list


    def cluster_tuples_by_mass_and_charge(self, tuples_list, ppm_tolerance=10):
        def calculate_ppm_diff(mass1, mass2):
            return abs((mass1 - mass2) / mass2) * 1e6

        clustered_tuples = []
        for tup in tuples_list:
            mass, charge = tup[0], tup[1]
            if not clustered_tuples:
                clustered_tuples.append([tup])
            else:
                found_cluster = False
                for cluster in clustered_tuples:
                    representative_tup = cluster[0]
                    cluster_mass, cluster_charge = representative_tup[0], representative_tup[1]
                    ppm_diff = calculate_ppm_diff(mass, cluster_mass)
                    if ppm_diff <= ppm_tolerance and charge == cluster_charge:
                        cluster.append(tup)
                        found_cluster = True
                        break
                if not found_cluster:
                    clustered_tuples.append([tup])
        return clustered_tuples


    def get_cluster_summary(self, clustered_tuples):
        summary_list = []

        for cluster in clustered_tuples:
            first_item = cluster[0]

            monoisotopic_mw = first_item[0]
            charge = first_item[1]

            detected_thrash_0_1 = any(item[5] == "thrash.0.1" for item in cluster)
            detected_thrash_0_2 = any(item[5] == "thrash.0.2" for item in cluster)
            detected_thrash_0_3 = any(item[5] == "thrash.0.3" for item in cluster)
            detected_thrash_0_4 = any(item[5] == "thrash.0.4" for item in cluster)
            detected_msdeconv = any(item[5] == "msdeconv" for item in cluster)
            detected_envcnn = any(item[5] == "envcnn" for item in cluster)

            cluster_summary = (
                monoisotopic_mw,
                charge,
                detected_thrash_0_1,
                detected_thrash_0_2,
                detected_thrash_0_3,
                detected_thrash_0_4,
                detected_msdeconv,
                detected_envcnn
            )
            summary_list.append(cluster_summary)

        df = pd.DataFrame(summary_list, columns=[
            "monoisotopic_mw",
            "charge",
            "detected_thrash_0_1?",
            "detected_thrash_0_2?",
            "detected_thrash_0_3?",
            "detected_thrash_0_4?",
            "detected_msdeconv?",
            "detected_envcnn?"
        ])

        df.to_csv(self.output_cluster_csv, index=False)



class Thrash():
    def __init__(self, file_path, directory_path, decon_params, file_base, multi=True):
        self.file_path = file_path
        self.directory_path = directory_path
        self.multi = multi
        self.file_base = file_base
        self.basename = os.path.basename(self.file_path).replace(".txt", "")
        self.thrash_file = ""
        self.thrash(decon_params)


    def thrash(self, decon_params):
        decontools_path = "./DeconTools_x64/DeconConsole.exe"
        signal_to_noise = decon_params[0]
        peak_background_ratio = decon_params[1]
        threshold = decon_params[2]
        peptide_background_ratio = decon_params[3]
        max_charge = decon_params[4]
        max_mw = decon_params[5]

        if self.multi:
            fits = ["0.1", "0.2", "0.3", "0.4"]
        else:
            fits = ["0.4"]

        for max_score in fits:
            print(f"Running THRASH deconvolution algorithim with a maximum score of {max_score}")
            self.create_thrash_param_file(
                signal_to_noise,
                peak_background_ratio,
                threshold,
                peptide_background_ratio,
                max_charge,
                max_mw,
                max_score
            )

            subprocess.run(
                [decontools_path, self.file_path, self.thrash_file, self.directory_path],
                check=False
            )

            os.remove(os.path.join(self.directory_path, f"{self.basename}_log.txt"))
            os.remove(os.path.join(self.directory_path, f"{self.basename}_peaks.txt"))
            os.remove(os.path.join(self.directory_path, f"{self.basename}_scans.csv"))
            try:
                os.remove(
                    os.path.join(
                        self.directory_path,
                        f"{self.file_base}.thrash.{max_score}.csv"
                    )
                )
            except:
                pass

            os.rename(
                os.path.join(self.directory_path,  f"{self.basename}_isos.csv"),
                os.path.join(self.directory_path,  f"{self.file_base}.thrash.{max_score}.csv")
            )
            os.remove(self.thrash_file)


    def create_thrash_param_file(self, signal_to_noise, peak_background_ratio, threshold,
                          peptide_background_ratio, max_charge, max_mw, max_score):
        with open("templateParameterFile.xml", "r") as template:
            template_param_file = template.read()
            replacements = {
                "$SIGNALTONOISETHRESHOLD$": signal_to_noise,
                "$PEAKBACKGROUNDRATIO$": peak_background_ratio,
                "$THRESHOLD$": threshold,
                "$PEPTIDEMINBACKGROUNDRATIO$": peptide_background_ratio,
                "$MAXCHARGE$": max_charge,
                "$MAXMW$": max_mw,
                "$MAXFIT$": max_score
            }

            for key, value in replacements.items():
                template_param_file = template_param_file.replace(key, value)

        self.thrash_file = os.path.join(
            self.directory_path,
            f"{self.basename}.thrashParameterFile.xml"
        )

        with open(self.thrash_file, "w") as output:
            output.write(template_param_file)



class TopFD():
    def __init__(self, file_path, directory_path, decon_params):
        self.file_path = file_path
        self.directory_path = directory_path
        self.decon_params = decon_params
        self.basename = os.path.basename(self.file_path).replace("_centroid.mzML", "")

        print("Running TopFD deconvolution algorithim with Ms-Deconv and EnvCNN scores.")

        for algo in ["msdeconv", "envcnn"]:
            self.topfd(algo)

            if algo == "envcnn":
                save_path = os.path.join(self.directory_path,  f"{self.basename}.envcnn.csv")
                if os.path.exists(save_path):
                    os.remove(save_path)
                os.rename(
                        os.path.join(self.directory_path,  f"{self.basename}_centroid_ms2.msalign"),
                        os.path.join(self.directory_path,  f"{self.basename}.envcnn.csv")
                    )
            elif algo == "msdeconv":
                save_path = os.path.join(self.directory_path,  f"{self.basename}.msdeconv.csv")
                if os.path.exists(save_path):
                    os.remove(save_path)
                os.rename(
                        os.path.join(self.directory_path,  f"{self.basename}_centroid_ms2.msalign"),
                        os.path.join(self.directory_path,  f"{self.basename}.msdeconv.csv")
                    )


    def topfd(self, algo):
        activation_type = self.decon_params[0]
        max_charge = self.decon_params[1]
        max_mass = self.decon_params[2]
        mz_error = self.decon_params[3]
        s2n = self.decon_params[4]

        if algo == "msdeconv":
            topfd_command = [
                "./topFD/topfd.exe",
                self.file_path,
                "-a",
                activation_type,
                "-s",
                s2n,
                "-m",
                max_mass,
                "-t",
                mz_error,
                "-c",
                max_charge,
                "-o",
                "-g",
                "-u",
                "1",
                "-n"
            ]

        elif algo == "envcnn":
            topfd_command = [
                "./topFD/topfd.exe",
                self.file_path,
                "-a",
                activation_type,
                "-s",
                s2n,
                "-m",
                max_mass,
                "-t",
                mz_error,
                "-c",
                max_charge,
                "-o",
                "-g",
                "-u",
                "1"
            ]

        subprocess.run(topfd_command, check=True)
