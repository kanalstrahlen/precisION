from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import wx
import seaborn as sns
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

class AnnotatorFunctions():
    def extract_protein_sequence_from_xml(up_xml_file):
        with open(up_xml_file, "r") as file:
            for record in SeqIO.parse(file, "uniprot-xml"):
                protein_sequence = str(record.seq)
            for feature in record.features:
                if feature.type == "disulfide bond":
                    if isinstance(feature.location.start, int):
                        if isinstance(feature.location.end, int):
                            sequence_list = list(protein_sequence)
                            sequence_list[feature.location.start] = "c"
                            sequence_list[feature.location.end-1] = "c"
                            protein_sequence = "".join(sequence_list)
        return protein_sequence


    def validation_plots(annotated_peak_file, calibration_file):
        peak_list = pd.read_csv(annotated_peak_file)
        peak_list = peak_list.dropna(subset=['name'])
        unique_names = peak_list['name'].unique()

        cal_data = np.loadtxt(calibration_file, skiprows=1)
        cal_x = cal_data[:, 0]
        cal_y = cal_data[:, 1]
        cal_regression_line = cal_data[:, 2]

        for name in unique_names:
            name_data = peak_list[peak_list['name'] == name]
            sequence = name_data['sequence'].tolist()[0]
            _, axs = plt.subplots(
                2,
                3,
                figsize=(12, 10),
                num=f"precisION - Validation Plots [{name}]"
            )

            axs[0, 0].scatter(cal_x, cal_y, color='black', marker='o', alpha=0.5)
            axs[0, 0].plot(cal_x, cal_regression_line, color='black')
            axs[0, 0].set_title("Calibration Curve")
            axs[0, 0].set_xlabel('m/z')
            axs[0, 0].set_ylabel('Error (Th)')

            sns.kdeplot(data=name_data, x='ppm_error', ax=axs[0, 1], color="black")
            sns.rugplot(data=name_data, x="ppm_error", ax=axs[0, 1], color="black")
            axs[0, 1].set_title("Post-calibration Errors")
            axs[0, 1].set_xlabel('Error (ppm)')

            axs[0, 2].hist(
                name_data['fit_score'],
                bins=round((1 / 0.05)),
                alpha=1.0,
                color="black",
                rwidth=0.85
            )
            axs[0, 2].set_title("Isotope Fitting Scores")
            axs[0, 2].set_xlabel('Fit score')
            axs[0, 2].set_ylabel('Frequency')

            ions = name_data['ion'].tolist()
            frag_sites = name_data['frag_site'].tolist()
            abundances = name_data['charge_scaled_abundance'].tolist()

            amino_acid_pairs = []
            corrected_abundances = []

            for i, ion in enumerate(ions):
                ion_type = ion[0]
                if ion_type != 'i':
                    amino_acid_pairs.append(frag_sites[i])
                    corrected_abundances.append(abundances[i])
                elif ion_type == "i":
                    frag_site_list = frag_sites[i].split(";")
                    frag_site_1 = frag_site_list[0]
                    frag_site_2 = frag_site_list[1]
                    amino_acid_pairs.append(frag_site_1)
                    corrected_abundances.append(abundances[i] / 2)
                    amino_acid_pairs.append(frag_site_2)
                    corrected_abundances.append(abundances[i] / 2)
    
            tuple_list = list(zip(amino_acid_pairs, corrected_abundances))
    
            all_amino_acids = [
                'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
            ]

            heatmap_data = np.zeros((len(all_amino_acids), len(all_amino_acids)))
            for pair, intensity in tuple_list:
                if type(pair) != float:
                    i, j = all_amino_acids.index(pair[0]), all_amino_acids.index(pair[1])
                    heatmap_data[i, j] += intensity

            for i, _ in enumerate(all_amino_acids):
                for j, _ in enumerate(all_amino_acids):
                    pair = all_amino_acids[i] + all_amino_acids[j]
                    num_pair = sequence.count(pair)
                    if num_pair >= 1:
                        heatmap_data[i, j] = heatmap_data[i, j] / sequence.count(pair)

            max_intensity = 0
            for i, _ in enumerate(all_amino_acids):
                for j, _ in enumerate(all_amino_acids):
                    if heatmap_data[i, j] >= max_intensity:
                        max_intensity = heatmap_data[i, j]
            for i, _ in enumerate(all_amino_acids):
                for j, _ in enumerate(all_amino_acids):
                    heatmap_data[i, j] = heatmap_data[i, j] / max_intensity * 100

            custom_cmap = LinearSegmentedColormap.from_list(
                'custom_colormap',
                ['white', 'black'],
                N=256
            )
            heat_map = axs[1, 0].imshow(heatmap_data, cmap=custom_cmap, vmin=0)
            axs[1, 0].set_xticks(range(len(all_amino_acids)), all_amino_acids)
            axs[1, 0].set_yticks(range(len(all_amino_acids)), all_amino_acids)
            axs[1, 0].set_xlabel('C-terminal residue')
            axs[1, 0].set_ylabel('N-terminal residue')
            axs[1, 0].set_title('Fragmentation Site Heatmap')

            cbar = plt.colorbar(heat_map, ax=axs[1, 0], orientation='horizontal', fraction = 0.05)
            cbar.set_label('Relative, normalised intensity (%)')

            b_ion_data = name_data[name_data['ion'].str.startswith(('b', 'c'), na=False)]
            axs[1, 1].scatter(
                b_ion_data["monoisotopic_mw"],
                b_ion_data["charge"],
                facecolor='#205283',
                edgecolors='black',
                marker="s",
                alpha=0.5
            )
            if len(b_ion_data["charge"])>=1:
                axs[1, 1].set_ylim(min(b_ion_data["charge"]) - 0.1, max(b_ion_data["charge"]) + 0.1)
                axs[1, 1].set_yticks(range(
                    int(min(b_ion_data["charge"])),
                    int(max(b_ion_data["charge"])+1)
                    )
                )
            axs[1, 1].set_xlabel('Mass (Da)')
            axs[1, 1].set_ylabel('Charge')
            axs[1, 1].set_title("b/c-type Ion Mass vs Charge")

            y_ion_data = name_data[name_data['ion'].str.startswith(('y', 'z'), na=False)]
            axs[1, 2].scatter(
                y_ion_data["monoisotopic_mw"],
                y_ion_data["charge"],
                facecolor='#b42920',
                edgecolors='black',
                marker="s",
                alpha=0.5
            )
            if len(y_ion_data["charge"])>=1:
                axs[1, 2].set_ylim(
                    min(y_ion_data["charge"]) - 0.1,
                    max(y_ion_data["charge"]) + 0.1
                )
                axs[1, 2].set_yticks(range(
                    int(min(y_ion_data["charge"])),
                    int(max(y_ion_data["charge"])+1)
                    )
                )
            axs[1, 2].set_xlabel('Mass (Da)')
            axs[1, 2].set_ylabel('Charge')
            axs[1, 2].set_title("y/z-type Ion Mass vs Charge")
            x_limits = axs[1, 2].get_xlim()

            axs[1, 2].set_xlim(x_limits[::-1])

            plt.tight_layout(rect=[0, 0, 1, 0.97])

            fig_manager = plt.get_current_fig_manager()
            fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

            plt.show()
