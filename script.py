import json
import time
from pprint import pprint

import numpy as np
import requests
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

class Analysis:
    def __init__(self, params_file=None, lambda_col=0, intensities_col=1):
        self.params = None
        self.lambdas = None
        self.intensities = None
        self.gammas = None
        self.elements = {}
        self.lorentzian_functions = []
        self.analysis_options = {
            "adjust_manually": False,
            "minimal_separation_between_maximum": 2,
            "iterations": 2,
            "intensity_detection_umbral": 200,
            "spectra_output_file": "spectra_output.txt",
            "maximum_output_file": "maximum_output.txt"
        }
        if params_file:
            self.read_params_file(params_file, lambda_col=lambda_col, intensities_col=intensities_col)

    def read_params_file(self, txt_file, **kwargs):
        with open(txt_file, "r") as f_params:
            self.params = json.loads(f_params.read())
            if 'file' in self.params:
                print("Reading experiment file...")
                data = self.read_measurement_file(self.params.get("file"))
                self.lambdas = data[:, kwargs.get("lambda_col", 0)]
                self.intensities = data[:, kwargs.get("intensities_col", 1)]

            if 'elements_from_file' in self.params:
                pass

            if 'elements_from_NIST' in self.params:
                print("Getting elements from NIST website...")
                for element in self.params.get('elements_from_NIST'):
                    self.elements[element] = self.scrape_data_from_NIST(element)
            if "analysis_options" in self.params:
                for option in self.params['analysis_options']:
                    self.analysis_options[option] = self.params['analysis_options'][option]
                print("Analysis configuration:")
                pprint(self.analysis_options)

            if "wavelength_min" in self.analysis_options or "wavelength_max" in self.analysis_options:
                self.limit_wavelengths()

                # .. more options

    def limit_wavelengths(self):
        if "wavelength_min" not in self.analysis_options:
            self.analysis_options["wavelength_min"] = min(self.lambdas)
        if "wavelength_max" not in self.analysis_options:
            self.analysis_options["wavelength_max"] = max(self.lambdas)

        wavelengths = []
        intensities = []
        for wavelength, intensity in zip(self.lambdas, self.intensities):
            if self.analysis_options["wavelength_min"] <= wavelength <= self.analysis_options["wavelength_max"]:
                wavelengths.append(wavelength)
                intensities.append(intensity)
        self.lambdas = wavelengths
        self.intensities = intensities

    def read_measurement_file(self, txt_file):
        data = np.loadtxt(txt_file)
        return data

    def scrape_data_from_NIST(self, element_str):
        element_arr = []
        element_URI = element_str.replace(" ", "+")
        url = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=" + element_URI
        url += "&low_wl=200&upp_wn=&upp_wl=900&low_wn=&unit=1&submit=Retrieve+Data&de=0&java_window=3&java_mult=&" \
               "format=0&line_out=0&en_unit=1&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&order_out=0" \
               "&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1" \
               "&forbid_out=1&min_accur=&min_intens=&enrg_out=on&J_out=on&g_out=on"
        headers = {
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Encoding": "gzip, deflate, sdch, br",
            "Accept-Language": "en",
            "Connection": "keep-alive",
            "Cookie": "__utmt_GSA_CP1=1; __utmt_GSA_CP2=1; __utma=259892043.1270432421.1502294478.1502678601.150267860"
                      "1.1; __utmb=259892043.8.10.1502678601; __utmc=259892043; __utmz=259892043.1502678601.1.1.utmcsr"
                      "=(direct)|utmccn=(direct)|utmcmd=(none); _gat_GSA_ENOR0=1; _ga=GA1.3.1270432421.1502294478; _gid"
                      "=GA1.3.513377169.1502678526; _gat_GSA_ENOR1=1",
            "Host": "physics.nist.gov",
            "Referer": "https://physics.nist.gov/PhysRefData/ASD/lines_form.html",
            "Upgrade-Insecure-Requests": "1",
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/57.0.2987.133 Safari/537.36"
        }
        r = requests.get(url, headers=headers)
        soup = BeautifulSoup(r.text, 'lxml')
        trows = soup.find_all('tr', 'evn') + soup.find_all('tr', 'odd')
        for trow in trows:
            data = {}
            for i, td in enumerate(trow.findAll('td')):
                if i == 0:
                    data['observed_wavelength_air'] = self._get_float_number(td.get_text())
                elif i == 1:
                    data['ritz_wavelength_air'] = self._get_float_number(td.get_text())
                elif i == 2:
                    data['relative_intensity'] = self._get_float_number(td.get_text())
                elif i == 3:
                    data['Aki'] = self._get_float_number(td.get_text())
                elif i == 5:
                    data['Ei'] = self._get_float_number(td.get_text())
                elif i == 7:
                    data['Ek'] = self._get_float_number(td.get_text())
                elif i == 8:
                    data['gi'] = self._get_float_number(td.get_text())
                elif i == 10:
                    data['gk'] = self._get_float_number(td.get_text())
            element_arr.append(data)
        return element_arr

    def _get_float_number(self, string):
        res = ""
        for ch in string:
            if ch.isdigit() or ch == ".":
                res += ch
        if not res:
            return None
        return float(res)

    def plot_result(self):
        # Create Figure (empty canvas)
        fig = plt.figure()

        # Add set of axes to figure
        ax1 = fig.add_subplot(211)

        res = self.params.get("result", {"resolution": 5}).get("resolution", 5) # times the original resolution
        x = np.linspace(self.lambdas[0], self.lambdas[-1], len(self.lambdas)*res)
        y_total = np.zeros(len(self.lambdas)*res)

        for lorentzian_function in self.lorentzian_functions:
            func_center = lorentzian_function["wavelength"]
            y = lorentzian_function["intensity"] / (((func_center - x) / lorentzian_function["gamma"]) ** 2 + 1)
            y_total += y
            ax1.plot(x, y, 'g')
        ax1.plot(x, y_total, ':b')

        # Plot on that set of axes
        ax1.set_title('Result')

        ax2 = fig.add_subplot(212, sharex=ax1, sharey=ax1)
        ax2.plot(self.lambdas, self.intensities, 'b')
        ax2.set_xlabel('Wavelength (nm)')
        ax2.set_ylabel('Intensity')

        plt.show()

    def _find_lorentzian_functions(self):

        # find maximum
        intensities = np.array(self.intensities)
        isAscending = False
        msbm = self.analysis_options["minimal_separation_between_maximum"]
        for i in range(len(self.intensities) - msbm):
            data_above = sum(intensities[i:i + msbm + 1] >= intensities[i])
            if data_above == msbm + 1:
                isAscending = True
            elif data_above == 1 and isAscending:
                isAscending = False
                if self.intensities[i] > self.analysis_options["intensity_detection_umbral"]:
                    lorentzian_function = {
                        "wavelength": self.lambdas[i],
                        "intensity": self.intensities[i],
                        "index": i,
                        "gamma": None
                    }
                    self.lorentzian_functions.append(lorentzian_function)

        # edit maximum found
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(211)
        self.axes.set_ylabel("Intenisity")
        self.axes.set_title("Maximum found")
        self.axes.plot(self.lambdas, self.intensities, ':og', picker=5, markersize=3)

        self.maximum_lines_added  = {}
        for maximum_found in self.lorentzian_functions:
            self.maximum_lines_added[maximum_found['index']] = self.axes.stem((maximum_found['wavelength'],), (maximum_found['intensity'],))

        # plot elements
        elements_axes = self.fig.add_subplot(212, sharex=self.axes, sharey=self.axes)
        elements_axes.set_xlabel("Wavelength (nm)")
        elements_axes.set_title("Elements database")
        max_intensity = max(self.intensities)
        min_wavelength = min(self.lambdas)
        max_wavelength = max(self.lambdas)
        for element in self.elements:
            wavelengths = []
            intensities = []
            for data in self.elements[element]:
                owa = data.get("observed_wavelength_air")
                rel_int = data.get("relative_intensity")
                if owa and rel_int:
                    if min_wavelength < owa < max_wavelength:
                        wavelengths.append(data.get("observed_wavelength_air"))
                        intensities.append(data.get("relative_intensity"))
            if intensities:
                max_local_intensity = max(intensities)
                intensities = np.array(intensities)/max_local_intensity*max_intensity

                intensities_tmp = intensities[intensities > self.params['analysis_options'].get('elements_min_intensity', 0)]
                wavelengths = np.array(wavelengths)[intensities > self.params['analysis_options'].get('elements_min_intensity', 0)]
                elements_axes.stem(wavelengths, intensities_tmp, ':b')

        self.fig.canvas.mpl_connect('pick_event', self._remove_or_add_lines_onpick)

        plt.show()

        # adjust a gamma value
        for i, lorentzian_function in enumerate(self.lorentzian_functions):
            index = lorentzian_function.get("index")
            maximum_intensity = lorentzian_function["intensity"]
            intensity_right = self.intensities[index + 1]
            intensity_left = self.intensities[index - 1]
            if intensity_right < intensity_left:
                xn = self.lambdas[index + 1]
                intensity_selected = intensity_right
            else:
                xn = self.lambdas[index - 1]
                intensity_selected = intensity_left

            if intensity_selected <= 0:
                # el 0.000000001 para evitar division entre 0
                maximum_intensity += abs(intensity_selected) + 0.000000001
                intensity_selected += abs(intensity_selected) + 0.000000001

            # From lorentzian function definition
            gamma = abs(self.lambdas[index] - xn) / np.sqrt((maximum_intensity / intensity_selected) - 1.0)
            self.lorentzian_functions[i]["gamma"] = gamma

        # calculate coefficient values to adjust lorentzian functions intensity
        coeff = self._calculate_coefficients()
        for _ in range(self.analysis_options["iterations"] - 1):
            coeff[np.isnan(coeff)] = 1
            self.lorentzian_functions = [self.lorentzian_functions[i] for i in range(len(self.lorentzian_functions)) if
                                         0 <= coeff[i] <= 1]
            coeff = self._calculate_coefficients()

        for i in range(len(self.lorentzian_functions)):
            self.lorentzian_functions[i]["intensity"] *= coeff[i]

    def _calculate_coefficients(self):
        A = []
        b = []
        l = len(self.lorentzian_functions)
        for lorentzian_function in self.lorentzian_functions:
            A.append(
                [self._calculate_functions_value_in_x0(self.lorentzian_functions[j], lorentzian_function["wavelength"])
                 for j in range(l)])
            b.append(lorentzian_function["intensity"])
        coeff = np.linalg.solve(A, b)
        return coeff

    def _remove_or_add_lines_onpick(self, event):
        ind = event.ind[0]
        if ind in self.maximum_lines_added:
            for line in self.maximum_lines_added[ind]:
                # if line is an array or tuple, iterate and remove all lines inside
                if hasattr(line, '__iter__'):
                    for l in line:
                        self.axes.lines.remove(l)
                else:
                    self.axes.lines.remove(line)

            # update maximum_lines_added
            self.maximum_lines_added.pop(ind)

            # erase from list of lorentzian functions found
            for i, lorentzian in enumerate(self.lorentzian_functions):
                if lorentzian['index'] == ind:
                    self.lorentzian_functions.pop(i)
                    break

        else:
            lorentzian_function = {
                "wavelength": self.lambdas[ind],
                "intensity": self.intensities[ind],
                "index": ind,
                "gamma": None
            }
            # add to lorentzian functions list
            self.lorentzian_functions.append(lorentzian_function)

            # update maximum_lines_added
            self.maximum_lines_added[ind] = self.axes.stem((lorentzian_function['wavelength'],), (lorentzian_function['intensity'],))

        # update plot
        self.fig.canvas.draw()

    def _calculate_functions_value_in_x0(self, lorentzian_function, x0):
        value = lorentzian_function["intensity"] / (
        ((lorentzian_function['wavelength'] - x0) / lorentzian_function["gamma"]) ** 2 + 1)
        return round(value, 4)

    def run_analysis(self):
        init = time.time()
        print("\nStarting analysis...\n")
        intensity_min = min([intensity for intensity in self.intensities if intensity >= 0])
        self.intensities = [intensity - intensity_min for intensity in self.intensities]

        self._find_lorentzian_functions()

        self.plot_result()

        print("%d lorentzian functions found" % len(self.lorentzian_functions))


if __name__ == "__main__":
    analysis = Analysis("params.txt")
    analysis.run_analysis()
