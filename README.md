# IntrinsicProperties
This code analyses the intrinsic properties, including 1 Input resistance, 2 Membrane time constant, 3 Sag ratio, 4 First spike latency, 5 Depolarizing hump amplitude, 6 AP amplitude, 7 AP threshold, 8 First AHP latency, 9 AHP half-width, of the recorded cells under current clamp configuration

It reads a .dat file generated from PatchMaster by using HEKA_Patchmaster_Importer (https://github.com/ChristianKeine/HEKA_Patchmaster_Importer). 

Specify the voltage traces for analysis in response to step currents pulses in the pop up window and the pre-defined parameters in the scripts (optional to view the traces, imshow = 1).

It generates an excel table (.xlsx) summarizing the analysed results.
