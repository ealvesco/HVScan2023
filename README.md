# HVScan2023
Repository for the 2023High-Voltage scans campaing for Resistive Plate Chambers @ CMS

This HVScan repository contains a set of scripts you can run depending on the desired input/output

Original Folder structure:
HVScan2023
|__data
|__macros

1. To start the analysis you first need to get the input files (.root) and the corresponding table of effective x set High Voltages (text file)

For example: in 2022, the input files where located in:

					ls -lht	/eos/cms/store/group/dpg_rpc/comm_rpc/Run-III/HVscan2022/Input_files/

with each file having a name standard such as 'AnalyzeEfficiency_HVscan_86_RPCMon2022.root. Copy these files to the '/data' directory. Also in the /data
directory, create a "hvEffective.txt" file, containing the corresponding high-voltages of the acquisiton. (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/RPCHVScan2022)

2. Run the Fit tool

After certifying that you have the necessary inputs, go to the '/macro' directory to run the "FitData.C" script, which must be run separately for the barrel and the endcap.

If you want the barrel, run: root -l -q 'FitData.C("barrel")'
If you want the endcap, run: root -l -q 'FitData.C("endcap")'

It creates the "results" directory and one subdirectory for each RPC Chamber, each one containing the fit results of efficiency (fitData.txt) and cluster size (fitDataCls.txt).

The meaning of each column is like the following example:

fitData.txt: fit parameters of efficiency X HV plots
9.36312 196.677 97.0909 8.89312 6.19747 94.9837 0.151692 0.102598 0.00234352 -8.10282
Wp      slope50 emax    hv50    chi2    effwp   emaxerr  slopeerr hv50err    slope

fitDataCls.txt: fit parameters of cluster size x HV
78.6441 -13.1354 -0.0142153 0.015002 46.1648 9.36312 2.01312
a       b        c          d        chi2    Wp      clswp


This macro takes the input files and makes the distributions: efficiency x HV (Sigmoid function) and efficiency x HV (Polynomial function).

3. Run the plot producer

The plots of distributions and the fit over it is produced by the pngProducer.C, that can be run separately as well for the barrel and endcap.

					root -l -b -q 'pngProducer.C("barrel")'
					root -l -b -q 'pngProducer.C("endcap")'

However, if you want to produce only selected plots, you can use, for instance: 

pngProducer_splitWpEndcapMaxCut.C - To produce plots filtered by working points in Endcap
pngProducer_splitWpBarrelMaxCut.C - To produce plots filtered by working points in Barrel
pngProducer_splitEmaxCut.C - To produce plots filtered by efficiency

4. Run a summary of fit results

You can make summary plots for all chambers running "MakeASummary.C":

					root -l 'MakeASummary.C("barrel", false)'
					root -l 'MakeASummary.C("endcap", false)'

where the bool is a parameter related to a blacklist you can create, containing a list of the problematic rolls you want to exclude from the analysis. This list of blacklisted 
chambers must be created after the guidance of an expert. If you need to run it without a blacklist, just switch the bool from 'false' to 'true'.

This will create a root file (for instance, endcap_summary_2022.root in the '/summary' directory.

To finish it, run 'root -l DrawingOUFlow.C' and you will have the plots in a special format.


