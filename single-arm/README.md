This repository holds the codes for comparing the Data Vs Simulation.
Data is Dummy Subtracted.
Simulation is mc-single-arm simulation.
...
Old_Codes repository holds the codes that I created for particular topic, i.e., ratio
(where plotting the histogram ratios were implemented first time), project (where ROOT
Project() method were used instead of Draw() method), etc.
...
DataVsSimPlot_MultiDataMultiDummy.C is the final developped code.
This naming is due to the fact that this code is capable to take multiple data run and 
dummy run to incorporate in the comparison.
...
You must have PDFs directory created prior to run the code. This directory holds the 
created plots. To run the code, run:
root -l DataVsSimPlot_MultiDataMultiDummy.C
