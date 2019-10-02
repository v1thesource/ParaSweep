# ParaSweep
A generalised tool for sensitivity analysis of defect equilibria

The tool goes through a series of steps to conduct a sensitivity analysis for a desired parameter:
  - A datafile of the defect equilibria data at each parameter increment is populated
  - The datafile is parsed and a plot of the defect equilibria at each increment is generated
  - The individual plots are stitched together to create a visualisation of a parameter sweep (ParaSweep!)

Example output:
<p align="center">
  <img src="https://github.com/v1thesource/ParaSweep/blob/master/tet_iodine_example.gif?raw=true">
</p>

In the above example, system temperature was chosen as the parameter to sweep across (500 < T < 2000). The sweeping parameter can, however, be any one of the inputs to the data generation script. 

For example, the band gap of an insulating material:

<p align="center">
  <img src="https://github.com/v1thesource/ParaSweep/blob/master/tet_bandgaps_example.gif?raw=true">
</p>
