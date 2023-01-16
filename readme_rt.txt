There are several files:
- input.py: all the input variables needed. If you want to change the output angle, you modify the variable "output_angle" from this file. In degrees and with respect to y-axis (the axis of symmetry)
Also in this file you can modify other variables: N (number of rays), frequency, and the shape of the dome.
- main.py: run this file to run the RT, it will call the other files.
- direct_rayTracing.py: the older version for the direct RT (I think it's not useful anymore)
- rayTracing.py: old function for the direct RT, used before implementig the recursive function
- rayTracingRecursive.py: direct RT implemented with vectors and the recursive function
- radPat.py: calculates the Efield for the rapat

- Phase_distribution_xdeg.xlsx: phase distribution for the array calculated with the direct RT, this is used as an input for COMSOL
- Reverse_anglesIn_x.xlsx: obtained from the revserse RT. Those are the angles we need for the direct RT, to obtain parallel rays at the output
- RT_radpat_xdeg.xlsx: the radiaton pattern calculated from RT, later is imported to matlab. 