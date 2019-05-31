Open RunExamples.m and choose the options and duration of the
simulation.

run RunExamples.m

save the output of RunExamples in a .mat file

Open PlotExample.m and rename the .mat file with the name you used to
save the output.

run PlotExample.m

You should see 3 figures: example of the voltage traces, a rastergram
with current input and LFPs (wide-band and filtered) and a histogram
of firing time binned in 5ms bins compared to the LFP trace.

use defaultparamsCA1.m to change the intrinsic properties of cells and
synapses.  use CountRipples.m to change the criteria to find a ripple
start and end (if you want to use a different threshold and see what
happens)

Note that CA1 will not receive non-zero current input in times earlier
than 1s. This is to let the first secod of simulation relax to a
network state before incoming CA3 input induces ripples.  Hence, set
starts to be always bigger than 1.
