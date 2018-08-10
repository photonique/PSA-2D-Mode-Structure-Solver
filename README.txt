
                 Phase-sensitive Amplfier 2D Mode Structure Solver
                        Author: Muthiah Annamalai

Solver
======
This package contains code for calculating the gain and the eigenmode spectrum
of X(2) [chi-2] non-linear amplifier configured for phase-sensitive amplification [in a degenerate mode]. 

You may find the details of the MATLAB solvers below.

Acknowledgements
================
This code is based upon work funded by DARPA’s Quantum Sensor Program under AFRL Contract No. FA8750-09-C-0194. Any opinions, findings, conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of DARPA or the U.S. Air Force.

I would like to also thank Dr. Michael Vasilyev, who supervised the dissertation work.

LICENSE:
========
The code in this package may be used and shared under terms of the MIT license,
unless otherwise clearly marked for use under a different license. Certain files in the Common directory are marked for use under terms of GPL.

Citation:
=========
If you use this code, I appreciate citations to any of the following research articles.

1) Thesis work of Muthiah Annamalai,
'Mode Structure of Noiseless Phase-Sensitive Image Amplifier,' (Dec, 2011).
 University of Texas at Arlington, USA.
   
2)  Muthiah Annamalai, Nikolai Stelmakh, Michael Vasilyev, and Prem Kumar, "Spatial modes of phase-sensitive parametric image amplifiers with circular and elliptical Gaussian pumps," Opt. Express 19, 26710-26724 (2011)
    URL: https://www.osapublishing.org/oe/abstract.cfm?uri=oe-19-27-26710

3) 'Compact representation of the spatial models of a phase-sensitive
    image amplifier,' M. Annamalai et-al, Optics Express (21) 23 pp. 28134-28153
    URL: https://www.osapublishing.org/oe/abstract.cfm?uri=oe-21-23-28134 

PACKAGE:
========

The following files are found in the package, and naming is generally
self-descriptive.

0. Solvers are found in folders with prefix 'HG00-' or 'HGmn' or 'LG00'.
    The spot size ratio between the spot size of signal:pump = Sqrt[2] or fs,
    which is indicated in the folder name.

  +++++++
  CAUTION
  +++++++
  0.1 Please make sure to addpath('Common') folder in MATLAB to pull utility functions for the solver.

  0.2 If you are using the higher order solvers you have to generate the 
    overlap integrals via Mathematica for added precision using the 
    symbolic capability. This is mandatory to proceed.

1. To run the codes you need to add the 'Common/' files to the path.
2. Some solvers 'HGmn' and 'HG00-pump-gen-scale/' require you to generate
    the data files from the 'Mathematica/' folder.
3. The folder 'Post-processing/' contains code to do miscellainous operations like
    visualization of eigenmodes in xy-domain, overlap of TEM00, and overlap with 
    H.O HGmn modes. Sometimes you may be required to generate the HGmn in Sqrt[2]
    basis (or LGpl in Sqrt[2] basis), using the Mathematica code for every order in calculations.
4. Known issues:
    -> All solvers are for (Mx, My) both even
    -> MV10 solver is slowly inaccurate for fs != Sqrt[2]
    -> Mathematica codes are required to produce data files for 
        some solvers and some post-processing tools.

============================= Contents =================================

TOP LEVEL (directory structure)
.:
README.txt
Mathematica/
Common/
Post-processing/
HG00-pump-gen-scale-solver/
HG00-pump-sqrt2-scale-solver/
HGmn-pump-sqrt2-scale-solver/
LG00-pump-sqrt2-scale-solver/

./Mathematica:
higher_order_pump_coupling_integral.nb
HO-pump-mode-overlap.nb
postproc-lgpoly-gen.nb
postproc-hgpoly-gen.nb

./Common:
load_higherorder_overlap.m
blockget.m
blockset.m
factorial2.m
hermitepoly.m
kronsum.m
laguerrepoly.m

./Post-processing:
calc_higherorder_overlap.m
elliptic_overlap.m
extractAmn.m
genhgmode.m
hgmode.m
modevis.m
overlap.m
overlap_sqrt2_unscaled.m
overlap_zoffset.m
postproc_bbnmodevis.m
postproc_coeff_extract.m
postproc.m
postproc_modevis.m
postprocmv11.m
postproc_MV11.m
postproc_renorm.m
spatialmodevis.m
spotsize_scalefactor.m
spotsize_scalefactor_woffset.m

./HG00-pump-gen-scale-solver:
MV10.m

./HG00-pump-sqrt2-scale-solver:
MV9_endsec.m
MV9_helper.m
MV9.m
MV9_parallel.m

./HGmn-pump-sqrt2-scale-solver:
MV11.m

./LG00-pump-sqrt2-scale-solver:
lgmodesolve.m

