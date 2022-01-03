# DESCRIPTION
# This file contains the aircraft parameters. Always re-run DI_design.jl after
# modifying the aircraft parameters.

# ------------------------- MASS AND AERO PARAMETERS ---------------------------

# CG location
XCG=0.35
# reference CG location
XCGR=0.35
# wing area [ft²]
S=300.0
# wing span [ft]
B=30.0
# wing chord [ft]
CBAR=11.32
# weight [lb]
WEIGHT=20500.0
# momentum of engine turbine [lb-ft]
HE=160.0
# inertia constants
C1=-0.770
C2=0.02755
C3=1.055E-4
C4= 1.642E-6
C5=0.9604
C6=1.759E-2
C7=1.792E-5
C8=-0.7336
C9=1.587E-5
# longitudinal location of vertical accelerometer [ft]
XA=15.0

# ------------------------ MISCELLANEOUS CONSTANTS------------------------------

# radians to degrees
R2D=180.0/pi
# degrees to radians
D2R=pi/180.0
# gravitational acceleration [ft²]
G=32.17

# ------------------- TRIM AND LINEARIZATION PARAMETERS ------------------------

# state vector perturbations
DELXLIN=[0.1;0.001;0.001;0.001;0.001;0.001;0.001;0.001;0.001;0.1;0.1;0.1;0.1]
# control vector perturbations
DELCLIN=[0.1;0.1;0.1;0.1]
# trim target tolerance
TOL=1e-5*[1.0;D2R;D2R;D2R;D2R;D2R;D2R;D2R;D2R;1.0;1.0;1.0;1.0;0.1]
# trim variables: [VT α β ϕ θ ψ P Q R pow controls]
TRIMVARS=[1:9; 13; 14; 15; 16; 17]
# trim targets
TRIMTARG=[collect(1:13); 15]
# number of states
NSTATES=13
# number of controls
NCTRLS=4
# number of outputs
NOUT=6
