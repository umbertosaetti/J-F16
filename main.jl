# TO RUN:
# include("main.jl")
#
# AUTHOR
# Dr. Umberto Saetti, Assistant Professor, Department of Aerospace Engineering,
# Auburn university
# saetti@auburn.edu
#
# DATE
# 1/3/2022
#
# DESCRIPTION
# Nonlinear simulation model of the rigid-body dynamics of an F-16
# aircraft. The simulation model is based on Refs 1 and 2. The trim,
# linearization routines, and overall architecture of the simulation
# follows the teachings of Dr. Joe Horn at Penn State. When using, please cite
# Ref. 3.
#
# REFERENCES
# 1) Stevens, B. L., Lewis, F. L., and Johnson E. N., “Aircraft Control
#    and Simulation: Dynamics, Controls, Design, and Autonomous Systems,”
#    Wiley, 3rd Edition, 2015.
#    ISBN: 978-1-118-87098-3
# 2) Garza, F. R., and Morelli, E. A., "A Collection of Nonlinear Aircraft
#    Simulations in MATLAB," NASA TM-2003-212145, 2003.
# 3) Saetti, U., and Horn, J. F., "Flight Simulation and Control using the Julia
#    Language", AIAA Scitech Forum, San Diego, CA, Jan 3-7, 2022.
#    DOI: https://arc.aiaa.org/doi/10.2514/6.2022-2354.
#
# STATES             | NINDICES    | DESCRIPTION
# _________________________________________________________________________
# V_T                | 1           | airspeed [ft/s]
# alpha              | 2           | angle of attach [rad]
# beta               | 3           | sideslip [rad]
# phi theta psi      | 4  5  6     | Euler angles [rad]
# p q r              | 7  8  9     | angula rates [rad]
# x                  | 10          | north position [ft]
# y                  | 11          | east position [ft]
# z                  | 12          | altitude [ft]
# P_a                | 13          | power [lb-ft^2/s]
#
# CONTROL INPUTS     | INDICED     | DESCRIPTION     | RANGE
# _________________________________________________________________________
# delta_t            | 1           | throttle [0-1]  |  (0, 1)
# delta_a            | 2           | aileron [deg]   |  (-21.5, 21.5) deg
# delta_e            | 3           | elevator [deg]  |  (-25, 25) deg
# delta_r            | 4           | rudder [deg]    |  (-30, 30) deg
#
# -------------------------------------------------------------------------

# include packages
#using DifferentialEquations
using LinearAlgebra
using SparseArrays
using ControlSystems
using MAT
using Interpolations
using ApproxFun
using Printf
#ENV["MPLBACKEND"]="tkagg"
using Plots
using LaTeXStrings
using JLD
#using PyPlot
#using GR

@printf("\n                     J-F16           \n")
@printf("\n          Author: Dr. Umberto Saetti\n")
@printf("\n-----------------------------------------------\n")

# include functions
include("atan2.jl")
include("linearize.jl")
include("trimmer.jl")
#include("finp_F16_loop.jl")
include("THRUST.jl")
include("DAMP.jl")
include("DNDR.jl")
include("DNDA.jl")
include("DLDR.jl")
include("DLDA.jl")
include("CN.jl")
include("CL.jl")
include("CM.jl")
include("CZ.jl")
include("CY.jl")
include("CX.jl")
include("RTAU.jl")
include("TGEAR.jl")
include("PDOT.jl")
include("ADC.jl")
include("fix.jl")
include("F16.jl")
include("rk4.jl")
include("finp.jl")
include("simulate.jl")
include("simulate_DI.jl")
include("F16_DI.jl")
# load F-16 constants
include("F16_constants.jl")

# ------------------------ TRIM FLIGHT DYNAMICS MODEL --------------------------

# set velocity [ft/s]
speed=600.0
# vertical flight path angle [deg]
gamma=0.0
# horizontal flight path angle [deg]
chi=0.0
# set initial altitude [ft]
alt=2000.0
# set initial North-South and East-West position [ft]
p_north=0.0
p_east=0.0
# turn rate [deg/s]
turnrate=0.0
# pull-up rate [deg/s]
pullup=0.0
# roll rate [deg/s]
rollrate=0.0
# initial guess for state and control vectors
phi_est=turnrate*D2R*speed/32.17
pndot=speed*cos(gamma*D2R)*cos(chi*D2R)
pedot=speed*cos(gamma*D2R)*sin(chi*D2R)
hdot=speed*sin(gamma*D2R)
x0=[speed; 0.0; 0.0; phi_est; 0; chi*D2R; 0.0; 0.0; 0.0; p_north; p_east; alt;
    0.0]
u0=[0.2; 0; 0; 0]
targ_des=[0; 0; 0; rollrate*D2R; pullup*D2R; turnrate*D2R; 0; 0; 0; pndot;
          pedot; hdot; 0; 0]

# trim aircraft
x0, u0, iter = trimmer(F16!,x0,u0,targ_des)
# linearize
A, Bmat, C, D = linearize(F16!,x0,u0)
# partition in longitudinal and lateral flight dynamics
Alon=A[[1,2,5,8],[1,2,5,8]]
Alat=A[[3,4,7,9],[3,4,7,9]]

# ----------------------------- SPECTRAL ANALYSIS ------------------------------

# eigenvalues
eigs=eigvals(A)
eigsLon=eigvals(Alon)
eigsLat=eigvals(Alat)
# store eigenvalues
save("eigs_J-SimpleHel.jld","eigs",eigs)
# plot eigenvalues
gr(size=(800,600))
plt=plot(real(eigs),imag(eigs),seriestype=:scatter,
      label="", legendfontsize=12, legend=:topleft,
      marker = (:circle, :royalblue, 6), markerstrokecolor = :royalblue)
xaxis!("Real",xguidefontsize=14)
yaxis!("Imag",yguidefontsize=14)
display(plt)
savefig(".\\Plots\\eigs_600fts.png")
savefig(".\\Plots\\eigs_600fts.svg")

# ------------------------------ TIME SIMULATION -------------------------------

# time step [s]
dt=0.01
# length of simulation [s]
Tsim=8
# run open-loop simulation
state_OL, time_OL = simulate(F16!,finp,Tsim,dt,x0,u0)
# number of dynamic inverse states
NDISTATES=9
# initial state vector of closed-loop simulation
x0_DI=[x0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; u0[2]; u0[3]; u0[4]]
# run closed-loop simulation
state_CL, time_CL = simulate_DI(F16_DI!,finp,Tsim,dt,x0,u0)

# plot attitude
gr(size=(800,600))
i=4
p1=plot(time_OL,state_OL[i,:]*180/pi,label="Open Loop",legend=:bottomright, legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,-state_CL[i,:]*180/pi,label="Closed Loop",legend=:bottomright, legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("φ [deg]",yguidefontsize=14)
xlims!((0,8))
i=5
p2=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",legend=:bottomright, legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("θ [deg]",yguidefontsize=14)
xlims!((0,8))
i=6
p3=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,-state_CL[i,:]*180/pi,label="",legend=:bottomright, legendfontsize=12,line=(3, :brown3, :dash))
xaxis!("Time [s]",xguidefontsize=14)
yaxis!("ψ [deg]",yguidefontsize=14)
xlims!((0,8))
p4=plot(p1, p2, p3, layout=(3,1), leftmargin=3Plots.mm)
display(p4)
savefig(".\\Plots\\att_F16.svg")

# plot angular rates
gr(size=(800,600))
i=7
p5=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,-state_CL[i,:],label="",line=(3, :brown3, :dash))
yaxis!("p [rad/s]",yguidefontsize=14)
xlims!((0,8))
i=8
p6=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",legend=:bottomright, legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("q [rad/s]",yguidefontsize=14)
xlims!((0,8))
i=9
p7=plot(time_OL,state_OL[i,:],label="Open Loop",legend=:bottomright, legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,-state_CL[i,:],label="Closed Loop",line=(3, :brown3, :dash))
xaxis!("Time [s]",xguidefontsize=14)
yaxis!("r [rad/s]",yguidefontsize=14)
xlims!((0,8))
p8=plot(p5, p6, p7, layout=(3,1), leftmargin=3Plots.mm)
display(p8)
savefig(".\\Plots\\ang_F16.svg")
