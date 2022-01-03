# To run:
# include("DI_design.jl")

# DESCRIPTION
# Trims and linearizes the aircraft model at incremental flight speeds, and
# generates the gains of the dynamic inversion controller. Always run after
# modifying the aircraft parameters in the F16_constants.jl file.
#
# ------------------------------------------------------------------------------

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
# load F-16 constants
include("F16_constants.jl")

# vector of trim speeds [kts]
speed_vec=range(200,stop=800,length=7)
# initialize DI matrices
CATAB = zeros(3,3,length(speed_vec))
CBinvTAB = zeros(3,3,length(speed_vec))
x0_mat = zeros(NSTATES,length(speed_vec))
u0_mat = zeros(NCTRLS,length(speed_vec))
# trim at different speeds
for iv=1:length(speed_vec)
    # set velocity [ft/s]
    speed=speed_vec[iv]*1.688
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
    x0, u0, iter = trimmer(F16,x0,u0,targ_des)
    # store trim state and controls
    x0_mat[:,iv] = x0
    u0_mat[:,iv] = u0
    # linearized aircraft dynamics
    A, Bmat, C, D = linearize(F16,x0,u0)
    # dynamic inverse matrices
    C_DI=[1 0 0
          0 1 0
          0 0 1]
    CATAB[:,:,iv]=C_DI*A[7:9,7:9]
    CBinvTAB[:,:,iv]=inv(C_DI*Bmat[7:9,[2,3,4]])
end

# command filters time constants
tau_p=1/10
tau_q=1/10
tau_r=1/4
# washout filter for controls
tau_u=5.0;
# disturbance rejection frequency and damping and integrator poles
wnroll_d=10.0
dmproll_d=1.0
wnpitch_d=10.0
dmppitch_d=1.0
wnyaw_d=4.;
dmpyaw_d=1.;
# feedback gains
kp_roll=2*wnroll_d*dmproll_d
ki_roll=wnroll_d^2
kp_pitch=2*wnpitch_d*dmppitch_d
ki_pitch=wnpitch_d^2
kp_yaw=2*dmpyaw_d*wnyaw_d
ki_yaw=wnyaw_d^2;
# pilot gains
k_a=(290*pi/180)/21.5
k_e=(9-1)*G/25 ./(500*1.688)
k_r=1/30
# save DI matrices and gains into a file
save("DI.jld","CATAB",CATAB,"CBinvTAB",CBinvTAB,"kp_roll",kp_roll,"ki_roll",
    ki_roll,"kp_pitch",kp_pitch,"ki_pitch",ki_pitch,"kp_yaw",kp_yaw,"ki_yaw",
    ki_yaw,"k_a",k_a,"k_e",k_e,"k_r",k_r,"tau_p",tau_p,
    "tau_q",tau_q,"tau_r",tau_r,"tau_u",tau_u)

# plot trim pitch attitude and controls
gr(size=(800,600))
p10=plot(speed_vec,x0_mat[5,:]*180/pi,label="",line=(3, :royalblue, :solid), marker = (:circ, :royalblue, 8), markerstrokecolor = :royalblue)
xaxis!((200, 800), 200:50:800)
yaxis!("θ [deg]",yguidefontsize=14)
p11=plot(speed_vec,u0_mat[1,:]*100,label="",line=(3, :royalblue, :solid), marker = (:circ, :royalblue, 8), markerstrokecolor = :royalblue)
xaxis!((200, 800), 200:50:800)
yaxis!(L"\delta_t [\%]",yguidefontsize=14)
p12=plot(speed_vec,u0_mat[3,:],label="",line=(3, :royalblue, :solid), marker = (:circ, :royalblue, 8), markerstrokecolor = :royalblue)
xaxis!((200, 800), 200:50:800)
yaxis!(L"\delta_e [\mathrm{deg}]",yguidefontsize=14)
p13=plot(p10, p12, p11, layout=(3,1), leftmargin=3Plots.mm)
display(p13)
savefig(".\\Plots\\trim_F16.svg")
