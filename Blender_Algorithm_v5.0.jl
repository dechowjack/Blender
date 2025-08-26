 # Constrained water and energy balance estimation 
# v4.2: adapted from v4.1 by treating SCF known, G estimated offline
# v4.3: adapted from v4.2 & 3.3 to handle to batch processing on OSC
# v4.4: adapted from v4.3 to update ipopt deprecation
# v4.5: moving σWRFG based on SCF 7/2/20 JLD
# v4.6: adding SDC to objective function: MD & JD: 9/28/20
# v4.7: adding "pseudo-valid prior" to starting poitns. remove SDC. MD: 1/7/21
# v4.8: tweak error parameters. MD & JD: 1/8/21
# v4.9 prior valid added to objective function. daily G capped in constraints. more outputs written 10/5/21
# v4.9.1 remove SWE from objective
# v4.9.2 added air temp
# v5.0  This is the final version used to generate data for my PhD dissertation (http://rave.ohiolink.edu/etdc/view?acc_num=osu1723774509328696) as well the data for
# the in preparation manuscript Dechow et al., 2025 (to be submitted to WRR)

## Add Julia Packages
using Pkg		# Package Manager
using JuMP		# Solver component
using Ipopt		# Solver
using DelimitedFiles	# Req to read csv files

## Assign DataDir
DataDir = ARGS[1]

# DataDir is a moving directory containing an argument given in the job file
# in the shell script the arg is increased by a step counter
# same arg moves the pix number directory in each seq run
# Example Data directory structure
# DataDir= "/users/projectDir/jldechow/TuolumneWY16/Pix3252" -> example of directory structure

# 0. read in data 
WRFSWE = readdlm(DataDir * "/WRFSWE.txt"); #[m]
WRFP   = readdlm(DataDir * "/WRFP.txt"); #[m]
WRFG   = readdlm(DataDir * "/WRFG.txt");
MSCF   = readdlm(DataDir * "/MODSCAG.txt");
AirT   = readdlm(DataDir * "/WRFT.txt");

AirT=AirT.-273.15 # Convert from C to K

## 1.  Define parameters
#  1.1 Experiment parameters
exp_dir = DataDir 
σWSWE=0.4 	# m
RelPUnc=0.3 	#[-]
σWRFG=15 	# W m^-2
σWMPmin=.001 	#minimum uncertainty in fluxes: 1 mm/day
nt=365		# force to run 365 days regardless of leap year
t=1:nt		# make day vector
SWEmin=1.0e-6 	#1/1000 mm
ρnew=100 	#density of new snow
Δt=86400 	#daily
Gmax=300 	# to prevent craziness. more than adequate for daily
Gmin=-300 	# to prevent craziness. more than adequate for daily

# NOT USED σSCF=.1 	# combined uncertainty in observed SCF and SDC # NOT USED
# NOT USED ρnew=100 	# density of new snow # 
# NOT USED z0=0.01 	# roughness length used in SDC 
# NOT USED mf=1.0 	# melt factor used in SDC 
# NOT USED ρWRF=350 	# note: previous versions used time-varying snow density but that is not currently exported 
 
# 1.1.1. Update σWRFG based on SCF

σWRFG=zeros(nt,1)
σWRFG_rel=0.5

for i=1:nt
  if MSCF[i]>0.1 && WRFSWE[i]>0.1  # if both are snow covered
    σWRFG[i]=abs(WRFG[i])*σWRFG_rel; 
  elseif MSCF[i]<0.1 && WRFSWE[i]<0.1 # if both not snowy
    σWRFG[i]=25; 
  else
    σWRFG[i]=500; #if they disagree, then don't use prior in cost function
  end
end

# 1.2 physical parameters
ρw=1000 	#density of water
Lf=0.334E6 	#Latent heat of fusion J/kg
cs_ice=2102 	#specific heat capacity for ice J/kg/K


# 1.3 Match up SWE and MSCF

for i=1:nt
  if MSCF[i]==0
    WRFSWE[i] = 0
  end
end

## 2. Compute useful variables and set up  arrays
#  2.1 Define SWEmax: upper limit of SWE, based on observed SCF
SWEmax=zeros(nt,1)
for i=1:nt
  if MSCF[i]==0
    SWEmax[i]=SWEmin
  else
    SWEmax[i]=5.0
  end
end

# NO LONGER USED
# NOT USED #2.2 SDC denominator: dependent only on snow density
# NOT USED DenomSDC=2.5.*z0.*(ρWRF./ρnew).^mf

# 2.3 Uncertainty for accumulation 
σWRFP=zeros(nt,1)
for i=1:nt
  if WRFP[i]==0.
    σWRFP[i]=σWMPmin
  else
    σWRFP[i]=WRFP[i]*RelPUnc
  end
end

## 3. Solve
#  3.1 Solve for G Prior Valid

Gmelt_pv=zeros(nt,1)
G_pv=zeros(nt,1)
U_pv=zeros(nt,1)
SWEpv=WRFSWE

for i=2:nt
  Gmelt_pv[i-1]=-(SWEpv[i]-SWEpv[i-1]-WRFP[i-1])*Lf*ρw/Δt
  if Gmelt_pv[i-1]<0.
    Gmelt_pv[i-1]=0
    SWEpv[i]=SWEpv[i-1]+WRFP[i-1]
  end
end
for i=2:nt
  if Gmelt_pv[i]>0. && MSCF[i] >0.
    G_pv[i]=Gmelt_pv[i]/MSCF[i]
    U_pv[i]=0.
  else
    G_pv[i]=WRFG[i]
    U_pv[i]=U_pv[i-1]+WRFP[i-1]*AirT[i-1]*ρw*cs_ice + G_pv[i-1]*MSCF[i-1]*Δt
    if U_pv[i]>0. || SWEpv[i]==0.
      U_pv[i]=0.
    end
  end
end

for i=2:nt
    if Gmelt_pv[i] > Gmax
        σWRFG[i]=1.0e9
    end
end

# 3.2 Solve for the posterior using prior valid
m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>5000))

@variable(m, SWEmin <= SWE[i=1:nt] <= SWEmax[i],start=SWEpv[i] )
@variable(m,  Precip[i=1:nt]>=0. ,start=WRFP[i])
@variable(m,  G[i=1:nt] , start=G_pv[i])
@variable(m, 0 <= Gmelt[i=1:nt] <= Gmax, start=Gmelt_pv[i])
@variable(m, Us[i=1:nt] <=0, start=U_pv[i])

@objective(m,Min,sum((Precip-WRFP).^2 ./σWRFP.^2)+ sum((Gmelt-Gmelt_pv).^2 ./σWRFG.^2) )
                   
for i in 1:nt-1
  @constraint(m,SWE[i+1]==SWE[i]+Precip[i]-Gmelt[i]*Δt/Lf/ρw)
  @NLconstraint(m,Us[i+1]==Us[i]+(1-(tanh(Us[i]/10000)+1))*G[i]*MSCF[i]*Δt+
                          Precip[i]*AirT[i]*ρw*cs_ice) #m x K x kg/m3 x J/kg/K
end
@constraint(m,Us[1]==0) 
@constraint(m,Us[nt]==0)
for i in 1:nt
  @NLconstraint(m,Gmelt[i]==G[i]*MSCF[i]*(tanh(Us[i]/10000)+1))
end
optimize!(m)

# 4. Output

SWEhat=JuMP.value.(SWE)
GmeltHat=JuMP.value.(Gmelt)
Ghat=JuMP.value.(G)
Ushat=JuMP.value.(Us)
Phat=JuMP.value.(Precip)

writedlm(exp_dir * "/SWE.txt",SWEhat)
writedlm(exp_dir * "/Gmelt.txt",GmeltHat)
writedlm(exp_dir * "/G.txt",Ghat)
writedlm(exp_dir * "/Precip.txt",Phat)
writedlm(exp_dir * "/Us.txt",Ushat)
writedlm(exp_dir * "/Gpv.txt",G_pv)
writedlm(exp_dir * "/Gmeltpv.txt",Gmelt_pv)
writedlm(exp_dir * "/Upv.txt",U_pv)
writedlm(exp_dir * "/SWEpv.txt",SWEpv)
