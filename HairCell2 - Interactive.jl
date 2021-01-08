"""
Based on Nieman 2011
Model of bullfrog saccular hair cell
Hodgkin-Huxely type
"""

using Unitful
using Makie
using MakieLayout
using Printf

const k_int = 112.0u"mmol/L" #Intracellular potassium conc. (millimolar)
const k_ext = 2.0u"mmol/L" #Extracellular potassium conc. (millimolar)
const F = 96485.3329u"s*A/mol" #Faraday constant
const R = 8.314u"J/(mol*K)" #universal gas constant
const T = 295.15u"K" #Temperature in Kalvin
const dt= 1.0e-5u"s" #timestep - theirs was 1e-5
const BI_time=2.0u"s" #Length of burn in
const Pulse_time=0.5u"s" #Time of current application for pulse

"""
    Structure modelling a bullfrog saccular hair cell
"""
struct Hair_Cell
    #state vector(s)
    x::Array{Any,1} #contains V, m_K1f, m_K1s, m_h, m_DRK, m_Ca
    BK::Array{Any,1} #contains h_BKT and state-probs of BK current (eg C2)
    PoMET::Array{Any,1} #contains open state prob. of mechanoceptors
    currents::Array{typeof(1.0u"nA"),1} #stores current current values
    #Parameters
    #
    #Reversal Potentials
    E_K::typeof(1.0u"mV")
    E_h::typeof(1.0u"mV") #for cation h-current
    E_Ca::typeof(1.0u"mV")
    E_L::typeof(1.0u"mV")
    E_MET::typeof(1.0u"mV")
    #Maximal Conductances
    g_K1::typeof(1.0u"nS") #g_K1 is one of the main control parameters
    g_h::typeof(1.0u"nS") #for cation h-current
    g_Ca::typeof(1.0u"nS")
    g_L::typeof(1.0u"nS")
    g_MET::typeof(1.0u"nS")
    #Max permeabilities of currents
    P_DRK::typeof(1.0u"L/s")#"L/s")
    P_BKS::typeof(1.0u"L/s") #calcium current
    P_BKT::typeof(1.0u"L/s") #calcium current
    P_A::typeof(1.0u"L/s")
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b::Float64
    #cell capacitance
    Cm::typeof(1.0u"pF")
end

###################################################################
#Construct GUI
###################################################################


scene, layout=layoutscene(nrows=1,ncols=3,resolution = (1800,850))
#set_window_config!(vsync = false) #something github said might help it not randomly close

#set up axes for burn-in plot
layout[1,1]=BIplots_layout=GridLayout(2,1)
BIplots_layout[1,1]=BI_axis=LAxis(scene, title="Burn-In", ylabel="memb. pot. (mV)", xlabel="time (sec)",xzoomlock=true)
#set up axes for burn-in current plot
BIplots_layout[2,1]=BI_c_axis=LAxis(scene, title="Currents", ylabel="amplitude (nA)", xlabel="time (sec)",xzoomlock=true)

#set up axes for pulse response plot
layout[1,3]=PRplots_layout=GridLayout(2,1)
PRplots_layout[1,1]=PR_axis=LAxis(scene, title="Pulse Response", ylabel="memb. pot. (mV)", xlabel="time (sec)",xzoomlock=true)
#axis for currents during pulse
PRplots_layout[2,1]=PR_c_axis=LAxis(scene, title="Currents", ylabel="amplitude (nA)", xlabel="time (sec)",xzoomlock=true)

#slider panel
layout[1,2]=slider_panel=GridLayout(11,1)
#EK slider
slider_panel[1,1]=EK_layout=GridLayout(1,3)
EK_layout[1,1]=LText(scene,"E_K", textsize=20)
EK_layout[1,2]=E_K_slider=LSlider(scene, range=-100.0:0.2:-70.0,startvalue=-95.0)
EK_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_K_slider.value))

#Eh slider
slider_panel[2,1]=Eh_layout=GridLayout(1,3)
Eh_layout[1,1]=LText(scene,"E_h", textsize=20)
Eh_layout[1,2]=E_h_slider=LSlider(scene, range=-50.0:0.1:-35.0,startvalue=-45.0)
Eh_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_h_slider.value))


#E_Ca slider
slider_panel[3,1]=ECa_layout=GridLayout(1,3)
ECa_layout[1,1]=LText(scene,"E_Ca", textsize=20)
ECa_layout[1,2]=E_Ca_slider=LSlider(scene, range=40.0:0.5:100.0,startvalue=42.5)
ECa_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_Ca_slider.value))

#EL slider
slider_panel[4,1]=EL_layout=GridLayout(1,3)
EL_layout[1,1]=LText(scene,"E_L", textsize=20)
EL_layout[1,2]=E_L_slider=LSlider(scene, range=-70.0:0.0,startvalue=-40.0)
EL_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), E_L_slider.value))

#gK1 slider
slider_panel[5,1]=gK1_layout=GridLayout(1,3)
gK1_layout[1,1]=LText(scene,"g_K1", textsize=20)
gK1_layout[1,2]=g_K1_slider=LSlider(scene, range=3.0:0.5:42.0,startvalue=20.0)
gK1_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_K1_slider.value))

#gh slider
slider_panel[6,1]=gh_layout=GridLayout(1,3)
gh_layout[1,1]=LText(scene,"g_h", textsize=20)
gh_layout[1,2]=g_h_slider=LSlider(scene, range=1.0:0.05:6.0,startvalue=2.2)
gh_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_h_slider.value))

#gCa slider
slider_panel[7,1]=gCa_layout=GridLayout(1,3)
gCa_layout[1,1]=LText(scene,"g_Ca", textsize=20)
gCa_layout[1,2]=g_Ca_slider=LSlider(scene, range=1.0:0.01:5.0,startvalue=1.2)
gCa_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_Ca_slider.value))

#gL slider
slider_panel[8,1]=gL_layout=GridLayout(1,3)
gL_layout[1,1]=LText(scene,"g_L", textsize=20)
gL_layout[1,2]=g_L_slider=LSlider(scene, range=0.1:0.01:0.9,startvalue=0.77)
gL_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), g_L_slider.value))

#b slider
slider_panel[9,1]=b_layout=GridLayout(1,3)
b_layout[1,1]=LText(scene,"b (BK)", textsize=20)
b_layout[1,2]=b_slider=LSlider(scene, range=0.01:0.01:1.0,startvalue=0.5)
b_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), b_slider.value))

#v0 slider
slider_panel[10,1]=v0_layout=GridLayout(1,3)
v0_layout[1,1]=LText(scene,"v0", textsize=20)
v0_layout[1,2]=v0_slider=LSlider(scene, range=-120.0:40.0,startvalue=-50.0)
v0_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), v0_slider.value))

#pulse size slider
slider_panel[11,1]=p_layout=GridLayout(1,3)
p_layout[1,1]=LText(scene,"Pulse", textsize=20)
p_layout[1,2]=p_slider=LSlider(scene, range=0.005:0.005:0.1,startvalue=0.005)
p_layout[1,3]=LText(scene, lift(x->@sprintf("%.2f",x), p_slider.value))


################################################################################
#Done constructing GUI
################################################################################

"""
    Contructor for Hair Cell, gets passed parameter values
"""
function Hair_Cell(V::typeof(1.0u"mV"),
                    #Reversal potentials
                    E_Kp::Float64, E_hp::Float64, E_Cap::Float64, E_Lp::Float64,
                    #Maximal Conductances
                    g_K1p::Float64, g_hp::Float64, g_Cap::Float64, g_Lp::Float64,
                    bp::Float64)
    #Reversal Potentials
    E_K=E_Kp*1.0u"mV"
    E_h=E_hp*1.0u"mV"
    E_Ca=E_Cap*1.0u"mV"
    E_L=E_Lp*1.0u"mV"
    E_MET=0.0u"mV"
    #Maximal Conductances
    g_K1=g_K1p*1.0u"nS"
    g_h=g_hp*1.0u"nS"
    g_Ca=g_Cap*1.0u"nS"
    g_L=g_Lp*1.0u"nS"
    g_MET=0.65u"nS"
    #Max permeabilities of currents
    P_DRK=2.4e-14u"L/s"
    P_BKS=2e-13u"L/s" #calcium current
    P_BKT=14e-13u"L/s" #calcium current
    P_A=1.08e-13u"L/s"
    #b is another key control parameter (dimensionless)
    #b controls the strength of the above currents (BKS and BKT)
    b=bp
    #cell capacitance
    Cm=10.0u"pF" # (a,c)
    #initial m values for state vector, set as equilibrium for V
    m_K1f=m_k1_inf(V)
    m_K1s=m_k1_inf(V)
    m_h=m_h_inf(V)
    m_DRK=m_DRK_inf(V)
    m_Ca=m_Ca_inf(V)
    m_A = m_A_inf(V)
    hA1 = h_A_inf(V)
    hA2 = h_A_inf(V)
    #initial values for params of BK-current stuff (another state vector)
    h_BKT=h_BKT_inf(V) #h for inactivation
    C0=0.6
    C1=0.1
    C2=0.1
    O2=0.1
    O3=0.1
    ca_conc=0.0u"μmol/L"

    PoMet=0.15
    return Hair_Cell([V,m_K1f, m_K1s, m_h, m_DRK, m_Ca, m_A, hA1, hA2], [h_BKT,C0,C1,C2,O2,O3,ca_conc],
    [PoMet],[0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA",0.0u"nA"],
    E_K,E_h,E_Ca,E_L,E_MET,g_K1,g_h,g_Ca,g_L,g_MET,P_DRK,P_BKS,P_BKT,P_A,b,Cm)
end

"""
function taking a hair cell and a current, mutates hair cell based on effect of current
"""
function update(cell::Hair_Cell,input=0.0u"nA")
    V=cell.x[1]
    thing=(V*F)/(R*T)
    pow=uconvert(Unitful.NoUnits,thing)
    IK1 = uconvert(u"nA",cell.g_K1*(V-cell.E_K)*(0.7*mk1f(V,cell.x[2])+0.3*mk1s(V,cell.x[3])))
    cell.currents[1]=IK1
    Ih = cell.g_h*(V-cell.E_h)*(3*mh(V,cell.x[4])^2*(1-mh(V,cell.x[4]))+mh(V,cell.x[4])^3)
    cell.currents[2]=Ih
    IDRK = uconvert(u"nA",cell.P_DRK*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(mDRK(V,cell.x[5])^2))
    cell.currents[3]=IDRK
    ICa = cell.g_Ca*(V-cell.E_Ca)*mCa(V,cell.x[6])^3
    cell.currents[4]=ICa
    IBKS = cell.b*cell.P_BKS*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(O_2(V,cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[7])+O_3(V,cell.BK[5],cell.BK[6],cell.BK[7]))
    cell.currents[5]=IBKS
    IBKT = cell.b*cell.P_BKT*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*(O_2(V,cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[7])+O_3(V,cell.BK[5],cell.BK[6],cell.BK[7]))*hBKT(V,cell.BK[1])
    cell.currents[6]=IBKT
    IL = cell.g_L*(V-cell.E_L)
    cell.currents[7]=IL
    IMET = cell.g_MET*cell.PoMET[1]*(V-cell.E_MET)
    cell.currents[8]=IMET
    aIA=a1(V)
    IA=cell.P_A*((V*F^2)/(R*T))*((k_int-k_ext*exp(-pow))/(1-exp(-pow)))*((mA(V,cell.x[7]))^3)*(aIA*hA1(V,cell.x[8])+(1-aIA)hA2(V,cell.x[9]))
    cell.currents[9]=IA
    ΔV=uconvert(u"mV",(-IK1-Ih-IDRK-ICa-IBKS-IBKT-IL-IMET-IA+input)*dt/cell.Cm)

    #update things
    cell.x[1]= V + ΔV #voltage/membrane potential
    cell.x[2]= mk1f(V,cell.x[2]) #m_K1f
    cell.x[3]= mk1s(V,cell.x[3]) #m_K1s
    cell.x[4]=   mh(V,cell.x[4]) #m_h
    cell.x[5]= mDRK(V,cell.x[5]) #m_DRK
    cell.x[6]=  mCa(V,cell.x[6]) #m_Ca
    cell.x[7]=   mA(V,cell.x[7]) #mA
    cell.x[8]=  hA1(V,cell.x[8]) #hA1
    cell.x[9]=  hA2(V,cell.x[9]) #hA2
    cell.BK[1]= hBKT(V,cell.BK[1]) #h_BKT
    cell.BK[2]= C_0(cell.BK[3],cell.BK[4],cell.BK[5],cell.BK[6]) #C0
    cell.BK[3]= C_1(V,cell.BK[2],cell.BK[3],cell.BK[4],cell.BK[7]) #C1
    cell.BK[4]= C_2(V,cell.BK[3],cell.BK[4],cell.BK[5],cell.BK[7]) #C2
    cell.BK[5]= O_2(V,cell.BK[4],cell.BK[5],cell.BK[6],cell.BK[7]) #02
    cell.BK[6]= O_3(V,cell.BK[5],cell.BK[6],cell.BK[7]) #03
    cell.BK[7]= Ca_conc(cell.BK[7],ICa)#calcium
end

#equations for K1 current
mk1f(V::typeof(1.0u"mV"),mk1f)=mk1f+Δm_k1f(V,mk1f)
Δm_k1f(V::typeof(1.0u"mV"),mk1f)=(m_k1_inf(V)-mk1f)*dt/τ_k1f(V)
τ_k1f(V::typeof(1.0u"mV"))=0.7u"ms"*exp(-(V + 120u"mV")/43.8u"mV")+0.04u"ms"

mk1s(V::typeof(1.0u"mV"),mk1s)=mk1s+Δm_k1s(V,mk1s)
Δm_k1s(V::typeof(1.0u"mV"),mk1s)=(m_k1_inf(V)-mk1s)*dt/τ_k1s(V)
τ_k1s(V::typeof(1.0u"mV"))=14.1u"ms"*exp(-(V + 120u"mV")/28u"mV")+0.04u"ms"

m_k1_inf(V::typeof(1.0u"mV"))=(1+exp((V + 110u"mV")/11u"mV"))^-1

#equations for h cation current
mh(V,mh)=mh+Δm_h(V,mh)
Δm_h(V,mh)=(m_h_inf(V)-mh)*dt/τ_h(V)
m_h_inf(V)=(1+exp((V+87u"mV")/16.7u"mV"))^(-1)
τ_h(V)=63.7u"ms"+135.7u"ms"*exp(-((V+91.4u"mV")/21.2u"mV")^2)

#equations for DRK current
mDRK(V,mDRK)= mDRK+ΔmDRK(V,mDRK)
ΔmDRK(V, mDRK) = (m_DRK_inf(V)- mDRK)*(dt/τ_DRK(V))
m_DRK_inf(V) = (1+exp((V+48.3u"mV")/4.19u"mV"))^(-1/2)

τ_DRK(V) = (α_DRK(V)+β_DRK(V))^(-1)
α_DRK(V) = (3.2u"ms"*exp(-V/20.9u"mV")+3.0u"ms")^(-1)
β_DRK(V) = (1467.0u"ms"*exp(V/5.96u"mV")+9.0u"ms")^(-1)

#equations for ICa current
mCa(V,mCa) = mCa+Δm_Ca(V,mCa)
Δm_Ca(V, mCa) = (m_Ca_inf(V)-mCa)*(dt/τ_Ca(V))
m_Ca_inf(V) = (1+exp(-(V+55u"mV")/12.2u"mV"))^(-1)
τ_Ca(V) = 0.046u"ms" + 0.325u"ms"*(exp(-((V+77u"mV")/51.67u"mV")^2))

#equations for BK currents:
#Params:
#Used directly in state-probability equations
const k_1=300.0u"1/s" #s^-1
const k_2=5000.0u"1/s" #s^-1
const k_3=1500.0u"1/s" #s^-1
const β_c=2500.0u"1/s"#1000.0u"1/s" #s^-1 #changed to OG value *******************************

#for auxiliary equations
const α_c_0=450.0u"1/s" #s^-1
const K1_0=6.0u"μmol/L" #micromolar
const K2_0=45.0u"μmol/L" #micromolar
const K3_0=20.0u"μmol/L" #micromolar
const z = 2 #sign of charge of Ca2+
const δ1 = 0.2 #fraction of the electric field experienced by Ca2 at the 1st binding site
const δ2 = 0.0 #fraction of the electric field experienced by Ca2 at the 2nd binding site
const δ3 = 0.2 #fraction of the electric field experienced by Ca2 at the 3rd binding site
const V_A = 30.0u"mV" #potential used to express the voltage dependence of ac

α_c(V)=α_c_0*exp(-V/V_A)

K1(V)=K1_0*exp(δ1*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))
K2(V)=K2_0*exp(δ2*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))
K3(V)=K3_0*exp(δ3*z*uconvert(Unitful.NoUnits,(V*F)/(R*T)))

k1(V)=k_1/(K1(V))
k2(V)=k_2/(K2(V))
k3(V)=k_3/(K3(V))

C_0(C1,C2,O2,O3) = 1.0 - (C1 + C2 + O2 + O3)
C_1(V,C0,C1,C2,Ca)= C1 + ΔC_1(V,C0,C1,C2,Ca)
ΔC_1(V,C0,C1,C2,Ca)=uconvert(Unitful.NoUnits,(k1(V)*Ca*C0+k_2*C2-(k_1+k2(V)*Ca)*C1)*dt)
C_2(V,C1,C2,O2,Ca)= C2 + ΔC_2(V,C1,C2,O2,Ca)
ΔC_2(V,C1,C2,O2,Ca)=uconvert(Unitful.NoUnits,(k2(V)*Ca*C1+α_c(V)*O2-(k_2+β_c)*C2)*dt)
O_2(V,C2,O2,O3,Ca)= O2 + ΔO_2(V,C2,O2,O3,Ca)
ΔO_2(V,C2,O2,O3,Ca)=uconvert(Unitful.NoUnits,(β_c*C2+k_3*O3-(α_c(V)+k3(V)*Ca)*O2)*dt)
O_3(V,O2,O3,Ca)= O3+ΔO_3(V,O2,O3,Ca)
ΔO_3(V,O2,O3,Ca)=uconvert(Unitful.NoUnits,(k3(V)*Ca*O2-k_3*O3)*dt)

Ca_conc(Ca,ICa)=Ca+ΔCa_conc(Ca,ICa)
const U=0.02 #free calcium proportion
const Cvol = 1256.637u"μm^3"#1.25u"pL"#cell volume
const ξ = 3.4e-5#proportion of cell hoarding calcium
const Ksca = 2800u"1/s" #rate constant

ΔCa_conc(Ca,ICa)=uconvert(u"μmol/L",(-(U*ICa)/(z*F*Cvol*ξ)-Ksca*Ca)*dt)


#inactivation of BKT current
hBKT(V,hBKT) = hBKT+Δh_BKT(V,hBKT)
Δh_BKT(V,hBKT)= (h_BKT_inf(V)-hBKT)*(dt/τ_BKT(V))
h_BKT_inf(V) = (1+exp((V+61.6u"mV")/3.65u"mV"))^(-1)
τ_BKT(V) = 2.1u"ms"+9.4u"ms"*exp(-((V+66.9u"mV")/17.7u"mV")^2)

#Equations for IA current
mA(V,mA)=mA+Δm_A(V,mA)
Δm_A(V,mA)=(m_A_inf(V)-mA)*dt/τ_A(V)
m_A_inf(V)=(1+exp(-(V+61u"mV")/10.7u"mV"))^(-1)
function τ_A(V)
    if ustrip(V)<-55
        a=173u"ms"
        b=17.9u"mV"
        c=5.4u"ms"
    else
        a=0.48u"ms"
        b=-23.1u"mV"
        c=1.14u"ms"
    end
    return a*exp(V/b)+c
end

#Inactivation of A-current
hA1(V,hA1) = hA1+Δh_A1(V,hA1)
Δh_A1(V,hA1)= (h_A_inf(V)-hA1)*(dt/τ_hA1(V))

hA2(V,hA2) = hA2+Δh_A2(V,hA2)
Δh_A2(V,hA2)= (h_A_inf(V)-hA2)*(dt/300u"ms")

h_A_inf(V) = (1+exp((V+83u"mV")/3.9u"mV"))^(-1)
τ_hA1(V) = 74u"ms"+321u"ms"*exp(-((V+82u"mV")/15.5u"mV")^2)

a1(V)=(1-0.54)/(1+exp((V+21.5u"mV")/15.6u"mV")) + 0.54

"""
    runs for a given time to let parameters settle into equilibrium state
"""
function burnin(cell, time=5.0u"s")
    t=0.0u"s":dt:time
    voltBI=Array{typeof(1.0u"mV")}(undef, length(t))
    currents_p=Array{typeof(1.0u"nA"),2}(undef,length(t),9)
    for i in 1:length(t)
        update(cell)
        voltBI[i]=cell.x[1]
        currents_p[i,:].=cell.currents
    end
    return ustrip.(voltBI),ustrip.(t),ustrip.(currents_p)
end

"""
   Generates a pulse waveform same length as t - called by pulseReponse
"""
function pulse(time=0.5u"s", start=2.0e-1u"s", len=2.0e-1u"s", amplitude=p_slider.value[]*1.0u"nA")
    t=0.0u"s":dt:time
    u = Array{typeof(1.0u"nA")}(undef, length(t))
    u.=0.0u"nA"
    u[findall( t-> (t>=start) & ( t<start+len), t)] .= amplitude

    return u
end

"""
    Models response of cell to a pulse of given amplitude and length
"""
function pulseResponse(cell,amplitude=0.005u"nA", time=Pulse_time, start=2.0e-1u"s", len=2.0e-1u"s")
    t=0.0u"s":dt:time
    input = pulse(time,start,len,amplitude) #in nA
    voltages=Array{typeof(1.0u"mV")}(undef, length(t))
    currents_p=Array{typeof(1.0u"nA"),2}(undef,length(t),9)

    for i in 1:length(input)
        update(cell,input[i])
        voltages[i]=cell.x[1]
        currents_p[i,:].=cell.currents
    end
    return voltages,t,currents_p
end

"""
  Calls other functions to allow animation
"""
function go(
            E_K::Float64, #E_K_slider.value[]::Float64,
            E_h::Float64, #E_h_slider.value[]::Float64,
            E_Ca::Float64, #E_Ca_slider.value[]::Float64,
            E_L::Float64, #E_L_slider.value[]::Float64,
            #Maximal Conductances
            g_K1::Float64,#g_K1_slider.value[]::Float64,
            g_h::Float64, #g_h_slider.value[]::Float64,
            g_Ca::Float64, #g_Ca_slider.value[]::Float64,
            g_L::Float64, #g_L_slider.value[]::Float64,
            #Other sliders
            b::Float64, #b_slider.value[]::Float64,
            v0::Float64, #v0_slider.value[]::Float64,
            amp::Float64, #p_slider.value[]::Float64,
            #Other params
            bi_time=BI_time
            )
    #initiate hair cell
    helga=Hair_Cell(v0*1.0u"mV",E_K, E_h,
                    E_Ca, E_L, g_K1, g_h,
                    g_Ca, g_L, b
                    )
    #burn in
    burn_arrays=burnin(helga,bi_time)
    burn_time=burn_arrays[2]
    burn_v=burn_arrays[1]
    burn_currents=burn_arrays[3]
    #Injected current reponse
    pulse_arrays=pulseResponse(helga,amp*1.0u"nA")
    pulse_V=ustrip.(pulse_arrays[1])
    pulse_time=ustrip.(pulse_arrays[2])
    pulse_currents=ustrip.(pulse_arrays[3])
    return burn_v,burn_time,pulse_V,pulse_time,burn_currents,pulse_currents
end

const all_arrays=go(E_K_slider.value[],E_h_slider.value[],
                      E_Ca_slider.value[], E_L_slider.value[],
                      g_K1_slider.value[], g_h_slider.value[],
                      g_Ca_slider.value[], g_L_slider.value[],
                      b_slider.value[], v0_slider.value[],
                      p_slider.value[])
const BI_volt0=all_arrays[1]
const BI_time0=all_arrays[2]
const PR_volt0=all_arrays[3]
const PR_time0=all_arrays[4]
const BI_currents0=all_arrays[5]
const PR_currents0=all_arrays[6]

#################################################################################
#Plot graphs based on initial slider values
#################################################################################
const BI_plothandle=lines!(BI_axis,BI_time0,BI_volt0, color=:black)

const BI_c_plothandle1=lines!(BI_c_axis,BI_time0,BI_currents0[:,1],color=:black)#   IK1
const BI_c_plothandle2=lines!(BI_c_axis,BI_time0,BI_currents0[:,2],color=:red)#     Ih
const BI_c_plothandle3=lines!(BI_c_axis,BI_time0,BI_currents0[:,3],color=:orange)#    IDRK
const BI_c_plothandle4=lines!(BI_c_axis,BI_time0,BI_currents0[:,4],color=:blue)#  ICa
const BI_c_plothandle5=lines!(BI_c_axis,BI_time0,BI_currents0[:,5],color=:green)# IBKS
const BI_c_plothandle6=lines!(BI_c_axis,BI_time0,BI_currents0[:,6],color=:green,linestyle=:dash)#   IBKT
const BI_c_plothandle7=lines!(BI_c_axis,BI_time0,BI_currents0[:,7],color=:red,linestyle=:dash)# IL
const BI_c_plothandle8=lines!(BI_c_axis,BI_time0,BI_currents0[:,8],color=:blue,linestyle=:dash)#  IMET
const BI_c_plothandle9=lines!(BI_c_axis,BI_time0,BI_currents0[:,9],color=:orange,linestyle=:dash)#  IA

const PR_plothandle=lines!(PR_axis,PR_time0,PR_volt0)

const PR_c_plothandle1=lines!(PR_c_axis,PR_time0,PR_currents0[:,1],color=:black)#   IK1
const PR_c_plothandle2=lines!(PR_c_axis,PR_time0,PR_currents0[:,2],color=:red)#     Ih
const PR_c_plothandle3=lines!(PR_c_axis,PR_time0,PR_currents0[:,3],color=:orange)#    IDRK
const PR_c_plothandle4=lines!(PR_c_axis,PR_time0,PR_currents0[:,4],color=:blue)#  ICa
const PR_c_plothandle5=lines!(PR_c_axis,PR_time0,PR_currents0[:,5],color=:green)# IBKS
const PR_c_plothandle6=lines!(PR_c_axis,PR_time0,PR_currents0[:,6],color=:green,linestyle=:dash)#   IBKT
const PR_c_plothandle7=lines!(PR_c_axis,PR_time0,PR_currents0[:,7],color=:red, linestyle=:dash)# IL
const PR_c_plothandle8=lines!(PR_c_axis,PR_time0,PR_currents0[:,8],color=:blue, linestyle=:dash)#  IMET
const PR_c_plothandle9=lines!(PR_c_axis,PR_time0,PR_currents0[:,9],color=:orange,linestyle=:dash)#  IA

# Node attached to all sliders - make a thing called stuff which stores the
# return values of go() when any of the sliders are altered
const stuff=lift(go, E_K_slider.value, E_h_slider.value,
                  E_Ca_slider.value, E_L_slider.value,
                  g_K1_slider.value, g_h_slider.value,
                  g_Ca_slider.value, g_L_slider.value,
                  b_slider.value, v0_slider.value,
                  p_slider.value)

# When a thing called stuff is made, do all the stuff in the lift block
@lift begin
    arrays=$stuff
    # Extract relevant arrays
    BIV=arrays[1]
    PRV=arrays[3]
    BIc=arrays[5]
    PRc=arrays[6]
    # Replace the y-value arrays in the plot handles with new ones
    BI_plothandle[2]=BIV
    PR_plothandle[2]=PRV
    BI_c_plothandle1[2]=BIc[:,1]
    BI_c_plothandle2[2]=BIc[:,2]
    BI_c_plothandle3[2]=BIc[:,3]
    BI_c_plothandle4[2]=BIc[:,4]
    BI_c_plothandle5[2]=BIc[:,5]
    BI_c_plothandle6[2]=BIc[:,6]
    BI_c_plothandle7[2]=BIc[:,7]
    BI_c_plothandle8[2]=BIc[:,8]
    BI_c_plothandle9[2]=BIc[:,9]
    PR_c_plothandle1[2]=PRc[:,1]
    PR_c_plothandle2[2]=PRc[:,2]
    PR_c_plothandle3[2]=PRc[:,3]
    PR_c_plothandle4[2]=PRc[:,4]
    PR_c_plothandle5[2]=PRc[:,5]
    PR_c_plothandle6[2]=PRc[:,6]
    PR_c_plothandle7[2]=PRc[:,7]
    PR_c_plothandle8[2]=PRc[:,8]
    PR_c_plothandle9[2]=PRc[:,9]
end

scene
