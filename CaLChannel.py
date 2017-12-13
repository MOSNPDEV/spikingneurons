from brian2 import *
import brian2.numpy_ as np
prefs.codegen.target = 'numpy'

start_scope()
# Implementing the independent CaL current.

duration = 100* ms
# Parameters
C = 1 * ufarad
gCaLmax = 0.001 * siemens
Cai = 50 * nmole
Cao = 2 * mmole
temp = 34
KTF = ((25 / 293.15) * (temp + 273.15))
f = KTF * mV / 2
ki = 1 * nmole
h2 = ki / (ki + Cai)
q10 = 5
qt = q10 ** ((temp - 25) / 10)
a0m = 0.1
zetam = 2
vhalfm = 4*mV
gmm = 0.1
V0 = -65 * mV


@check_units(z=1, result=1)
def efun(z):
    if (np.absolute(z) < 0.0001):
        return 1 - z / 2
    else:
        return z / (exp(z) - 1)


@check_units(V=volt, result=1)
def nu(V):
    return V / f


@check_units(V=volt, result=volt)
def ghk(V):
    return -f * (1 - (Cai / Cao) * exp(nu(V))) * efun(nu(V))

#eqs = Equations("""
                   # dV/dt = (ICaL + I)/C : volt
                   # #ICaL = gCaLmax*mCaL*mCaL*h2*ghk(V) : amp
                   # ICaL = gCaL*ghk(V) : amp (constant over dt)
                   # gCaL = gCaLmax*mCaL*mCaL*h2 : siemens
                   # dmCaL/dt = (infmCaL - mCaL)/tau : 1
                   # infmCaL = alp * 1/(alp + bet) : 1
                   # alp = 15.69*(mV**-1)*(-V+81.5*mV)/(exp((-V+81.5*mV)/mV*10)-1) : 1
                   # bet = 0.29 * exp(-V /mV*10.86)
                   # tau = betmt(V)*ms/(qt*a0m*(1+alpmt(V))) : second (constant over dt)
                   # I : amp
                   # """)


eqs = Equations('''
                    dV/dt = (ICaL + I)/C : volt
                    #ICaL = gCaLmax*mCaL*mCaL*h2*ghk(V) : amp (constant over dt)
                    ICaL = gCaL*ghk(V) : amp (constant over dt)
                    gCaL = gCaLmax*mCaL*mCaL*h2 : siemens
                    dmCaL/dt = (infmCaL - mCaL)/tau : 1
                    infmCaL = alp * 1/(alp + bet) : 1
                    alp = 15.69*(mV**-1)*(-V+81.5*mV)/(exp((-V+81.5*mV)/(10*mV))-1) : 1
                    bet = 0.29 * exp(-V/(10.86*mV)) : 1
                    tau = betmt*ms/(qt*a0m*(1+alpmt)) : second 
                    alpmt = exp(0.0378*(mV**-1)*zetam*(V-vhalfm)) : 1
                    betmt = exp(0.0378*(mV**-1)*zetam*gmm*(V-vhalfm)) : 1
                    I : amp                    
                    ''')


group = NeuronGroup(1, eqs, threshold='V > -40*mV', refractory='V > -40*mV', method='rk4')


alp0 = 15.69*(mV**-1)*(-V0+81.5*mV)/(exp((-V0+81.5*mV)/(10*mV))-1)
bet0 = 0.29 * exp(-V0/(10.86*mV))
print(alp0)
print(bet0)
group.mCaL = alp0 * 1/(alp0 + bet0)

M = StateMonitor(group, variables=True, record=True)

#run(1*ms)

store()
# Plot tau and the channel current for different voltages
voltages = np.linspace(-50, 10, 70) * mV
vs = []
Is = []
for volt in voltages:
    # Restore the original state of the network
    restore()
    group.V = volt
    # Run it with the new value of tau
    run(1 * ms)
    vs.append(M.V[0])
    Is.append(M.ICaL[0])

figure(1)
plot(voltages, vs)
xlabel('Voltage')
ylabel('Tau');

figure(2)
plot(voltages, Is)
xlabel('Volt')
ylabel('Amps');
show()