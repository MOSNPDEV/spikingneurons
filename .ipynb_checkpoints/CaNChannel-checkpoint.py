from brian2 import *
import brian2.numpy_ as np
prefs.codegen.target = 'numpy'

start_scope()
# Implementing the independent CaN current.

duration = 100 * ms
# Parameters
C = 1 * ufarad
gCaNmax = 0.001 * siemens
Cai = 50 * nmole
Cao = 2 * mmole
q10 = 5
a0m = 0.03 / ms
zetam = 2 / mV
vhalfm = -14 * mV
gmm = 0.1
temp = 34
f = ((25 / 293.15) * (temp + 273.15)) *mV/ 2
V0 = -65 * mV
ki = 1 * nmole
h2 = ki / (ki + Cai)
qt = q10 ** ((temp - 25) / 10)


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


eqs = Equations('''
                    dV/dt = (ICaN + I)/C : volt
                    ICaN = gCaN*ghk(V) : amp (constant over dt)
                    gCaN = gCaNmax*mCaN*mCaN*hCaN*h2 : siemens
                    dmCaN/dt = (infmCaN - mCaN)/taum : 1
                    dhCaN/dt = (infhCaN - hCaN)/tauh : 1
                    infmCaN = alpm*(1/alpm+betm) : 1
                    infhCaN = alph*(1/alph+beth) : 1
                    taum = betmt/(qt*a0m*(1+alpmt)) : second
                    tauh = (1*ms/(alph+beth))/qt : second
                    alpm = 0.1967*(mV**-1)*(-V+19.88*mV)/(exp((-V+19.88*mV)/(10.0*mV))-1) : 1
                    alph = 0.00016*exp(-V/(48.4*mV)) : 1
                    alpmt = exp(0.0378*zetam*(V-vhalfm)) : 1 
                    betm = 0.046*exp(-V/(20.73*mV)) : 1
                    beth = 1/(exp((-V+39*mV)/(10*mV))+1.) : 1
                    betmt = exp(0.0378*zetam*gmm*(V-vhalfm)) : 1
                    I : amp
                    ''')

group = NeuronGroup(1, eqs, threshold='V > -40*mV', refractory='V > -40*mV', method='rk2')

alpm0 = 0.1967 * (mV ** -1) * (-V0 + 19.88 * mV) / (exp((-V0 + 19.88 * mV) / (10 * mV)) - 1)
print(alpm0)
betm0 = 0.046 * exp(65 / 20.73)
print(betm0)
print(alpm0 * (1 / alpm0 + betm0))
group.mCaN = alpm0 * (1 / alpm0 + betm0)
group.mCaN = 0.1
group.hCaN = 0.1
group.hCaN = 0.00016 * exp(-V0 / (48.4 * mV)) * (
1 / 0.00016 * exp(-V0 / (48.4 * mV)) + 1 / (exp((-V0 + 39 * mV) / (10 * mV)) + 1))

M = StateMonitor(group, variables=True, record=True)

run(1 * ms)

figure(1)
plot(M.t / ms, M.mCaN[0])
xlabel('ms')
ylabel('act')

show()