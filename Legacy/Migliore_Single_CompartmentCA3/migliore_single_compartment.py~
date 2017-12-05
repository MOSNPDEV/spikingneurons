############################################################################
# Test program for replicating Migliore NEURON firing patterns simulation in a single-compartment neuron (Ca3 pyramidal cell)
# Adrian Gutierrez agutie@ini.uzh.ch
# Modified by Saray Soldado ssaray@ini.uzh.ch 
############################################################################  
from neuron import h, gui
import numpy as np
from matplotlib import pyplot
import random 
import scipy.io as sio

##Create soma
soma = h.Section(name='soma')

##Modify soma's properties
soma.L = 50.0      #Length of soma (um)
soma.diam = 50.0   #Diameter of soma (um)
h.celsius = 34
h.steps_per_ms= 40 #Default neuron integration
h.dt=0.025

##Insert passive biophysical mechanisms to soma
soma.insert('pas')
soma.Ra = 150.0            #Membrane axial resistance (Ohm/cm^2) 
soma.cm = 1.41             #Membrane capacitance 
soma.g_pas = 1.0/25370.0   #Leak maximal conductance 


##Insert active biophysical mechanisms to soma (inserting .mod files)

for channel in ['na3', 'kdr', 'kap', 'km', 'kd', 
		'cacum', 'cal', 'can', 'cat', 'cagk']:
	soma.insert(channel)


###############################################
# Setup model conductances
################################################3

#Which pattern are we testing
pattern = 'iBurst'

soma.gbar_na3= 0.04
soma.gkdrbar_kdr=0.01
soma.gkabar_kap=0.07
soma.gcatbar_cat=0.001
soma.gcanbar_can=0.001
soma.gcalbar_cal=0.001
soma.gbar_cagk=0.0001
soma.gbar_km=0.0006
soma.gkdbar_kd=0.00045

#Reversal potentials and calcium
soma.depth_cacum = 50.0/2.0     #Set the point where calcium mechanisms act (Soma diameter/2)
soma.ek = -90.0             #Potassium reversal potential (mV)
soma.ena = 55.0             #Sodium reversal potential (mV)
#soma.cai0_cacum=50e-6 #this is from migliore, but we did not include it, so no calcium inside the cell at t=0...

###############################################
# Setup recording vectors
################################################3

v_vec = h.Vector()   #Membrane potential vector
t_vec = h.Vector()   #Time stamp vector
i_vec = h.Vector()   #Leak current vector 

i_vec_ct = h.Vector()   #CaT channel current vector
i_vec_kd = h.Vector()   #Kd channel current vector
i_vec_km = h.Vector()   #Km channel current vector

v_vec.record(soma(0.5)._ref_v)  #Recording from the soma the desired quantities
t_vec.record(h._ref_t)
i_vec.record(soma(0.5)._ref_i_pas)

i_vec_ct.record(soma(0.5)._ref_ica_cat)
i_vec_kd.record(soma(0.5)._ref_ik_kd)
i_vec_km.record(soma(0.5)._ref_ik_km)

###############################################
# Current Clamp
################################################

st = h.IClamp(0.5)   #Choose where in the soma to point-stimulate 
	
st.dur = 1000    #Stimulus duration (ms) 
st.delay = 200    #Stimulus delay (ms)
st.amp = 1    #Stimulus amplitude (nA)
h.tstop = 1400  #stop the simulation (ms)


##########################################################################################################
# Run the simulation
##########################################################################################################

h.v_init = -65         #Set initializing simulation voltage (mV) at t0
h.finitialize(-65)     #Set initializing voltage for all mechanisms in the section
h.fcurrent()
if h.ismembrane('cal'):  #start the 'vClamp' by setting e_pas to match v_init (NEURON Book, ch8, p11) 
	soma.e_pas = soma.v + (soma.ina + soma.ik + soma.ica)/soma.g_pas
else:
	soma.e_pas = soma.v + (soma.ina + soma.ik)/soma.g_pas

h.run()	


###################################################
# PLOTTING
###################################################

#Initialize plot
fig1 = pyplot.figure(figsize=(12,7))
fig1.suptitle(pattern + ' firing pattern (IClamp)', fontweight='bold')
ax1 = fig1.add_subplot(321)
ax2 = fig1.add_subplot(322)
ax3 = fig1.add_subplot(323)
ax4 = fig1.add_subplot(324)
ax5 = fig1.add_subplot(325)


#Plot time course of signals of interest
ax1.plot(t_vec, v_vec)
ax2.plot(t_vec, i_vec)
ax3.plot(t_vec, i_vec_ct)
ax4.plot(t_vec, i_vec_km)
ax5.plot(t_vec, i_vec_kd)

#Plotting settings for Fig1		
ax1.set_title('Membrane voltage (mV)')
ax2.set_title('Passive (membrane) current (mA/cm2)')
ax3.set_title('Current through CaT channel (mA/cm2)')
ax4.set_title('Current through Km channel (mA/cm2)')
ax4.set_xlabel('Time (ms)')
ax5.set_title('Current through Kd channel (mA/cm2)')
ax5.set_xlabel('Time (ms)')

ax1.axes.get_xaxis().set_visible(False)
ax2.axes.get_xaxis().set_visible(False)
ax3.axes.get_xaxis().set_visible(False)

pyplot.show()

