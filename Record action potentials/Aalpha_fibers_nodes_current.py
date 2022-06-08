##-----------------------------------------------------------------------------
# General info: generate the action potetial related membrane current from a 
# myelinated Aalpha fiber
##-----------------------------------------------------------------------------
# The scripts is based on the Gaines’s implementation of the classical MRG fiber:
# J. L. Gaines, K. E. Finn, J. P. Slopsema, L. A. Heyboer, and K. H. Polasek, 
# “A model of motor and sensory axon activation in the median nerve using 
# surface electrical stimulation,” J Comput Neurosci, vol. 45, no. 1, pp. 
# 29–43, Aug. 2018, doi: 10.1007/s10827-018-0689-5.
# Repository: https://senselab.med.yale.edu/ModelDB/ShowModel?model=243841#tabs-1
##-----------------------------------------------------------------------------
# Import the neuron module into Python
##-----------------------------------------------------------------------------
import neuron
from neuron import h
from neuron.units import ms, mV
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import scipy.io
##-----------------------------------------------------------------------------
# Create the class
##-----------------------------------------------------------------------------
class MyelinatedFiber(object):
    ##-------------------------------------------------------------------------
    # morphological properties
    length_mysa = 3.0 # um 
    length_ranvier = 1.0 # um 
    periwidth_mysa = 0.002
    periwidth_flut = 0.004
    periwidth_stin = 0.004
    ##-------------------------------------------------------------------------
    # electrical properties
    cm = 2  # [uF/cm2]
    rhoa = 0.7e6  # [Ohm-um]
    cm_my = 0.1  # [uF/cm2-lamella]
    gm_my = 0.001  # [S/cm2-lamella]
    ##-------------------------------------------------------------------------
    # boundary conditions
    h.celsius = 37  # [°C]
    v_init = -85.9411  # [mV]
    ##-------------------------------------------------------------------------
    def __init__(self, diameter, n_internodes):
        # properties depending from self.n_internodes
        self.n_internodes = n_internodes
        self.n_ranvier = n_internodes + 1
        self.n_mysa = n_internodes * 2
        self.n_flut = n_internodes * 2
        self.n_stin = n_internodes * 6
        ##---------------------------------------------------------------------
        self.diameter = diameter
        ##---------------------------------------------------------------------
        # properties depending from self.diameter
        self.g = 0.0172 * diameter + 0.5076
        self.diam_stin = 0.889 * diameter - 1.9104  # diameter of the axon
        self.diam_ranvier = 0.3449 * diameter - 0.1484  # diameter of the node
        self.diam_mysa = 0.3527 * diameter - 0.1804  # diameter of mysa
        self.diam_flut = 0.889 * diameter - 1.9104  # diameter of flut
        self.length_internode = 969.3 * np.log(diameter) - 1144.6  # total length between nodes (including 1/2 the node on each side)
        self.length_flut = 2.5811 * diameter + 19.59  # length of flut
        self.n_lamellae = 65.897 * np.log(diameter) - 32.666  # number of lamellae
        ##---------------------------------------------------------------------
        self.length_stin = (self.length_internode - self.length_ranvier - (2 * self.length_mysa) - (2 * self.length_flut)) / 6
        ##---------------------------------------------------------------------
        # periaxonal resistivity computation
        ##---------------------------------------------------------------------
        area_ext_ranvier = np.pi * (self.diam_ranvier / 2 + self.periwidth_mysa)**2
        area_int_ranvier = np.pi * (self.diam_ranvier / 2)**2
        periarea_ranvier = area_ext_ranvier - area_int_ranvier # um^2
        ##---------------------------------------------------------------------
        area_ext_mysa = np.pi * (self.diam_mysa / 2 + self.periwidth_mysa)**2
        area_int_mysa = np.pi * (self.diam_mysa / 2)**2
        periarea_mysa = area_ext_mysa - area_int_mysa # um^2
        ##---------------------------------------------------------------------
        area_ext_flut = np.pi * (self.diam_flut / 2 + self.periwidth_flut)**2
        area_int_flut = np.pi * (self.diam_flut / 2)**2
        periarea_flut = area_ext_flut - area_int_flut # um^2
        ##---------------------------------------------------------------------
        area_ext_stin = np.pi * (self.diam_stin / 2 + self.periwidth_stin)**2
        area_int_stin = np.pi * (self.diam_stin / 2)**2
        periarea_stin = area_ext_stin - area_int_stin # um^2
        ##---------------------------------------------------------------------
        self.R_pa_ranvier = (self.rhoa * 0.01) / periarea_ranvier # Ohm
        self.R_pa_mysa = (self.rhoa * 0.01) / periarea_mysa # Ohm
        self.R_pa_flut = (self.rhoa * 0.01) / periarea_flut # Ohm
        self.R_pa_stin = (self.rhoa * 0.01) / periarea_stin # Ohm
    ##-------------------------------------------------------------------------
    def generate_fiber(self):
        # generate sections
        self.Ranvier = [h.Section() for i in range(self.n_ranvier)]
        self.MYSA = [h.Section() for i in range(self.n_mysa)]
        self.FLUT = [h.Section() for i in range(self.n_flut)]
        self.STIN = [h.Section() for i in range(self.n_stin)]
        ##---------------------------------------------------------------------
        self.nodes = []
        self.x = []
        self.x.append(0)
        for i in range(self.n_internodes):
            if i > 0:
                self.x.append((self.length_internode+self.x[i]))
            else:
                self.x.append(self.length_internode)
            self.nodes.append(self.Ranvier[i])
            self.nodes.append(self.MYSA[2*i])
            self.nodes.append(self.FLUT[2*i])
            self.nodes.append(self.STIN[6*i])
            self.nodes.append(self.STIN[6*i + 1])
            self.nodes.append(self.STIN[6*i + 2])
            self.nodes.append(self.STIN[6*i + 3])
            self.nodes.append(self.STIN[6*i + 4])
            self.nodes.append(self.STIN[6*i + 5])
            self.nodes.append(self.FLUT[2*i + 1])
            self.nodes.append(self.MYSA[2*i + 1])
        self.nodes.append(self.Ranvier[self.n_internodes])
        self.x = 1e-6*(self.x - (self.x[n_internodes-1]/2))
        ##---------------------------------------------------------------------
        # set biophysics
        for sec in h.allsec():
            sec.insert('extracellular')
            sec.e_extracellular = 0
        ##---------------------------------------------------------------------
        for ranvier in self.Ranvier:
            ranvier.nseg = 1
            ranvier.diam = self.diam_ranvier
            ranvier.L = self.length_ranvier
            ranvier.Ra = self.rhoa / 1e4
            ranvier.cm = self.cm
            ranvier.insert('node_motor')
            ranvier.xraxial[0] = self.R_pa_ranvier # longitudinal resistances MOhms/cm 
            ranvier.xg[0] = 1e10 # net radial conductance mho/cm2
            ranvier.xc[0] = 0 # net radial capacitance uF/cm2 
        ##---------------------------------------------------------------------
        for mysa in self.MYSA:
            eta_diam = self.diam_mysa / self.diameter
            mysa.nseg = 1
            mysa.diam = self.diameter
            mysa.L = self.length_mysa
            mysa.Ra = (self.rhoa / eta_diam**2) / 1e4
            mysa.cm = self.cm * eta_diam
            mysa.insert('mysa_motor')
            mysa.xraxial[0] = self.R_pa_mysa
            mysa.xg[0] = self.gm_my / (2 * self.n_lamellae)
            mysa.xc[0] = self.cm_my / (2 * self.n_lamellae)
        ##---------------------------------------------------------------------    
        for flut in self.FLUT:
            eta_diam = self.diam_flut / self.diameter
            flut.nseg = 1
            flut.diam = self.diameter
            flut.L = self.length_flut
            flut.Ra = (self.rhoa / eta_diam**2) / 1e4 
            flut.cm = self.cm * eta_diam
            flut.insert('flut_motor')
            flut.xraxial[0] = self.R_pa_flut
            flut.xg[0] = self.gm_my / (2 * self.n_lamellae)
            flut.xc[0] = self.cm_my / (2 * self.n_lamellae)
        ##---------------------------------------------------------------------    
        for stin in self.STIN:
            eta_diam = self.diam_stin / self.diameter
            stin.nseg = 1
            stin.diam = self.diameter
            stin.L = self.length_stin
            stin.Ra = (self.rhoa / eta_diam**2) / 1e4
            stin.cm = self.cm * eta_diam
            stin.insert('stin_motor')
            stin.xraxial[0] = self.R_pa_stin
            stin.xg[0] = self.gm_my / (2 * self.n_lamellae)
            stin.xc[0] = self.cm_my / (2 * self.n_lamellae)
        ##---------------------------------------------------------------------    
        # set connectivity
        for i in range(self.n_internodes):
            self.MYSA[2*i].connect(self.Ranvier[i](1), 0)
            self.FLUT[2*i].connect(self.MYSA[2*i](1), 0)
            self.STIN[6*i].connect(self.FLUT[2*i](1), 0)
            self.STIN[6*i + 1].connect(self.STIN[6*i](1), 0)
            self.STIN[6*i + 2].connect(self.STIN[6*i + 1](1), 0)
            self.STIN[6*i + 3].connect(self.STIN[6*i + 2](1), 0)
            self.STIN[6*i + 4].connect(self.STIN[6*i + 3](1), 0)
            self.STIN[6*i + 5].connect(self.STIN[6*i + 4](1), 0)
            self.FLUT[2*i + 1].connect(self.STIN[6*i + 5](1), 0)
            self.MYSA[2*i + 1].connect(self.FLUT[2*i + 1](1), 0)
            self.Ranvier[i+1].connect(self.MYSA[2*i + 1](1), 0)
    ##-------------------------------------------------------------------------    
    def set_intra_stim(self, stim_node, amp, delay, dur, tstop):
        stim = h.IClamp(self.Ranvier[stim_node](0.5))
        stim.amp = amp
        stim.delay = delay
        stim.dur = dur
        ic_rec = []        
        for i in range(self.n_ranvier):
            ic_rec.append(h.Vector().record(self.Ranvier[i](0.5)._ref_i_cap)) # mA/cm2
        t_rec = h.Vector().record(h._ref_t)
        ##---------------------------------------------------------------------
        # model parameters
        h.load_file('stdrun.hoc')
        h.finitialize(self.v_init * mV)
        h.dt = 0.005
        h.continuerun(tstop * ms)
        return ic_rec, t_rec
##-----------------------------------------------------------------------------
# Main
##-----------------------------------------------------------------------------    
for d in range(9,21):
    diameter = d
    length_internode = 969.3 * np.log(diameter) - 1144.6  # total length between nodes (including 1/2 the node on each side)
    total_length = 10000
    n_internodes = int(round(total_length/length_internode)+2*30)
    dx = length_internode
    delay = 50
    dur = 3 # ms
    amp = 1 # mA
    fiber = MyelinatedFiber(diameter, n_internodes)
    fiber.generate_fiber()
    ic, t = fiber.set_intra_stim(0, amp, delay, dur, 80)
    area = fiber.Ranvier[1](0.5).area() # square microns (2*pi*r*Δx, where r is the axonic radius and Δx is the segment length)
    ##-------------------------------------------------------------------------
    for i in range(fiber.n_ranvier):
        if i == 0:
            arr1 = ic[i]
            arr2 = ic[i+1]
            data = np.concatenate((arr1,arr2), axis=0)
        if i == 1:
            arr2 = ic[i]
        else:
            arr2 = ic[i]
            data = np.concatenate((data,arr2), axis=0)
    data = data.reshape(int(len(data)/len(arr2)), len(arr2))
    ##-------------------------------------------------------------------------
    scipy.io.savemat('Aalpha_fiber_' + str(diameter) + '.mat', {'data': data, 'area': area, 'num_node': fiber.n_ranvier, 'dt': h.dt, 'diameter': diameter, 'x': fiber.x})
