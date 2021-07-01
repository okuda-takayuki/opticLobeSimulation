import sys
import numpy as np
import json
with open("meta.json","r") as f:
    meta = json.load(f)
    sys.path.append(meta["nrnpy-path"])
import neuron


class HHneuron:
    def __init__(self, index,opt={},params={}):
        self.index = index
        self.cell={}
        self.cell["soma"] = neuron.h.Section(name="soma")
        self.axon = neuron.h.Section(name="axon")
        self.axon.nseg = 20
        self.axon.diam = 0.2
        self.axon.L = 10
        self.axon.insert("hh")
        self.ap_dend = neuron.h.Section(name="ap_dend")
        self.ap_dend.L = 10
        self.ap_dend.diam = 2
        self.ap_dend.nseg = 1
        self.ap_dend.insert("hh")
        self.ap_dend.gnabar_hh = 0.012
        self.ap_dend.gkbar_hh = 0.0036
        self.ap_dend.gl_hh = 0.00003
        self.soma.connect(self.axon, 1)
        self.ap_dend.connect(self.soma, 1)
        neuron.h.psection()
        self.esyn = neuron.h.Exp2Syn(self.ap_dend(0.5))
        self.esyn.tau1 = 0.5
        self.esyn.tau2 = 1.0
        self.esyn.e = 0

    def synapticConnection(self, target, setting=[-10, 1, 10]):
        netcon = neuron.h.NetCon(self.axon(0.5)._ref_v, target.esyn, sec=self.axon)
        netcon.threshold = setting[0]
        netcon.weight[0] = setting[1]
        netcon.delay = setting[2]
        return netcon
