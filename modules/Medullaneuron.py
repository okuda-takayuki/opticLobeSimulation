import sys
import numpy as np
import json
with open("meta.json","r") as f:
    meta = json.load(f)
    sys.path.append(meta["nrnpy-path"])
import neuron

center = 0.6
tm9 = 1600000
ct1 = 1
dif = 0


class Medullaneuron:
    def __init__(self, index, opt={"type":"Tm1"},params={}):
        self.index = index
        self.celltype = opt["type"]
        self.params = params
        self.cell={}
        self.generateCell(self.celltype,self.params)
        self.synlist = []

    def generateCell(self, celltype="Tm1", params={}):
        if self.celltype == "Tm1":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.4 - dif
            self.cell["axon"].L = 140
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = center
            self.cell["axon"].Ra = 0.01

            #self.cell["ap_dend"] = neuron.h.Section(name="ap_dend")
            #self.cell["ap_dend"].nseg = 9
            #self.cell["ap_dend"].L = 100
            #self.cell["ap_dend"].diam = 0.2
            #self.cell["ap_dend"].insert("mole")
            #self.cell["ap_dend"].cm = center
            #self.cell["ap_dend"].Ra = 0.01

            #self.cell["ap_dend"].connect(self.cell["axon"], 0.6)


        if self.celltype == "Tm2":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.4 - dif
            self.cell["axon"].L = 140
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = center
            self.cell["axon"].Ra = 0.01

            #self.cell["ap_dend"] = neuron.h.Section(name="ap_dend")
            #self.cell["ap_dend"].nseg = 9
            #self.cell["ap_dend"].L = 100
            #self.cell["ap_dend"].diam = 0.2
            #self.cell["ap_dend"].insert("mole")
            #self.cell["ap_dend"].cm = center
            #self.cell["ap_dend"].Ra = 0.01

            #self.cell["ap_dend"].connect(self.cell["axon"], 0.6)

        if self.celltype == "Tm4":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.4 - dif
            self.cell["axon"].L = 140
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = center
            self.cell["axon"].Ra = 0.01

            #self.cell["ap_dend"] = neuron.h.Section(name="ap_dend")
            #self.cell["ap_dend"].nseg = 9
            #self.cell["ap_dend"].L = 100
            #self.cell["ap_dend"].diam = 0.2
            #self.cell["ap_dend"].insert("mole")
            #self.cell["ap_dend"].cm = center
            #self.cell["ap_dend"].Ra = 0.01

            #self.cell["ap_dend"].connect(self.cell["axon"], 0.6)

        if self.celltype == "Tm9":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.4
            self.cell["axon"].L = 140
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = center
            self.cell["axon"].Ra = 0.01


        if self.celltype == "Mi1":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.52 - dif
            self.cell["axon"].L = 134
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = center
            self.cell["axon"].Ra = 0.01
        if self.celltype == "CT1":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.3 - dif
            self.cell["axon"].L = 78
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = ct1
            self.cell["axon"].Ra = 0.01
        if self.celltype == "C3":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.378 - dif
            self.cell["axon"].L = 157
            self.cell["axon"].insert("mole")
            self.cell["axon"].cm = center
            self.cell["axon"].Ra = 0.01
                
        if self.celltype == "T4":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.38
            self.cell["axon"].L = 106
            self.cell["axon"].insert("T4")
            self.cell["axon"].gl_T4 = 8000
            self.cell["axon"].Ra = 0.01
        if self.celltype == "T5":
            self.cell["axon"] = neuron.h.Section(name="axon")
            self.cell["axon"].nseg = 9
            self.cell["axon"].diam = 0.38
            self.cell["axon"].L = 97
            self.cell["axon"].insert("T4")
            self.cell["axon"].gl_T4 = 8000
            self.cell["axon"].Ra = 0.01


    def synapticConnection(self, connection_gid=0,type="E",pc=None,position=["soma",0.5]):
        syn = self.generateSynapse(type=type,position=position)
        pc.target_var(syn,syn._ref_vpre, connection_gid)

    def generateSynapse(self,type="E",position=["soma",0.5]):
        if position[1] > 1 or position[1] < 0:
            print("Synaptic position out of range." + str(position[1]))
            exit()
        syn = neuron.h.gsyn(self.cell[position[0]](position[1]))
        if type == "E":
            # gsat/k equals dynamic range
            syn.vth = -58
            syn.gsat = 0.08
            syn.k = 0.02
            syn.n = 1
            syn.numsyn = 10
            syn.vre = -40
        elif type == "I":
            syn.vth = -58
            syn.gsat = 0.08
            syn.k = 0.02
            syn.n = 1
            syn.numsyn = 10
            syn.vre = -70

        self.synlist.append(syn)
        return syn
