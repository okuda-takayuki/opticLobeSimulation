import sys
import numpy as np
import json
with open("meta.json","r") as f:
    meta = json.load(f)
    sys.path.append(meta["nrnpy-path"])
import neuron

class Medullaneuron:
    def __init__(self, index, opt={"type":"Tm1"},params={}):
        self.index = index
        self.celltype = opt["type"]
        self.params = params
        self.cell={}
        self.generateCell(self.celltype,self.params)
        self.synlist = []

    def generateCell(self, celltype="Tm1", params={}):

        if self.celltype =="Tm1":
            medullafile = "modules/swc/Tm1 home_226.swc"
        if self.celltype =="Tm2":
            medullafile = "modules/swc/Tm2 home_227.swc"
        if self.celltype =="Tm4":
            medullafile = "modules/swc/Tm4-A-post_36.swc"
        if self.celltype =="Tm9":
            medullafile = "modules/swc/Tm9 home_173.swc"
        if self.celltype =="C3":
            medullafile = "modules/swc/C3 home_10.swc"
        if self.celltype =="Mi1":
            medullafile = "modules/swc/Mi1 home_212.swc"
        if self.celltype =="T4":
            medullafile = "modules/swc/T4a-03_167.swc"
        if self.celltype =="T5":
            medullafile = "modules/swc/T5a-06 rep_571.swc"

        mecell = np.loadtxt(medullafile)
        mecell[:,2:6]=0.001*8*mecell[:,2:6]
        interim=mecell.shape
        me_nofcomps=interim[0]
        compdiam=mecell[:,5]*2.0

        for i in range(1, me_nofcomps, 1):
            aind = int(mecell[i,0]-1)
            bind = int(mecell[i,6]-1)
            axyz = mecell[aind,2:5]
            bxyz = mecell[bind,2:5]

            complength=np.sqrt(np.sum((axyz-bxyz)**2))

            meandiam = (compdiam[aind]+compdiam[bind])*0.5

            self.cell[str(aind)] = neuron.h.Section()
            self.cell[str(aind)].diam = meandiam
            self.cell[str(aind)].L = complength
            self.cell[str(aind)].insert ("mole")
            self.cell[str(aind)].cm = 0.6
            self.cell[str(aind)].Ra = 0.1

            if bind == 0:
                continue

            self.cell[str(aind)].connect(self.cell[str(bind)], 1)

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
