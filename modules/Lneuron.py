import sys
import numpy as np
sys.path.append("/home/hayato/lib/python")
import neuron


class Lneuron:
    def __init__(self,index,opt={"type":"L1"},params={}):
        self.index = index
        self.celltype = opt["type"]
        self.params = params
        self.cell = {}
        self.generateCell(self.celltype, self.params)
        self.synlist = []

    def generateCell(self, celltype="L1", params={}):

        if self.celltype =="L1":
            laminafile = "modules/swc/L1 home_3529395.swc"
        if self.celltype =="L2":
            laminafile = "modules/swc/L2 home_10312912.swc"
        if self.celltype =="L4":
            laminafile = "modules/swc/L4 home_10317362.swc"
        if self.celltype =="L5":
            laminafile = "modules/swc/L5 home_3537643.swc"

        lacell = np.loadtxt(laminafile)
        lacell[:,2:6]=0.001*8*lacell[:,2:6]
        interim=lacell.shape
        la_nofcomps=interim[0]
        compdiam=lacell[:,5]*2.0

        for i in range(1, la_nofcomps, 1):
            aind = int(lacell[i,0]-1)
            bind = int(lacell[i,6]-1)
            axyz = lacell[aind,2:5]
            bxyz = lacell[bind,2:5]


            complength=np.sqrt(np.sum((axyz-bxyz)**2))

            meandiam = (compdiam[aind]+compdiam[bind])*0.5

            self.cell[str(aind)] = neuron.h.Section()
            self.cell[str(aind)].diam = meandiam
            self.cell[str(aind)].L = complength
            self.cell[str(aind)].insert ("LMC")
            self.cell[str(aind)].cm = 0.6
            self.cell[str(aind)].Ra = 0.1
            self.cell[str(aind)].gl_LMC = 97

            if bind == 0:
                continue

            self.cell[str(aind)].connect(self.cell[str(bind)], 1)
