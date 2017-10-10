import matplotlib.pyplot as plt
import pandas as pd
resAll = pd.read_hdf("./AO4ELT3.h5", "resAll")

"""
To save vectoriel files:
plt.savefig("plot.pdf", transparent=True, bbox_inches="tight")


"""

"""
SH VS PYR
"""

PYR01_GSMAG = [10.997543 ,12.009828 ,13.503686 ,14.02457 ,14.506143 ,15.027027 ,15.498772 ,16.019655 ,16.501228 ,17.002457, 17.5, 17.995087, 18.5, 19.0]
PYR01 = [60.5625,60.1875,59.71875 ,58.875,58.03125,55.96875,54.09375,49.59375,45.65625,40.03125,31.1 ,19.5, 8, 0.84375]

PYR3_GSMAG = [11.007371,12.0,13.012285,13.58231,14.014742,14.51597,15.007371,15.518428,16.009829,16.501228,17.002457,17.513514,18.02457,18.51597]
PYR03 = [60.1875 ,59.90625,59.625,58.03125,56.4375,54.28125,49.3125,40.03125,18.65625,4.40625,0.1875,0.1875,0.09375,0.09375]

SH3_GSMAG = [11.007371,12.0,13.012285,13.503686,14.014742,14.51597,15.007371,15.498772,16.019655,16.992628,17.513514,18.004913]
SH3 = [57.4687,56.8125,54.09375,51.1875,45.65625,35.34375,21.84375,7.6875,1.40625,0.09375,0.0,0]


plt.figure(18)
plt.clf()
magoffset = 1
plt.plot(np.array(PYR01_GSMAG)-magoffset, PYR01, color="green", marker= "d", ls = "--",label="PYRAMID 0.1e- RON")
plt.plot(np.array(PYR3_GSMAG)-magoffset, PYR03, color="green", marker= "s", ls = "-", label="PYRAMID 3e- RON")
plt.plot(np.array(SH3_GSMAG)-magoffset, SH3, color="blue", marker= "v", ls = "-", label="Shack-Hartmann 3e- RON")
plt.ylabel("Strehl Ratio (H band)")
plt.xlabel("GS Magnitude")
plt.xlim(10, 19)
plt.ylim(0, 65)
plt.grid()
plt.legend()





"""
CN2 profiles

"""

zenithAngle = 30.
altESO = np.array([30,90,150,200,245,300,390,600,1130,1880,2630,3500,4500,5500,6500,7500,8500,9500,10500,11500,12500,13500,14500,15500,16500,17500,18500,19500,20500,21500,22500,23500,24500,25500,26500])/np.cos(zenithAngle*2*np.pi/360)
altESO = altESO.astype(int)

fracmed = [24.2,12,9.68,5.9,4.73,4.73,4.73,4.73,3.99,3.24,1.62,2.6,1.56,1.04,1,1.2,0.4,1.4,1.3,0.7,1.6,2.59,1.9,0.99,0.62,0.4,0.25,0.22,0.19,0.14,0.11,0.06,0.09,0.05,0.04]
fracQ1 = [22.6,11.2,10.1,6.4,4.15,4.15,4.15,4.15,3.1,2.26,1.13,2.21,1.33,0.88,1.47,1.77,0.59,2.06,1.92,1.03,2.3,3.75,2.76,1.43,0.89,0.58,0.36,0.31,0.27,0.2,0.16,0.09,0.12,0.07,0.06]
fracQ2 = [25.1,11.6,9.57,5.84,3.7,3.7,3.7,3.7,3.25,3.47,1.74,3,1.8,1.2,1.3,1.56,0.52,1.82,1.7,0.91,1.87,3.03,2.23,1.15,0.72,0.47,0.3,0.25,0.22,0.16,0.13,0.07,0.11,0.06,0.05]
fracQ3 = [25.5,11.9,9.32,5.57,4.5,4.5,4.5,4.5,4.19,4.04,2.02,3.04,1.82,1.21,0.86,1.03,0.34,1.2,1.11,0.6,1.43,2.31,1.7,0.88,0.55,0.36,0.22,0.19,0.17,0.12,0.1,0.06,0.08,0.04,0.04]
fracQ4 = [23.6,13.1,9.81,5.77,6.58,6.58,6.58,6.58,5.4,3.2,1.6,2.18,1.31,0.87,0.37,0.45,0.15,0.52,0.49,0.26,0.8,1.29,0.95,0.49,0.31,0.2,0.12,0.1,0.09,0.07,0.06,0.03,0.05,0.02,0.02]

plt.figure(10)
plt.clf()
plt.plot(fracQ1, altESO, color="green", marker= "H", label="Q1")
plt.plot(fracQ2, altESO, color="blue", marker= "o", label="Q2")
plt.plot(fracQ3, altESO, color="orange", marker= "d", label="Q3")
plt.plot(fracQ4, altESO, color="red", marker= "s", label="Q4")
plt.plot(fracmed, altESO, color="black", marker= "p", label="Median")
plt.ylabel("Altitude (m)")
plt.xlabel("layer fraction (%)")
plt.grid()
plt.legend()

"""
SR VS GSMAG (J, H, K)
"""


sizetext = 15
resMedian = resAll[resAll.comment == "SRVsGSmag_Median"]
# Finding best gain for perf for each gsmag:
ind = []
#for i in range(len(resMedian)):
for i in list(set(resMedian.gsmag)):
    ind.append((resMedian["SR_2.20"][resMedian.gsmag == int(i)]).idxmax())
resMedianBest = resMedian.loc[ind]
plt.ion()
plt.figure(0)
plt.clf()
plt.plot(resMedianBest["gsmag"], resMedianBest["SR_2.20"], color= "red", marker= "s", label = "K Band")
plt.plot(resMedianBest["gsmag"], resMedianBest["SR_1.65"], color= "blue", marker= "o", label = "H Band")
plt.plot(resMedianBest["gsmag"], resMedianBest["SR_1.20"], color= "green", marker= "v", label = "J Band")
plt.legend()
plt.grid()
plt.ylabel("Strehl ratio")
plt.xlabel("Guide star Magnitude")
plt.ylim(0, 0.9)
plt.text(12, 0.1, r'Median seeing conditions (0.702")', size=sizetext)


resQ1 = resAll[resAll.comment == "SRVsGSmag_Q1"]
# Finding best gain for perf for each gsmag:
ind = []
#for i in range(len(resMedian)):
for i in list(set(resQ1.gsmag)):
    ind.append((resQ1["SR_2.20"][resQ1.gsmag == int(i)]).idxmax())
resQ1Best = resQ1.loc[ind]

plt.ion()
plt.figure(1)
plt.clf()
plt.plot(resQ1Best["gsmag"], resQ1Best["SR_2.20"], color= "red", marker= "s", label = "K Band")
plt.plot(resQ1Best["gsmag"], resQ1Best["SR_1.65"], color= "blue", marker= "o", label = "H Band")
plt.plot(resQ1Best["gsmag"], resQ1Best["SR_1.20"], color= "green", marker= "v", label = "J Band")
plt.legend()
plt.grid()
plt.ylabel("Strehl ratio")
plt.xlabel("Guide star Magnitude")
plt.ylim(0, 0.9)
plt.text(12, 0.1, r'Q1 turbulence conditions (0,471")', size=sizetext)



resQ2 = resAll[resAll.comment == "SRVsGSmag_Q2"]
# Finding best gain for perf for each gsmag:
ind = []
#for i in range(len(resMedian)):
for i in list(set(resQ2.gsmag)):
    ind.append((resQ2["SR_2.20"][resQ2.gsmag == int(i)]).idxmax())
resQ2Best = resQ2.loc[ind]
plt.ion()
plt.figure(2)
plt.clf()
plt.plot(resQ2Best["gsmag"], resQ2Best["SR_2.20"], color= "red", marker= "s", label = "K Band")
plt.plot(resQ2Best["gsmag"], resQ2Best["SR_1.65"], color= "blue", marker= "o", label = "H Band")
plt.plot(resQ2Best["gsmag"], resQ2Best["SR_1.20"], color= "green", marker= "v", label = "J Band")
plt.legend()
plt.grid()
plt.ylabel("Strehl ratio")
plt.xlabel("Guide star Magnitude")
plt.ylim(0, 0.9)
plt.text(12, 0.1, r'Q2 turbulence conditions (seeing = 0,619")', size=sizetext)



resQ3 = resAll[resAll.comment == "SRVsGSmag_Q3"]
# Finding best gain for perf for each gsmag:
ind = []
#for i in range(len(resMedian)):
for i in list(set(resQ3.gsmag)):
    ind.append((resQ3["SR_2.20"][resQ3.gsmag == int(i)]).idxmax())
resQ3Best = resQ3.loc[ind]
plt.ion()
plt.figure(3)
plt.clf()
plt.plot(resQ3Best["gsmag"], resQ3Best["SR_2.20"], color= "red", marker= "s", label = "K Band")
plt.plot(resQ3Best["gsmag"], resQ3Best["SR_1.65"], color= "blue", marker= "o", label = "H Band")
plt.plot(resQ3Best["gsmag"], resQ3Best["SR_1.20"], color= "green", marker= "v", label = "J Band")
plt.legend()
plt.grid()
plt.ylabel("Strehl ratio")
plt.xlabel("Guide star Magnitude")
plt.ylim(0, 0.9)
plt.text(12, 0.1, r'Q3 turbulence conditions (seeing = 0,793")', size=sizetext)



resQ4 = resAll[resAll.comment == "SRVsGSmag_Q4"]
# Finding best gain for perf for each gsmag:
ind = []
#for i in range(len(resMedian)):
for i in list(set(resQ4.gsmag)):
    ind.append((resQ4["SR_2.20"][resQ4.gsmag == int(i)]).idxmax())
resQ4Best = resQ4.loc[ind]
plt.ion()
plt.figure(4)
plt.clf()
plt.plot(resQ4Best["gsmag"], resQ4Best["SR_2.20"], color= "red", marker= "s", label = "K Band")
plt.plot(resQ4Best["gsmag"], resQ4Best["SR_1.65"], color= "blue", marker= "o", label = "H Band")
plt.plot(resQ4Best["gsmag"], resQ4Best["SR_1.20"], color= "green", marker= "v", label = "J Band")
plt.legend()
plt.grid()
plt.ylabel("Strehl ratio")
plt.xlabel("Guide star Magnitude")
plt.ylim(0, 0.9)
plt.text(11.5, 0.65, r'Q4 turbulence conditions (seeing = 1,136")', size=sizetext)





"""
Modulation size

"""
resModu = resAll[resAll.comment == "SRVsGsMagVsModRadius_Q3"]
plt.figure(11)
plt.clf()
mods = resModu["pyr_ampl"].values
vals, ind = np.unique(mods, return_index=True)
ind = ind.tolist()
ind.append(len(resModu["gsmag"]))
mk = ['o', 'v', '^', '<', '>', '8', 's', 'p']
ls =["--", "-", "-.", ":", "-", "--", "-.", ":"]
for i in range(len(ind)-1):
    st=ind[i]
    end=ind[i+1]
    plt.plot(resModu["gsmag"][st:end], resModu["SR_2.20"][st:end], marker= mk[i], ls=ls[i], label = 'Pyr Mod='+str(int(vals[i]))+r'$\lambda/D$')
    #plt.plot(resModu["gsmag"][st:end], resModu["rmsError"][st:end], marker= "x", label = "Modulation="+str(vals[i][0]))
plt.legend(loc='lower left')
plt.ylabel("K Band Strehl ratio")
plt.xlabel("Guide star Magnitude")
plt.grid()


"""

SR Vs Nsubap

"""
plt.figure(17)
plt.clf()
resNssp = resAll[resAll.comment == "SRVsNSSP"][2:]
resNssp = resNssp.reset_index()
vals, ind = np.unique(resNssp["gsmag"], return_index=True)
ind = ind.tolist()
ind.append(len(resModu["gsmag"]))
markerList = ['o', 's', 'v']
for i in range(len(ind)-1):
    st=ind[i]
    end=ind[i+1]
    plt.plot(resNssp["nxsub"][st:end], resNssp["SR_2.20"][st:end], marker= markerList[i], label = 'GS Mag='+str(int(vals[i])))
    #plt.plot(resModu["gsmag"][st:end], resModu["rmsError"][st:end], marker= "x", label = "Modulation="+str(vals[i][0]))
plt.legend(loc='lower right')
plt.ylabel("K Band Strehl ratio")
plt.xlabel("Pyramid Nsubp (npix per Pupil)")
plt.grid()
plt.ylim(0, 0.9)



"""

SR VS off axis FOV

"""
resAll = pd.read_hdf("./AO4ELT5_fov.h5", "resAll")
plt.figure(18)
plt.ion()
plt.clf()
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][0][0:13], color="red", label= "K Band @ Magnitude 11", ls="-", marker="o")
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][0][13:13*2], color="blue", label= "H Band @ Magnitude 11",ls="-", marker="s")
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][0][13*2:13*3], color="green", label= "J Band @ Magnitude 11", ls="-", marker="v")


plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][1][0:13], color="red", label= "K Band @ Magnitude 14", ls="--", marker="o")
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][1][13:13*2], color="blue", label= "H Band @ Magnitude 14",ls="--", marker="s")
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][1][13*2:13*3], color="green", label= "J Band @ Magnitude 14", ls="--", marker="v")

plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][2][0:13], color="red", label= "K Band @ Magnitude 17",ls="-.", marker="o")
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][2][13:13*2], color="blue", label= "H Band @ Magnitude 17",ls="-.", marker="s")
plt.plot([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], resAll["srir"][2][13*2:13*3], color="green", label= "J Band @ Magnitude 17",ls="-.", marker="v")
plt.legend()
plt.xlabel("Off-axis distance from Guide Star ('')")
plt.ylabel("Strehl ratio")
plt.grid()




"""

SR VS # Controlled modes @ mag 11, 15, 17

"""
plt.figure(19)
plt.ion()
plt.clf()
resAll = pd.read_hdf("./AO4ELT3_NModeFilt.h5", "resAll")
resNModes = resAll[resAll.comment == "SRVsGSVsNControlledModes"]
resNModes = resNModes.reset_index()

vals, ind = np.unique(resNModes["gsmag"], return_index=True)
ind = ind.tolist()
ind.append(len(resNModes["gsmag"]))
markerList = ['o', 's', 'v']
for i in range(len(ind)-1):
    st=ind[i]
    end=ind[i+1]
    plt.plot(resNModes["NklFilt"][st:end], resNModes["SR_2.20"][st:end], marker= markerList[i], label = 'GS Mag='+str(int(vals[i])))
    #plt.plot(resModu["gsmag"][st:end], resModu["rmsError"][st:end], marker= "x", label = "Modulation="+str(vals[i][0]))
plt.legend(loc='lower right')
plt.ylabel("K Band Strehl ratio")
plt.xlabel("# of filtered modes (total=4473 modes)")
plt.grid()
plt.xticks(resNModes["NklFilt"][st:end].tolist())
