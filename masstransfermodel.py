import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import log, exp

#given constants
timestep = 0.1 #time step in sec
t_min = 3 #time in minutes
t_sec = int(t_min*60) #time in seconds
kFe = 0.0025 #mass tranfer coeffiC_i_net of Fe
kCr = 0.0025 #mass transfer coeffiC_i_ent of Cr
Area = 4 #area of reaC_i_on
t = [] #t stack value
for i in range(164+1): #(t_sec/timestep + 1):
    t.append(i*timestep)
mmass_initial = 150000
Fe_initial_percent = 84.9
Cr_initial_percent = 15
O_initial_percent = 100 - Fe_initial_percent - Cr_initial_percent
total_Oflow = 300
smass_initial = 100
FeO_initial = 50
Cr2O3_initial = 100-FeO_initial
mDen = 7000 #density of metal = 7000kg/m^3
sDen = 3500 #density of slag = 3500kg/m^3
temp_initial = 1600 #initial temperature

# list for changing quantities in each time step
mass_slag = [] #mass of slag
mass_metal = [] #mass of metal
mass_Fe = [] #mass of Fe
mass_Cr = [] #mass of Cr
mass_O = [] #mass of O
mass_FeO = [] #mass of FeO
mass_Cr2O3 = [] #mass of Cr2O3
kmolM =[] #kilomoles of metal
kmolS = [] #kilomoles of slag
kmolFe = [] #kilomoles of Fe
kmolCr = [] #kilomoles of Cr
kmolO = [] #kilomoles of O
kmolCr2O3 = [] #kilomoles of Cr2O3
volM = [] #volume of metal
volS = [] #volume of slag
kmolFeO = [] #kilomoles of FeO
X_i_Fe = [] #interface molefraction of Fe
C_i_Fe = [] #interface concentration of Fe
X_i_Cr = [] #interface molefraction of Cr
C_i_Cr = [] #interface concentration of Cr
x_Fe = [] #bulk molefraction of Fe
c_Fe = [] #bulk concentration of Fe
x_Cr = [] #bulk molefraction of Cr
c_Cr = [] #bulk concentration of Cr
x_O = [] #bulk molefraction of O
c_O = [] #bulk concentration of O
x_FeO = [] #bulk molefraction of FeO
c_FeO = [] #bulk concentration of FeO
x_Cr2O3 = [] #bulk molefraction of Cr2O3
c_Cr2O3 = [] #bulk concentration of Cr2O3
w_percentFe = [] #weight percent of Fe
w_percentCr = [] #weight percent of Cr
w_percentO = [] #weight percent of O
w_percentFeO = [] #weight percent of FeO
w_percentCr2O3 = [] #weight percent of Cr2O3
sys_Temp = [] #temperature of the system

#initializing metal phase data
mass_metal.append(mmass_initial) #mass of metal = 150 tons
#print('mass_metal',mass_metal[0])
volM.append(mass_metal[0]/mDen) #volume of metal in m^3
#print('volM',volM[0])
mass_Fe.append(Fe_initial_percent*mmass_initial/100) #mass of Fe in kg
#print('mass_Fe',mass_Fe[0])
mass_Cr.append(Cr_initial_percent*mmass_initial/100) #mass of Cr in kg
#print('mass_Cr',mass_Cr[0])
kmolFe.append(mass_Fe[0]/56) #kilomoles of Fe
#print('kmolFe',kmolFe[0])
kmolCr.append(mass_Cr[0]/52) #kilomoles of Cr
#print('kmolCr',kmolCr[0])
molOaddedpersec = 2*total_Oflow*101325/(8.314*298*60)# moles of O added to metal per sec
#print('molOaddedpersec',molOaddedpersec)
mass_O.append(O_initial_percent*mmass_initial/100 + molOaddedpersec*16*timestep/1000) #mass of O in kg
#print('mass_O',mass_O[0])
kmolO.append(mass_O[0]/16) #kilomoles of O
#print('kmolO',kmolO[0])
kmolM.append(kmolO[0]+kmolCr[0]+kmolFe[0]) #kilomoles of metal
#print('kmolM',kmolM[0])
molvolM = kmolM[0]*1000/(volM[0]) #molar volume of metal in mol/m^3
#print('molvolM',molvolM)
w_percentFe.append(mass_Fe[0]*100/mass_metal[0])
#print('w_percentFe',w_percentFe[0])
w_percentCr.append(mass_Cr[0]*100/mass_metal[0])
#print('w_percentCr',w_percentCr[0])
w_percentO.append(mass_O[0]*100/mass_metal[0])
#print('w_percentO',w_percentO[0])


#initializing about slag phase
mass_slag.append(smass_initial) #mass of slag = 100kg
#print('mass_slag', mass_slag[0])
volS.append(mass_slag[0]/sDen) #volume of slag in m^3
#print('volS', volS[0])
mass_Cr2O3.append(Cr2O3_initial*smass_initial/100) #mass of Cr2O3 in kg
#print('mass_Cr2O3', mass_Cr2O3[0])
mass_FeO.append(FeO_initial*smass_initial/100) #mass of FeO in kg
#print('mass_FeO',mass_FeO[0])
kmolCr2O3.append(mass_Cr2O3[0]/(2*52 + 3*16)) #kilomoles of Cr2O3
#print('kmolCr2O3',kmolCr2O3[0])
kmolFeO.append(mass_FeO[0]/(56 + 16)) #kilomoles of FeO
#print('kmolFeO',kmolFeO[0])
kmolS.append(kmolCr2O3[0]+kmolFeO[0]) #kilomoles of slag
#print('kmolS',kmolS[0])
molvolS = kmolS[0]*1000/volS[0] #molar volume of slag
#print('molvolS',molvolS)
w_percentFeO.append(mass_FeO[0]*100/mass_slag[0])
#print('w_percentFeO',w_percentFeO[0])
w_percentCr2O3.append(mass_Cr2O3[0]*100/mass_slag[0])
#print('w_percentCr2O3',w_percentCr2O3[0])

# mass transfer model
for i in range(164):#(int(t_sec/timestep)):
    # updating bulk molefractions of all species and temperature
    x_Fe.append(kmolFe[i]/kmolM[i])
    x_Cr.append(kmolCr[i]/kmolM[i])
    x_O.append(kmolO[i]/kmolM[i])
    x_FeO.append(kmolFeO[i]/kmolS[i])
    x_Cr2O3.append(kmolCr2O3[i]/kmolS[i])
    sys_Temp.append(temp_initial+0.2*i*timestep)
    #print(sys_Temp[i])
    #print('kmolFe[i]',kmolFe[i],'kmolCr[i]',kmolCr[i],'kmolO[i]',kmolO[i])
    #print('x_FeO',x_FeO[i],'\nx_Cr2O3',x_Cr2O3[i], '\nx_O',x_O[i], '\nx_Fe', x_Fe[i], '\nx_Cr', x_Cr[i])


    # solving for interface mol-fractions of Fe and Cr
    deltaG_Fe = -121009.9 + 53.114*sys_Temp[i] + 8.314*sys_Temp[i]*log(0.5585/(16))
    #print('deltaG_Fe',deltaG_Fe)
    b = exp(-deltaG_Fe/(8.314*sys_Temp[i]))
    #print('eq constant Fe',b)
    X_i_Fe.append((x_FeO[i])/(x_O[i]*b)) #got mol-fraction of Fe at interface
    #print('X_i_Fe[i]',X_i_Fe[i])
    deltaG_Cr = -274347 + 120.55*sys_Temp[i] + 8.314*sys_Temp[i]*log(0.5585/(16))
    #print('deltaG_Cr',deltaG_Cr)
    b = (exp(-deltaG_Cr/(8.314*sys_Temp[i])))
    #print('eq constant Cr',b)
    X_i_Cr.append(((x_Cr2O3[i]**(1/3))/(x_O[i]*b))**(3/2)) #got mol-fraction of Cr at interface
    #print('X_i_Cr[i]',X_i_Cr[i])

    #finding concentrations from interface molefractions
    c_Fe.append(x_Fe[i]*molvolM)
    C_i_Fe.append(X_i_Fe[i]*molvolM)
    c_Cr.append(x_Cr[i]*molvolM)
    C_i_Cr.append(X_i_Cr[i]*molvolM)
    c_O.append(x_O[i]*molvolM)
    c_FeO.append(x_FeO[i]*molvolS)
    c_Cr2O3.append(x_Cr2O3[i]*molvolS)
    #print('c_Fe',c_Fe[i],'\nC_i_Fe',C_i_Fe[i],'\nc_Cr',c_Cr[i],'\nC_i_Cr',C_i_Cr[i],'\nc_O',c_O[i],'\nc_FeO',c_FeO[i],'\nc_Cr2O3',c_Cr2O3[i])

    #solving flux
    JFe = -kFe*(c_Fe[i]-C_i_Fe[i])
    JCr = -kCr*(c_Cr[i]-C_i_Cr[i])
    JFeO = -JFe
    JCr2O3 = -(1/2)*JCr
    JO = JFe+(3/2)*JCr

    #finding new kilo-moles of species for next loop
    kmolFe.append(JFe*Area*timestep/1000 + kmolFe[i])
    #print('kmol of Fe added',JFe*Area*timestep/1000)
    kmolCr.append(JCr*Area*timestep/1000 + kmolCr[i])
    #print('kmol of Cr added',JCr*Area*timestep/1000)
    kmolO.append(JO*Area*timestep/1000 + kmolO[i] + timestep*molOaddedpersec/1000)
    #print('kmol of O added',JO*Area*timestep/1000)
    #print('kmol of O added from outside',timestep*molOaddedpersec/1000)
    kmolFeO.append(JFeO*Area*timestep/1000 + kmolFeO[i])
    #print('kmol of FeO added',JFeO*Area*timestep/1000)
    kmolCr2O3.append(JCr2O3*Area*timestep/1000 + kmolCr2O3[i])
    #print('kmol of Cr2O3 added',JCr2O3*Area*timestep/1000)
    kmolM.append(kmolFe[i+1]+kmolCr[i+1]+kmolO[i+1])
    kmolS.append(kmolFeO[i+1]+kmolCr2O3[i+1])
    #print('new kmolFe',kmolFe[i+1],'\nnew kmolCr',kmolCr[i+1],'\nnew kmolO',kmolO[i+1],'\nnew kmolFeO',kmolFeO[i+1],'\nnew kmolCr2O3',kmolCr2O3[i+1],'\nnew kmolM',kmolM[i+1],'\nnew kmolS',kmolS[i+1])
    mass_Fe.append(kmolFe[i+1]*56)
    mass_Cr.append(kmolCr[i+1]*52)
    mass_O.append(kmolO[i+1]*16)
    mass_FeO.append(kmolFeO[i+1]*(16+56))
    mass_Cr2O3.append(kmolCr2O3[i+1]*(3*16+2*52))
    mass_metal.append(mass_Fe[i+1]+mass_Cr[i+1]+mass_O[i+1])
    mass_slag.append(mass_FeO[i+1]+mass_Cr2O3[i+1])
    #print('new metal mass',mass_metal[i+1],'\nnew slag mass',mass_slag[i+1])
    #print('mass_Fe',mass_Fe,'mass_Cr',mass_Cr,'kmolO',mass_O,'kmolFeO',mass_FeO,'kmolCr2O3',mass_Cr2O3,'kmolM',mass_metal,'kmolS',mass_slag)
    w_percentFe.append(mass_Fe[i+1]*100/mass_metal[i+1])
    w_percentCr.append(mass_Cr[i+1]*100/mass_metal[i+1])
    w_percentO.append(mass_O[i+1]*100/mass_metal[i+1])
    w_percentFeO.append(mass_FeO[i+1]*100/mass_slag[i+1])
    w_percentCr2O3.append(mass_Cr2O3[i+1]*100/mass_slag[i+1])
    #print('w_percentFe',w_percentFe[i+1],'\nw_percentCr',w_percentCr[i+1],'\nw_percentO',w_percentO[i+1],'\nw_percentFeO',w_percentFeO[i+1],'\nw_percentCr2O3',w_percentCr2O3[i+1])


plt.plot(t, w_percentFe, 'g-', label='Fe mass %')
plt.plot(t, w_percentCr, 'r-', label='Cr mass %')
plt.plot(t, w_percentO, 'c-', label='O mass %')
plt.plot(t, w_percentFeO, 'k-', label='FeO mass %')
plt.plot(t, w_percentCr2O3, 'b-', label='Cr2O3 mass %')
plt.legend(loc='best', bbox_to_anchor=(0.35, -0.05, 0.5, 0.5), prop = {'size': 5})
plt.xlabel('time in seconds')
plt.ylabel('mass percent')

ax = plt.gca()
ax.set_facecolor('y')

fig = plt.gcf()
fig.set_size_inches(4.3, 5)
fig.show()
plt.show()

#making all plots
#plot1 = plt.figure(1)
#plt.plot(t, mass_metal, label="Mass of Metal")
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass of metal vs time")

#plot2 = plt.figure(2)
#plt.plot(t, mass_slag, label="Mass of Slag")
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass of slag vs time")

#plot3 = plt.figure(3)
#plt.plot(t, w_percentFe)
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass percent of Fe vs time")

#plot4 = plt.figure(4)
#plt.plot(t, w_percentCr)
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass percent of Cr vs time")

#plot5 = plt.figure(5)
#plt.plot(t, w_percentO)
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass percent of O vs time")

#plot6 = plt.figure(6)
#plt.plot(t, w_percentFeO)
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass percent of FeO vs time")

#plot7 = plt.figure(7)
#plt.plot(t, w_percentCr2O3)
#plt.xlabel("time in seconds")
#plt.ylabel("Mass in kg")
#plt.title("Mass percent of Cr2O3 vs time")
#t.pop()

#plot8 = plt.figure(8)
#plt.plot(t, c_FeO)
#plt.xlabel("time in seconds")
#plt.ylabel("Conc in mol/m^3")
#plt.title("Concentration of FeO vs time")

#plot9 = plt.figure(9)
#plt.plot(t, c_Cr2O3)
#plt.xlabel("time in seconds")
#plt.ylabel("Conc in mol/m^3")
#plt.title("Concentration of Cr2O3 vs time")





