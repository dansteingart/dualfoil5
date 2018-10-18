from pithy import *
import dualfoil

df = dualfoil.dualfoil("user")


#dualfoil script

#clear 

df.clear_cycles()
df.set_ocv(3) #set ocv for 3 minutes
df.set_current(15,120) #discharge at 15 A/m^2 for 120 minutes
df.set_ocv(10) #set ocv for 10 minutes
df.set_current(-6,10) #charge for 10 minutes
#df.set_potential(4,10) #set to 4V for 10 minutes 


pp = df.parts.keys()
pp.sort()
#for p in pp: print p,df.parts[p]
df.parts['ep1'] = .2 #set 
df.parts['h1'] = 500e-6
df.parts['h3'] = 500e-6
#df.parts['tmmax'] = 1
df.cycles = df.cycles.strip()
df.writeOut()
df.runDualFoil(debug=True)
data = df.readOutput()

df.ivtplot()

#df.surfplot()


profs = df.readProfiles()

print "Plot the current distribution"
df.surfplot(profs,"Distance (um)","Liq Cur (A/m^2)")


print "Plot the electrolyte conentration"
df.surfplot(profs,"Distance (um)","C Elec (mol/m3)")

print "Plot the solid conentration"
df.surfplot(profs,"Distance (um)","C Sol Surf (x/y)")


