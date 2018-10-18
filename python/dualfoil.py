##OK TO COPY
from pithy import *
from commands import getoutput as go
import subprocess
import atexit
import os
from glob import glob
import signal
from StringIO import StringIO
import pandas as pd
from scipy.interpolate import interp1d
import filecmp
from scipy.interpolate import interp2d as i2d

class dualfoil():
    def __init__(self,name,dfdir="/dualfoil5/df5.1/"):
        self.user = name
        self.df = dfdir
        self.dualbase = dfdir+"dualfoil5.in"
        self.parts = self.readin()
        self.cycles = self.getcycles()

    def clear_cycles(self): self.cycles = ""

    def set_ocv(self,mins): 
        self.cycles += "0 %f 1 2.0 4.70\n" % mins
    
    def set_current(self,cur,mins,minv=2.0,maxv=4.7): 
        self.cycles += "%f %f 1 %f %f\n" % (cur,mins,minv,maxv)




    def readin(self):
        data = open(self.dualbase).read()
        parts = {}
        lines = data.split("\n")
        for l in lines[0:74]:
            p = l.split("!")
            try:
                key = p[1].split(",")[0].strip().split(" ")[0]
                parts[key] = p[0]
            except Exception as err:
                print err
        go("mkdir df_"+self.user)
        go("cp -n "+self.df+"/* df_"+self.user+"/")
        return parts
    
    def getcycles(self):
        data = open(self.dualbase).read()
        lines = data.split("\n")
        return lines[74]
    
    def writeOut(self,endless=False):
        self.cycles = self.cycles.strip()
        data = open(self.dualbase).read()
        self.parts['lcurs'] = str(len(self.cycles.split('\n')))
        if endless: self.parts['lcurs'] = 101
        lines = data.split("\n")
        pull = self.cycles.replace('\n', 'd', 101).find('\n')
        if pull != -1:
            self.cycles = self.cycles[0:pull]
            print "too many changes, going endless"
            self.parts['lcurs'] = 101
    
        fn =  "df_"+self.user+"/dualfoil5.in"
        fnc =  "df_"+self.user+"/lastin"
        out = ""
        for l in lines[0:74]:
            p = l.split("!")
            key = p[1].split(",")[0].strip().split(" ")[0]
            out += str(self.parts[key]) + " ! "+p[1]+"\n"
        out += self.cycles +"\n"
        try:
            fil = open(fn).read()
            if fil == out:
                doo = "foo"
                #print "same as the old file, not changing"
            else:
                #print "changed!  writing"
                open(fn,"w").write(out)
        except:    
            open(fn,"w").write(out)
    
    def runDualFoil(self,debug = False,filename=None):

        try:
            fni =  "df_"+self.user+"/dualfoil5.in"
            fno =  "df_"+self.user+"/dualfoil5.out"
            ti = os.path.getmtime(fni)
            to = os.path.getmtime(fno)
        except:
            #if something farts above, just make it so we have to 
            #run the simulation
            ti = 1
            to = 0
        
        runs = glob("files/df_%s/*/dualfoil5.in" % self.user)
        same = None
        for run in runs:
            if filecmp.cmp(fni,run) == True:
                same = run
                rd = run.replace("dualfoil5.in","")
                go("cp %s/*.out df_%s/" % (rd,self.user))
                print "same as %s" % run
                print "not running, pulling previous run"
                return 0
        
        if debug: 
            print "Time Difference Between last input and output:", ti-to," s"

        if ti > to:
            if debug: 
                print "input older that output: analysis is not current: running simulation"
            
            a = go("rm df_"+self.user+"/*.out")
            if debug: print a
            b = go("cd df_"+self.user+";./dualfoil")
            if debug: print b
            c = go("mkdir files/"+"df_"+self.user)
            if debug: print c
        
        
            if filename == None: filename = "run_"+str(time.time())
            dd = {}
            dd['user'] = self.user
            dd['name'] = "files/df_"+self.user+"/"+filename
            d = go("mkdir -p %(name)s; cp df_%(user)s/*.out %(name)s/; cp df_%(user)s/*.in %(name)s/" % dd) 
            if debug: print d
        else: 
            doo = "foo"
            if debug: 
                print "input younger that output: analysis is current: no need to run simulation"
    
    def readProfiles(self,filename=None):
        if filename == None:
            fn = "df_%s/profiles.out" % self.user
        f = open(fn).read()

        Header = "Distance (um),C Elec (mol/m3),C Sol Surf (x/y),Liq Pot (V),Solid Pot (V),Liq Cur (A/m^2),i main (A/m^2),j side 1 (A/m^2),j side 2 (A/m^2), j side 3 (A/m^2)"
        
        #print f
        data = f.split("\n  \n  \n")
        
        profiles = {}
        
        for i in range(len(data)):
            try:
                d = data[i].split("\n")
                step_time = float(d[3].split("=")[-1].replace("min",""))
                profs = Header+"\n"
                for l in d[4:]: profs += l+"\n"
                df = pd.read_csv(StringIO(profs))
                profiles[step_time] = df
            except: weare = "moving on"
        return profiles

    def readOutput(self,debug = False,filename=None):
        #these are the values we care about
        header = ['t','nutil','putil','vcell',"ocp",'i','T','Q']

        #If we don't ask for a specific run, see if there's a recent run
        if filename == None:
            try: data = open("df_"+self.user+"/dualfoil5.out").read()
            #if there's not a recent run, no biggie.
            except: data = ""
        else:
            try:
                data = open("files/df_"+self.user+"/dualfoil_"+filename+".out").read()
            except:
                print "going to video tape"
                print filename
                data = open(filename).read()
    
    
        start = data.find("DUAL INSERTION CELL")
        rest = data
        data = data[start:]
        data = data.split("\n")
        out = {}
        for h in header:
            out[h] = []
        try:
            data.pop(0)
            data.pop(0)
            data.pop(0)
            data.pop(0)
        except:
            print "can't pop"
        for d in data:
            p = d.strip().split(",")   
            for i in range(0,len(header)):
                try:
                    out[header[i]].append(float(p[i]))
                except Exception as Err:
                    if debug: print Err
                    nope = True 
        for h in header:
            out[h] = array(out[h])
    
        sestart = rest.find("specific energy segment")
        seend = rest.find("W-h/kg")
        try:
            out['se'] = float(rest[sestart:seend].split("=")[1])
        except Exception as err: 
            out['se'] = -1
        
        spstart = rest.find("specific power segment")
        spend = rest.find("W/kg")
        try:
            out['sp'] = float(rest[spstart:spend].split("=")[1])
        except:
                out['sp'] = -1
        return out

    def surfplot(self,profs,xval,yval,dpi=150):
            times = profs.keys()
            times.sort()
            app = []
            for t in times:
                x = profs[t][xval]
                y = profs[t][yval]
                x = array(x)
                y = array(y)
                app.append(y)
                
            
            app = array(app)
            times = array(times)
            z = i2d(x,times,app)
            xx = linspace(0,max(x),1000)
            yy = linspace(0,max(times),1000)
            zz = z(xx,yy)
            
            
            imshow(zz,
                    extent=(0,max(x),max(times),0),
                    aspect='auto',
                    cmap=cm.jet)
            colorbar()
            title(yval)
            xlabel("Position ($\mu m$)")
            ylabel("Time (minutes)")
            showme(dpi=dpi)
            clf()

    def ivtplot(self,dpi=150):
        data = self.readOutput()
        #plot stuff
        subplot(2,1,1)
        plot(data['t'],data['vcell'],'k',label="loaded potential")
        plot(data['t'],data['ocp'],'grey',label='eq potential',linewidth=.5)
        legend(loc='best',fontsize=8)
        xticks([])
        ylabel("Potential (V)")
        subplot(2,1,2)
        xlabel("Time (m)")
        plot(data['t'],data['i'],'k')
        ylabel("Current ($mA/cm^2$)")
        xlabel("Time (m)")
        showme(dpi=dpi)
        clf()


if __name__ == "__main__":
    
    dpi = 300 #pretty pictures
    
    df = dualfoil("user")
    
    #dualfoil script
    
    #clear 

    df.clear_cycles()
    df.set_ocv(3) #set ocv for 3 minutes
    df.set_current(15,60) #discharge at 15 A/m^2 for 30 minutes
    df.set_ocv(10) #set ocv for 10 minutes
    df.set_current(-6,10) #charge for 10 minutes
    #df.set_potential(4,10) #set to 4V for 10 minutes 

    
    pp = df.parts.keys()
    pp.sort()
    #for p in pp: print p,df.parts[p]
    df.parts['ep1'] = .2 #set 
    df.parts['h1'] = 200e-6
    df.parts['h3'] = 200e-6
    #df.parts['tmmax'] = 1
    df.cycles = df.cycles.strip()
    df.writeOut()
    df.runDualFoil(debug=True)
    data = df.readOutput()


    #plot stuff
    subplot(2,1,1)
    plot(data['t'],data['vcell'],label="loaded potential")
    plot(data['t'],data['ocp'],'grey',label='eq potential',linewidth=.5)
    legend(loc='best',fontsize=8)
    xticks([])
    ylabel("Potential (V)")
    subplot(2,1,2)
    xlabel("Time (m)")
    plot(data['t'],data['i'])
    ylabel("Current ($mA/cm^2$)")
    xlabel("Time (m)")
    showme(dpi=dpi)
    clf()
    
    profs = df.readProfiles()
    
    times = profs.keys()
    times.sort()

    #for coloration
    count = 0.0
    app = []
    for t in times:
        fac = count*.02
        #coloring trick
        f = count/len(times)
        c = (f,0,1-f)  
                #make perty
        x = profs[t]['Distance (um)']
        y = profs[t]['C Sol Surf (x/y)']
        x = array(x)
        y = array(y)
        y[y == 0] = 'nan' # or use np.nan
        annotate("%1.f min" % t,xy=(x[-1],y[-1]+fac))
        plot(x,y+fac,color=c)
        count +=1
    
    yticks([])
    ylabel("Solid Lithium Fraction")
    xlabel("Position ($\mu m$)")
    showme(dpi=dpi)
    clf()
    

    count = 0.0
    app = []
    for t in times:
        #coloring trick
        fac = count+1000
        f = count/len(times)
        c = (f,0,1-f)  
        #make perty
        x = profs[t]['Distance (um)']
        y = profs[t]['C Elec (mol/m3)']
        x = array(x)
        y = array(y)
        annotate("%1.f min" % t,xy=(x[-1],y[-1]+fac))
        plot(x,y+fac,color=c)
        count +=1
    
    yticks([])
    ylabel("Electrolyte Lithium Concentration")
    xlabel("Position ($\mu m$)")
    showme(dpi=dpi)
    clf()##OK TO COPY
