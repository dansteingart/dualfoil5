

                                                                     
                                                                     
                                                                     
                                             
                                                                     
                                                                     
                                                                     
                                             
c***********************************************************************
c     dualfoil.f  (version 5.1)  January 31, 2009

c     Dual lithium ion insertion cell
c  Copyright Marc Doyle and John Newman 1998.
c  You may make a copy of this program which you may
c  personally and freely use in its unaltered form.
c  You may distribute this program subject to the
c  conditions that it be made freely available and
c  that any duplication of this program must be
c  essentially unaltered and must include this notice.
c
c  We make no warranties, express or implied, that
c  this program is free of errors or that it will
c  meet the requirements of your application.  The
c  author and publisher disclaim all liablility for
c  direct or consequential damages resulting from
c  use of this program.
c
c  Note: For lflag=0, the model works only for initially zero current.
c
c***********************************************************************
c Based on the original articles:
c Marc Doyle, Thomas F. Fuller, and John Newman,
c J. Electrochem. Soc., 140 (1993), 1526-1533.
c Thomas F. Fuller, Marc Doyle, and John Newman,
c J. Electrochem. Soc., 141 (1994), 1-10.
c Thomas F. Fuller, Marc Doyle, and John Newman,
c J. Electrochem. Soc., 141 (1994), 982-990.

c     This program was developed for lithium and lithium-ion
c     batteries, with the flexibility that a number of positive
c     and negative active materials of the insertion or intercalation
c     type can be selected or added, in addition to a lithium foil
c     electrode.  Several electrolytes are also provided, and
c     these can be added to by providing properties in subroutine
c     prop.  An aqueous system, such as a nickel/metal hydride
c     battery in a KOH solution, could also be added, since it
c     follows the pattern of being a dual-insertion system.  Primary
c     cells, as well as rechargeable cells, can be simulated with
c     this program.

c     This program should be very similar to Marc
c     Doyle's impedance program of 1995 to 1998.
c     See Doyle, Meyers, Darling, Newman, JES, 147 (2000),
c     99-110 and 147 (2000), 2930-2940.
c     The program has been modified to include a
c     zero time step in order to get the potentials
c     and currents correct and to have the impedance
c     independent of the initial values set in
c     subroutine guess.  It should run even with
c     a nonzero current density, but this is not
c     recommended and is not proper in the sense
c     of having a steady state as a basis for the
c     impedance simulation (one of the criteria for
c     the validity of the Kramers-Kronig relations).
c     The program has been cleaned up a lot by
c     rearranging equations ii2, i2div, and kp1 and
c     simplifying the impedance statements.
c
c     Includes FOIL and impedance option 7-1-95
c     Some modifications made to Bellcore stuff on 1-12-96
c     Error in EIS b.c.'s for eqns i2div and ip1 fixed 1-23-96
c     Some additional changes up to 7/96 incl. film resistance
c
c     Modified impedance section by Newman on 4/4/98
c     HALF-CELL version:  split off from impnew1.f on 4/19/98
c     Current newest version June 29, 1998
c     Has Jeremy Meyers modifications 6/98
c
c     Revised June, 1998, to include double-layer capacitance in
c     each electrode and to correct a factor of two in Ohm's law.
c     See J. Electrochem. Soc., 146 (1999), 4360-4365.
c
c     Revised by Karen Thomas, Feb. 12, 1999:
c     - if n1 = 0, then code treats the negative electrode as metal foil.
c     - To run, simply type "dualfoil", then enter
c       input and output file names when prompted.
c     - double-layer capacitance is not currently calculated at
c       a foil electrode
c     July 16, 2000:
c     -  Program modified to be compatible with g77 fortran compiler.
c     July 6, 2001:
c     - Prints total and irreversible heat generation; positive heat
c       is exothermic
c     August 2002:
c     - corrected double precision, added graphite, LiCoO2, and LiNi0.8Co0.2O2
c     May 2003:
c     - made compatible with Compaq visual fortran compiler, other minor changes
c     April, 2005:
c      -The version of Dualfoil dated (July 15, 2003) used a single cutoff
c      potential, which caused problems when a user attempted to charge the cell
c      following a discharge.  For example, if the global cutoff potential was
c      set to 3.0 V and the cell discharged beneath that value, upon charging the
c      cutoff potential would be triggered when the cell potential exceeded 3.0 V.
c      The revised version of dualfoil.f uses two global cutoff potentials,
c      one that serves as an upper bound and one that serves as a lower bound.
c      dualfoil.in has been modified to allow the user to specify two cutoff
c      potentials.
c     July, 2007:
c      -A series of major modifications have been made since the April, 2005 version.
c       See the guide posted on this page for more information.  Note that
c       all previous functionality is maintained by the current version.  Following is
c       an abbreviated list of the changes that have been made:
c        - Addition of specified load and power modes
c        - Ability to specify an upper and lower cutoff potential for each segment
c        - Addition of NiMH materials (MH anode, NiOOH cathode, KOH electrolyte)
c        - Addition of side reactions to model gas-phase recombination in NiMH system
c          or side reactions in another system (e.g., solvent reduction in a Li-ion cell).
c        - Ability to use a variable solid-phase diffusion coefficient in a
c          "pseudo-2d" arrangement
c        - Use of a new input file with activation energies for temperature dependences
c          of the transport properties (solid and electrolyte diffusion, electrolyte
c          conductivity), rate constants, and film resistances.  Where appropriate
c          there is a separate activation energy in each electrode.  The activation
c          energies are stored in the input file ebar.in
c        - Detailed profiles are now written as the program goes along (rather
c          than after the code has finished running) to a separate file, called
c          profiles.out, if the user selects il1=1
c        - A problem with the identity/capacity/density of the cathode in the
c          default input file has been corrected
c        - The unused subroutine util has been removed
c        - An additional equation has been added to improve convergence
c          when running constant potential, power, and load.
c        - The ability to have up to five particles per electrode has been added
c        - The variables and equations have been given names rather than numbers
c        - The energy balance is now done according to the method given
c          in Rao and Newman, J. Electrochem. Soc., 144 (1997), 2697-2704.
c          Cp(dT/dt)-Q=-Integral(sum over reactions (a*in,i*Uh,i))dx-IV
c      February, 2008:
c        The following changes have been made in going from version 5.0 to 5.1
c        - Subroutine guess has been removed
c        - Convergence has been improved and the dV/dt trigger has been removed
c        - fj(j) for cutting off part of the electrode has been removed
c        - The variables have been renamed
c        - A grid resistance, RG, has been added
c        - An error in the impedance calculation for anode and cathode has
c          been corrected; the original full-cell impedance was correct
c        - The equation structure for impedance mode has been made general,
c          facilitating future additions
c        - We have improved the output formatting
c        - We have improved the logic of the multiple-particle selection section

c***********************************************************************
      implicit real*8(a-h,o-z)
      character *30 filin, filout
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/power/ ed,Vold,ranodesave,rcathdesave,heat,qlosstot
      common/ssblock/ xp0(16),xx0(16,221),term(221)
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/imp/ vreal(110),vimag(110),phi(110),zmag(110),omega(510),
     1vreala(110),vrealc(110),vimaga(110),vimagc(110),
     1phia(110),phic(110),zmaga(110),zmagc(110)
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/tprop/df(221),cd(221),tm(221),ddo2(221),ddh2(221),
     1ddf(221),dcd(221),dtm(221),dfu(221),d2fu(221),do2(221),dh2(221)
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
      common/side/rksc1,c1init,c2init,rksa1,term_s1(221),vol,
     &rksa2,rksa3,rksc3,rksc2,UsO2,UsH2,cn2,term_s2(221),nside
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      common/gas/ epg1,epg2,epg3
      common /maxpow/ lcount
      common/RG/ RG,RGn,RGp,RGext
      dimension terms(221),tt(200),cu(200),tot(200),mc(200),taper(200)
      dimension vcutlo(200),vcuthi(200),RGin(200)
      dimension term_s1s(221)
      dimension rad1paf(5),rad3paf(5)
   44 format(/'            mass  = ',f7.4,' kg/m2')
   45 format(' specific energy whole run  = ',f8.2,' W-h/kg')
   46 format(' specific power whole run  = ',f8.2,' W/kg')
   47 format(' total heat = ', f15.4,' W-h/m2')
   48 format(' specific energy segment = ', f8.2,' W-h/kg')
   49 format(' specific power segment  = ', f8.2,' W/kg')
   57 format(e15.6,', ',e15.6,', ',e15.6)
   89 format(f8.1,',',f8.1,',',f8.1,',',i3,',',g10.5,',',f5.2)
   90 format (g10.3,',')

c
      open(3,file='halfcells.out',status='unknown')
      write (3,*) ' Time (min)    V neg      V pos       Cur (A/m2)   '
     &,' Temp (C)   Heat Gen (W/m3)'
c     print *, 'Enter input file name, press return'
c     read *, filin
c     open (1, FILE = 'dualfoil5.in', status = 'old')
      open (1, FILE = 'li-ion.in', status = 'old')
c     open (1, FILE = '2Rload=0.0003ohm.m2.in', status = 'old')
c     open (1, FILE = filin, status = 'old')
c     print *, 'Enter output file name, press return'
c     read *, filout
      open (2, file = 'dualfoil5.out', status = 'unknown')
c     open (2, file = filout, status = 'unknown')
      open (4,file='profiles.out',status='unknown') ! for detailed profiles
      open (5,file='li-ion-ebar.in', status='old') ! activation energies
      open (6,file='resistances.out',status='unknown') ! for cell resistance
      open (7,file='solidprof.out',status='unknown') ! for profile inside solid
      open (8,file='Ds',status='unknown') ! for Diffusion cofficient
c
      data fc/96487.0d0/, r/8.314d0/, pi/3.141592653589d0/
      data ed/0.d0/, Vold/0.d0/

      npa=2 ! Number of particle sizes
      if (npa.ne.1.and.imp.eq.1) print *,'npa must be 1 for imp=1'
      
c     Variable #      Description
c     1               Solid phase,matrix 
c     2               Liquid phase
      
      ii2=3+npa ! equation number for current density i2
      ki2=3+npa ! variable number for current density i2
      i2div=4+npa ! equation number for current balance.
      kj=4+npa ! variable number for transfer current (flux density).
      ip1=5+npa ! equation number for Ohm's law in the matrix.
      kp1=5+npa ! variable number for PHI_1, matrix potential.
      ip2=2 ! equation number for Ohm's law in solution.
      kp2=2 ! variable number for PHI_2, solution potential.
      
      imb1=6+npa ! equation number (side reaction 1 matl balance)
      kS1=6+npa ! variable number (side reaction 1 matl balance)
      iSk1=7+npa ! equation number (side reaction 1 kinetics)
      kj1=7+npa ! variable number (side reaction 1 kinetics)
      iSk2=8+npa ! equation number (side reaction 2 kinetics)
      kj2=8+npa ! variable number (side reaction 2 kinetics)
      imb2=9+npa ! equation number (side reaction 2 matl balance)
      kS2=9+npa ! variable number (side reaction 2 matl balance)
      iSk3=10+npa ! equation number (side reaction 3 kinetics)
      kj3=10+npa ! variable number (side reaction 3 kinetics)
      lim2=100  ! limit for iterations to converge load/power/poten


c***********************************************************************
c     read in parameters and boundary conditions

c
      read (1,*) lim !limit on number of iterations
      read (1,*) h1  !thickness of negative electrode (m)
      read (1,*) h2  !thickness of separator (m)
      read (1,*) h3  !thickness of positive electrode (m)
      read (1,*) hcn !thickness of negative electrode current collector (m)
      read (1,*) hcp !thickness of positive electrode current collector (m)
      thk=h1+h2+h3
      read (1,*) n1  !number of nodes in negative electrode
c     If negative electrode is metal foil, let n1 = 0
      read (1,*) n2  !number of nodes in separator
      read (1,*) n3  !number of nodes in positive electrode
      read (1,*) n4  !number of nodes in solid particle
      read (1,*) mvdc1 !flag for variable solid diff coeff, anode
      read (1,*) mvdc3 !flag for variable solid diff coeff, cathode
      read (1,*) lims  !limit on number of iterations in solid phase
      read (1,*) t   !temperature (K)
      write (2,1101) lim,1.d6*h1,1.d6*h2,1.d6*h3,1.d6*hcn,1.d6*hcp
     &,n1,n2,n3,t
      n2=n2+1
      nj=n1+n2+n3 !total number of nodes
      nnj=n4 !number of nodes for particle radius
c
      read (1,*) xx(1,n1+2) ! initial electrolyte concentration (mol/m3)
      read (1,*) csx  !initial stochiometric parameter for negative particle
      read (1,*) csy  !initial stochiometric parameter for positive particle
      read (1,*) tmmax !maximum time step size (s)
      read (1,*) dfs1 !diffusion coefficient in negative solid (m2/s)
      read (1,*) dfs3 !diffusion coefficient in positive solid (m2/s)
      dfs1save=dfs1 !save the original values because of temp dependence
      dfs3save=dfs3
      read (1,*) Rad1 !radius of negative particles (m)
      read (1,*) Rad3 !radius of positive particles (m)
      write (2,1102) xx(1,n1+2),csx,csy,tmmax,dfs1,dfs3,
     &1.d6*Rad1,1.d6*Rad3
c     If negative electrode is metal foil, let ep1=epp1=epf1=0.0
      read (1,*) ep1  !volume fraction of electrolyte in negative electrode
      read (1,*) epp1 !volume fraction of polymer phase in negative electrode
      read (1,*) epf1 !volume fraction of inert filler in negative electrode
      read (1,*) epg1 !volume fraction of gas in negative
      read (1,*) ep2  !volume fraction of electrolyte in separator
      read (1,*) epp2 !volume fraction of polymer phase in separator
      read (1,*) epg2 !volume fraction gas in separator
      read (1,*) ep3  !volume fraction of electrolyte in positive electrode
      read (1,*) epp3 !volume fraction of polymer phase in positive electrode
      read (1,*) epf3 !volume fraction of inert filler in positive electrode
      read (1,*) epg3 !volume fraction of gas in positive
      read (1,*) sig1 !conductivity of solid negative matrix (S/m)
      read (1,*) sig3 !conductivity of solid positive matrix (S/m)
      read (1,*) rka1 !reaction rate constant for negative insertion reaction
      read (1,*) rka3 !reaction rate constant for positive insertion reaction
      rka1save=rka1
      rka3save=rka3
      read (1,*) ranode  !anode film resistance (out of place)
      read (1,*) rcathde !cathode film resistance (out of place)
      ranodesave=ranode
      rcathdesave=rcathde
      read (1,*) cot1 !coulombic capacity of negative material (mAh/g)
      read (1,*) cot3 !coulombic capacity of positive material (mAh/g)
      write (2,1103) ep1,epp1,epf1,ep2,epp2,ep3,epp3,epf3,sig1,
     & sig3,cot1, cot3,rka1,rka3
      read (1,*) re  ! density of electrolyte (kg/m3)
      read (1,*) rs1 ! density of negative insertion material (kg/m3)
      read (1,*) rs3 ! density of positive insertion material (kg/m3)
      read (1,*) rf  ! density of inert filler (kg/m3)
      read (1,*) rpl ! density of polymer phase (kg/m3)
      read (1,*) rc  ! density of separator material (kg/m3)
      read (1,*) rcn ! density of negative current collector (kg/m3)
      read (1,*) rcp ! density of positive current collector (kg/m3)
      write (2,1104) re,rs1,rs3,rf,rpl,rc,rcn,rcp
      read (1,*) htc  !heat-transfer coefficient with external medium (W/m2K)
      read (1,*) Cp   !heat capacity of cell (J/kgK)
      read (1,*) Tam  !ambient temperature (K)
      read (1,*) ncell !number of cells in a cell stack
      read (1,*) lht  !0 uses htc, 1 calcs htc, 2 isothermal
      read (1,*) il1   !1 for long print-out 0 for short print-out
      read (1,*) il2   !1/il2 = fraction of nodes in long print-out
      read (1,*) il3   !1/il3 = fraction of time steps in long print-out
      read (1,*) lflag ! 0 for electrolyte in separator only, 1 for uniform
      read (1,*) imp  ! 0 for no impedance, 1 for impedance
      read (1,*) capp1 ! capacitance of negative material (F/m2)
      read (1,*) capp3 ! capacitance of positive material (F/m2)
      read (1,*) lpow  ! 0 for no power peaks, 1 for power peaks
      read (1,*) jsol  ! calculate solid profiles if 1<jsol<nj
      read (1,*) nside ! flag for side reaction
      read (1,*) rksa1 ! rate constant side reaction anode
      read (1,*) rksc1 ! rate constant side reaction cathode
      read (1,*) rksa2 ! rate constant side reaction anode
      read (1,*) rksc2 ! rate constant side reaction cathode
      read (1,*) rksa3 ! rate constant side reaction anode
      read (1,*) rksc3 ! rate constant side reaction cathode
      c1init=0.1d0   ! init conc of side reaction species 1 (mol/m^3)
      c2init=3.0d0     ! init conc of side reaction species 2 (mol/m^3)
      cn2=45.0d0       ! background N2 concentration for NiMH (mol/m^3)
      residm=0.0d0     ! Residual mass not in cell sandwich (kg/m^2)
      write (2,1105) ranode,rcathde,htc,dudt,Cp,residm,tam,ncell,lht
      write (2,1116) nside
      if (nside.ge.1) then
      write (2,1111) rksa1,rksc1,rksa2,rksc2,rksa3,rksc3
      endif
      read (1,*) nneg  ! designates negative electrode system
      read (1,*) nprop ! designates electrolyte system
      read (1,*) npos  ! designates positive electrode system
      read (1,*) lcurs ! number of current changes
      write (2,1106) il1,il2,il3,lflag,imp,capp1,capp3,lcurs
      read (5,*) EbarS1  !activation energy, anode solid diffusion
      read (5,*) EbarS3  !activiation energy, cathode solid diffusion
      read (5,*) Ebarkap !activation energy electrolyte conductivity
      read (5,*) EbarD   !activation energy electrolyte diffusion
      read (5,*) Ebarka  !activation energy negative kinetics
      read (5,*) Ebarkc  !activation energy positive kinetics
      read (5,*) Ebarks1a  !activation energy O2 side rxn
      read (5,*) Ebarks1c  !activation energy O2 side rxn
      read (5,*) Ebarks2a  !activation energy H2 side rxn
      read (5,*) Ebarks2c  !activation energy H2 side rxn
      read (5,*) Ebarks3a  !activation energy shuttle side rxn
      read (5,*) Ebarks3c  !activation energy shuttle side rxn
      read (5,*) Ebarr1  !activation energy, film resistance anode
      read (5,*) Ebarr3  !activation energy, film resistance cathode

      nn=1 ! used for variable solid-phase diffusion coefficient.
c     n is number of equations, always include equation that 
c     helps converge for load, power, potential
c     side reactions only work with one particle type (npa=1)

      if (nside.eq.0) n=5+npa  !no side reaction, use six eqns
      if (nside.eq.1) n=7+npa  ! first side reaction, mass balance & kinetic
      if (nside.eq.2) n=9+npa  ! second side reaction, mass balance & kinetic
      if (nside.eq.3) n=10+npa  ! third side reaction, only kinetic

      if (nside.eq.1) kj2=kj2+1 !to prevent variable overlap

      do 555 i = 1, lcurs
      read (1,*) cu(i),tt(i),mc(i),vcutlo(i),vcuthi(i),RGin(i)
      
c     Split of RG into its components;Test for first leg. WHT 3/10/14
      RG=RGin(1) !total resistance in foils, leads, and contacts, ohm-m2
      RGn=RG/3.0d0 !resistance affecting negative half cell
      RGp=RG-RGn !resistance affecting positive half cell
      RGext=RG/4.0d0  !resistance outside cell
c     heat produced by RGext donot affect the internal cell temperature

  555 continue
c     cu(i)      operating current density (A/m2) or cell potential (V)
c     tt(i)      time (min) or cutoff potential (V)
c     mc(i) The mode of discharge
c           0 for potentiostatic
c           1 for galvanostatic for a given time
c           2 for galvanostatic to a cutoff potential
c          -1 for galvanostatic for a given time with tapered
c             charge/discharge upon reaching cutoff potential
c          -2 for specified power (in W/m2)
c          -3 for specified load (in ohm-m2)
c     vcutlo(i)  lower cutoff potential
c     vcuthi(i)  upper cutoff potential
c     RGin(i) internal resistance of foils, tabs, etc., Ohm-m2

      print *, 'Now running DUAL...'
 1101 format (i7,'  lim, limit on number of iterations'
     &/1x,f6.2,'  h1,  thickness of negative electrode (microns)'
     &/1x,f6.2,'  h2,  thickness of separator (microns)'
     &/1x,f6.2,'  h3,  thickness of positive electrode (microns)'
     &/1x,f6.2,'  hcn, ',
     &'thickness of negative electrode current collector (microns)'
     &/1x,f6.2,'  hcp, thickness of positive electrode current'
     &,' collector (microns)'
     &/i7,'  n1,  number of nodes in negative electrode'
     &/i7,'  n2,  number of nodes in separator'
     &/i7,'  n3,  number of nodes in positive electrode'
     &/1x,f6.2, '  T,   temperature (K)')
 1102 format (/1x,f6.1,'  xx(1,n1+2), initial concentration (mol/m3)'
     &/1x,f6.4,'  csx,   initial stoichiometric parameter for negative'
     &/1x,f6.4,'  csy,   initial stoichiometric parameter for positive'
     &/1x,f6.1,'  tmmax, maximum time step size (s)'
     &/1x,e6.1,'  dfs1, diffusion coefficient in negative solid (m2/s)'
     &/1x,e6.1,'  dfs3, diffusion coefficient in positive solid (m2/s)'
     &/1x,f6.2,'  Rad1,  radius of negative particles (microns)'
     &/1x,f6.2,'  Rad3,  radius of positive particles (microns)')
 1103 format (/1x,f6.3,'  ep1,'
     &,'   volume fraction of electrolyte in negative electrode'
     &/1x,f6.3,'  epp1,'
     &,'  volume fraction of polymer phase in negative electrode'
     &/1x,f6.3,'  epf1,'
     &,'  volume fraction of inert filler in negative electrode'
     &/1x,f6.3,'  ep2,   volume fraction of electrolyte in separator'
     &/1x,f6.3,'  epp2,  volume fraction of polymer phase in separator'
     &/1x,f6.3,'  ep3,'
     &,'   volume fraction of electrolyte in positive electrode'
     &/1x,f6.3,'  epp3,'
     &,'  volume fraction of polymer phase in positive electrode'
     &/1x,f6.3,'  epf3,'
     &,'  volume fraction of inert filler in positive electrode'
     &/1x,g8.2,'  sig1,  conductivity of negative matrix (S/m)'
     &/1x,g8.2,'  sig3,  conductivity of positive matrix (S/m)'
     &/1x,f7.2,'  cot1,  coulombic capacity of negative material'
     &,' (mAh/g)'
     &/1x,f7.2,'  cot3,  coulombic capacity of positive material'
     &,' (mAh/g)'
     &/1x,e6.1,'  rka1,  reaction rate constant for negative reaction '
     &/1x,e6.1,'  rka3,  reaction rate constant for positive reaction')
 1104 format (/1x,f6.1,'  re,   density of electrolyte (kg/m3)'
     &/1x,f6.1,'  rs1,  density of negative insertion material (kg/m3)'
     &/1x,f6.1,'  rs3,  density of positive insertion material (kg/m3)'
     &/1x,f6.1,'  rf,   density of inert filler (kg/m3)'
     &/1x,f6.1,'  rpl,  density of polymer phase (kg/m3)'
     &/1x,f6.1,'  rc,   density of separator material (kg/m3)'
     &/1x,f6.1,'  rcn,  density of negative current collector (kg/m3)'
     &/1x,f6.1,'  rcp,  density of positive current collector (kg/m3)')
 1105 format (/1x,f10.6,'  ranode,   anode film resistance (ohm-m2)'
     &/1x,f10.6,'  rcathde,  cathode film resistance (ohm-m2)'
     &/1x,f6.2,'  htc,   heat-transfer coefficient with'
     &,' external medium (W/m2K)'
     &/1x,f10.6,'  dUdT,  temperature coefficient of EMF (V/K)'
     &/1x,f6.1,'  Cp,    heat capacity of cell (J/kg-K)'
     &/1x,f6.2,'  residm,    residual mass (kg/m2)'
     &/1x,f6.2,'  Tam,   ambient temperature (K)'
     &/i7,'  ncell, number of cells in a cell stack'
     &/i7,'  lht,   0 uses htc,  1 calcs htc,  2 isothermal')
 1106 format (/i7,'  il1,   1 for long print-out  0 for short print-out'
     &/i7,'  il2,   prints every il2 th node in long print-out'
     &/i7,'  il3,   prints every il3 th time step in long print-out'
     &/i7,'  lflag, 0 for electrolyte in separator only, 1 for uniform'
     &/i7,'  imp, 0 for no impedance, 1 for impedance'
     &/1x,f6.2,'  capp1,  capacitance of negative material'
     &,' (F/m2)'
     &/1x,f6.2,'  capp3,  capacitance of positive material'
     &,' (F/m2)'
     &/i7,'  lcurs, number of current changes')
 1107 format ('   Time     Util N  Util P  Cell Pot   Uocp      Curr',
     &'      pH2     pO2   Total P  Temp   heatgen')
 1108 format ('   (min)       x       y      (V)       (V)      (A/m2)',
     %'   (bar)   (bar)   (bar)   (C)     (W/m2)')
 1109  format ('   Time     Util N  Util P  Cell Pot   Uocp      Curr',
     &'      Temp   heatgen')
 1110 format ('   (min)       x       y      (V)       (V)      (A/m2)',
     %'    (C)    (W/m2)')
 1111 format (
     &g10.3,'  rksa1,   rate constant 1 for negative side reaction'
     &/1x,g10.3,'  rksc1,   rate constant 1 for positive side reaction'
     &/1x,g10.3,'  rksa2,   rate constant 2 for negative side reaction'
     &/1x,g10.3,'  rksc2,   rate constant 2 for positive side reaction'
     &/1x,g10.3,'  rksa3,  rate constant 3 for negative side reaction'
     &/1x,g10.3,'  rksc3r,  rate constant 3 for positive side reaction')
 1114 format ('   Time     Util N  Util P  Cell Pot   Uocp      Curr',
     &'    Temp    Heat Gen')
 1115 format ('   (min)       x       y      (V)       (V)      (A/m2)',
     %'  (C)    (W/m2)')
 1116  format(/i7, '    nside, side reaction flag ')
      write (2,*) ' '

      go to (131,132,133,134,135,136,137,138),nneg
  131 write (2,*) 'Li foil'
      go to 147
  132 write (2,*) 'Carbon (petroleum coke)'
      go to 147
  133 write (2,*) 'MCMB 2528 Graphite (Bellcore)'
      go to 147
  134 write (2,*) 'TiS2'
      go to 147
  135 write (2,*) 'Tungsten oxide (LixWO3 with 0<x<0.67)'
      go to 147
  136 write (2,*) 'Lonza KS6 graphite (Bellcore)'
      go to 147
  137 write (2,*) 'Metal Hydride'
      go to 147
  138 write (2,*) 'Add your own negative electrode'
  147 go to (101,102,103,104,105,106,107,108,109,110,111,112,
     & 113,114),nprop
  101 write (2,*) 'AsF6 in methyl acetate'
      go to 200
  102 write (2,*) 'Perchlorate in PEO'
      go to 200
  103 write (2,*) 'Sodium Triflate in PEO'
      go to 200
  104 write (2,*) 'LiPF6 in PC (Sony cell simulation)'
      go to 200
  105 write (2,*) 'Perchlorate in PC (West simulation)'
      go to 200
  106 write (2,*) 'Triflate in PEO'
      go to 200
  107 write (2,*) 'LiPF6 in EC/DMC and p(VdF-HFP) (Bellcore)'
      go to 200
  108 write (2,*) 'LiPF6 in EC/DMC and p(VdF-HFP) (Bellcore) cell #2'
      go to 200
  109 write (2,*) 'Ion exchange membrane, t+=1.0'
      go to 200
  110 write (2,*) 'LiTFSI in PEMO (data from Steve Sloop)'
      go to 200
  111 write (2,*) 'LiPF6 in EC:DMC'
      go to 200
  112 write (2,*) 'LiTFSI in PEO at 85 C (data from Ludvig Edman)'
      go to 200
  113 write (2,*) 'KOH in H2O'
      go to 200
  114 write (2,*) 'add your own electrolyte'
  200 go to (201,202,203,204,205,206,207,208,209,210,211,212,213),npos
  201 write (2,*) 'TiS2'
      go to 300
  202 write (2,*) 'Spinel Mn2O4 (lower plateau)'
      go to 300
  203 write (2,*) 'NaCoO2:  Sodium Cobalt Oxide'
      go to 300
  204 write (2,*) 'Spinel Mn2O4 (upper plateau)'
      go to 300
  205 write (2,*) 'Tungsten oxide (LixWO3 with 0<x<0.67)'
      go to 300
  206 write (2,*) 'LiCoO2 (Cobalt dioxide)'
      go to 300
  207 write (2,*) 'V2O5 (Vanadium oxide)'
      go to 300
  208 write (2,*) 'LiNi0.8Co0.2O2 (Nickel dioxide Gen 1)'
      go to 300
  209 write (2,*) 'Spinel Mn2O4 (Bellcore)'
      go to 300
  210 write (2,*) 'V6O13 (Vanadium oxide)'
      go to 300
  211 write (2,*) 'LiAl0.2Mn1.8O4F0.2 spinel from Bellcore'
      go to 300
  212 write (2,*) 'NiOOH'
      go to 300
  213 write (2,*) 'Add your own'
  300 continue

      if (mvdc1.eq.1) write (2,*)
     &'Variable solid-phase diff coeff, anode'
      if (mvdc3.eq.1) write (2,*)
     &'Variable solid-phase diff coeff, cathode'

c     Convert coulombic capacity to total concentrations in ekin
      ct1=3.6d03*cot1*rs1/fc
      ct3=3.6d03*cot3*rs3/fc
c     write (2,*) 'N/P capacity ratio',(h1*ct1*csx*(1-ep1-epp1-epf1-epg1))/
c    &(h3*ct3*(1-csy)*(1-ep3-epp3-epf3-epg3))
c     write (2,*) 'Ah/m2, positive',(h3*cot3*(1-ep3-epp3-epf3-epg3)*rs3*
c    &(1-csy))
c Calculate the amount of gas in the headspace (1.4 factor for 2005 Prius battery)
c This is necessary only for a gas-phase material balance
      vol=(epg1*h1+epg2*h2+epg3*h3)/((epg1*h1+epg2*h2+epg3*h3)+
     &0.5d0*(h1+h2+h3))/1.4d0
      if (nside.ge.1.and.nprop.ne.13) vol=1.0d0 !vol=1 for liquid-phase side rxn
      lcount=1 !counter for max power

      shape1=3.0d0 ! negative particle,3. for spherical, 1. for planar, 2. for cylindrical
      shape3=3.0d0 ! positive particle,3. for spherical, 1. for planar, 2. for cylindrical
      if (imp.eq.0) then 
      cap1=capp1 ! F/m2, capacitance for negative
      cap3=capp3 ! F/m2, capacitance for positive
      else ! For impedance set cap1 and cap3 to zero.
      cap1= 0.0d0 ! F/m2, capacitance for negative
      cap3= 0.0d0 ! F/m2, capacitance for positive
      endif

c     Keeps track of tapered discharge/charge mode
      do i=1,lcurs
         taper(i) = 0
         if(mc(i).eq.-1) then
            mc(i) = 1
            taper(i) = 1
         endif
      enddo
c
c     Convert times to seconds and sum up times of mode changes
      if (mc(1).lt.2) then
      tot(1)=6.0d01*tt(1)
      else
      tot(1)=0.0d0
      end if
      do 51 i=2,lcurs
      if (mc(i).lt.2) then
      tot(i)=tot(i-1)+6.0d01*tt(i)
      else
      tot(i)=tot(i-1)
      end if
   51 continue
      
c     ** This section defines the particle size and distribution **
c     Single particle in each electrode
c     Specific area calculated from geometry
      area3=shape3*(1.0d0-ep3-epf3-epp3-epg3)/Rad3
      if (n1 .gt. 0) then
      area1=shape1*(1.0d0-ep1-epf1-epp1-epg1)/Rad1 
      else
      area1 = 1.0d0/h1
      endif

c Set up the radius and area of each particle.  The following six lines
c apply to electrodes with a single particle, or electrodes with
c multiple particles all of the same radius and area.
      if (npa.eq.1) then
      do mpa=1,npa
      Rad1pa(mpa)=Rad1 ! start with uniform particles
      area1pa(mpa)=area1/dble(npa)
      Rad3pa(mpa)=Rad3 ! start with uniform particles
      area3pa(mpa)=area3/dble(npa)
      enddo ! mpa
      vf1(1)=1.0d0
      vf3(1)=1.0d0
      else

c For nonuniform particles use this section.
c We use volume fractions to specify our distribution; the user should
c be able to calculate how to turn other distributions, such as number
c or area fraction, into volume fractions.

      vol1=(1.0d0-ep1-epf1-epp1-epg1) ! Solid vol frac in negative
      vol3=(1.0d0-ep3-epf3-epp3-epg3) ! Solid vol frac in positive

c Here input the radii of each particle.  This is set up for a two-particle
c distribution, but can be extended to as many as five particles by
c adding additional particle radii here.

      Rad1pa(1)=10.0d-6
      Rad1pa(2)=10.0d-6
      Rad1pa(3)=1.0d-6
      Rad1pa(4)=1.0d-6
      Rad1pa(5)=1.0d-6

      Rad3pa(1)=5.0d-6
      Rad3pa(2)=5.0d-6
      Rad3pa(3)=1.0d-6
      Rad3pa(4)=1.0d-6
      Rad3pa(5)=1.0d-6

c Here input the relative volume fraction of each particle.  Note that
c these need to add to 1.0, not to vol1 or vol3!  For example, a
c relative volume fraction of 0.2 means that 20% of the solid volume of an
c electrode is composed of particles of a given radius.

      vf1(1)=0.2d0
      vf1(2)=0.8d0
      vf1(3)=0.2d0
      vf1(4)=0.2d0
      vf1(5)=0.2d0

      vf3(1)=0.2d0
      vf3(2)=0.8d0
      vf3(3)=0.2d0
      vf3(4)=0.2d0
      vf3(5)=0.2d0

c Reset area1 and area3 for nonuniform particles.
      area1=0.0d0
      area3=0.0d0
c
      do mpa=1,npa
      area1pa(mpa)=vol1*vf1(mpa)/(4.d0/3.d0*pi*Rad1pa(mpa)**3.0d0)*
     &4.d0*pi*Rad1pa(mpa)**2.d0
      area3pa(mpa)=vol3*vf3(mpa)/(4.d0/3.d0*pi*Rad3pa(mpa)**3.0d0)*
     &4.d0*pi*Rad3pa(mpa)**2.d0
      area1=area1+area1pa(mpa)
      area3=area3+area3pa(mpa)
      enddo !mpa
      endif
c     Get effective solid-phase conductivity with Bruggeman
      sig3=sig3*((1.0d0-ep3-epp3-epg3)**1.5d0)
      sig1=sig1*((1.0d0-ep1-epp1-epg3)**1.5d0)
c
      h2=h2/dble(n2-1)
      h3=h3/dble(n3)
      if (n1 .gt. 1) h1=h1/dble(n1)
      h=h2
      hha=rad1/dble(n4-1)
      hhc=rad3/dble(n4-1)
      frt=fc/(r*t)
c
      mpa=1 ! only one size for initialization
c     Find initial solid-phase potential guesses
c     from initial solid concentrations:
      call ekin(1,1,1,csx,mpa,mvdc1,mvdc3)
      write (2,*) 'open-circuit potential, negative ',g0
      write (2,*) 'delta S negative ', dudt
      Ua=g0
      call ekin(nj,1,1,csy,mpa,mvdc1,mvdc3)
      write (2,*) 'open-circuit potential, positive ',g0
      write (2,*) 'delta S positive ', dudt
      Uc=g0
      OCP=Uc-Ua ! open-circuit cell potential

      write (2,*) ' '
      write (2,*) '     DUAL INSERTION CELL VERSION 5.1'
      write (2,*) ' '
      if (imp.eq.1) then
      write (2,*) '         IMPEDANCE OPTION   '
      write (2,*) ' '
      write (2,*) '  Omega                Z real         Z imag  '
      write (2,*) ' (Rad/s)              (ohm cm2)      (ohm cm2)  '
      go to 54
      endif ! end of impedance section (output)

        if (nside.ne.0.and.nneg.eq.7) then
        write (2,1107) !This output includes pressures
        write (2,1108)
        else
        write (2,1109) !This output excludes pressures
        write (2,1110)
        endif
      write (2,*) ' '

      write (6,*) 'Resistances in each portion of the cell, from a',
     1' simple (delta V) / I'
      write (6,*) 'Time       Cur Den      R Anode       R Separator'
     &,'      R Cathode        R Total'
      write (6,*) '(min)      (A/m^2)     (ohm-m^2)       (ohm-m^2)'
     &,'       (ohm-m^2)      (ohm-m^2)'
      write (6,*) ' '
   54 continue
c
c***********************************************************************
c
c     initialize time counting variables
      k=1
      time=0.0d0
      time2=0.0d0
      rr=0.0d0
      coul=0.d0
      L=1
      ts(1)=0.0d0

      do j=1,nj
        xx(1,j)=xx(1,n1+2) ! initial concentration (mol/m3)
c     Uniform initial concentration if lflag=1.
c     Step-function initial concentration if lflag=0.
        if(lflag.eq.0 .and. (j.le.n1+1 .or. j.ge.n1+n2))
     &  xx(1,j)=1.0d-01
      enddo

c     assume current density linear in electrodes
      call prop(nj,n2,n1)
      curold=0.d0 ! open circuit
      vvold=OCP
      cur=1.d0 ! guess current to get internal resistance
c      cur=1000.d0 ! If high currents are expected, use a high current
c                  here to get an appropriate resistance
      if(mc(1).ge.1.or.mc(1).eq.-1) cur=cu(1)
      if (cap1.ne.0.0d0.or.cap3.ne.0.0d0) cur=0.0d0

      do j=1,n1+1 ! negative electrode
        do mpa=1,npa
          xx(2+mpa,j)=ct1*csx
          utz(mpa,j)=csx
        enddo ! mpa
      enddo
      do j=n1+2,n1+n2-1 ! separator
        do mpa=1,npa
          xx(2+mpa,j)=0.d0
          utz(mpa,j)=0.d0
        enddo ! mpa
      enddo
      do j=n1+n2,nj ! positive electrode
        do mpa=1,npa
          xx(2+mpa,j)=ct3*csy
          utz(mpa,j)=csy
        enddo ! mpa
      enddo

c     n1+1 is the last node in the negative
c     n1+n2 is the first node in the positive
c     etaneg=xx(kj,j)*fc*ranode
c     etapos=xx(kj,j)*fc*rcathde

      xx(kp2,nj)=0.d0
      do 78 j=nj,n2+n1,-1 ! positive
        xx(ki2,j)=cur*dble(nj-j)/dble(n3)
        xx(kj,j)=-cur/fc/h3/area3/dble(n3)
        if(j.ne.nj)
     &  xx(kp2,j)=xx(kp2,j+1)+h3*xx(ki2,j)/cd(j) ! ohmic drop is not accurate
   78 xx(kp1,j)=Uc+xx(kp2,j)+xx(kj,j)*fc*rcathde

      do j=n1+n2-1,n1+2,-1 ! separator
        xx(ki2,j)=cur
        xx(kj,j)=0.d0
        xx(kp2,j)=xx(kp2,j+1)+h2*xx(ki2,j)/cd(j) ! ohmic drop is not accurate
        xx(kp1,j)=0.0d0
      enddo

      do 79 j=n1+1,1,-1 ! negative electrode
          if (n1.gt.0) then
          xx(ki2,j)=cur*dble(j-1)/dble(n1)
          xx(kj,j)=cur/fc/h1/area1/dble(n1)
          xx(kp2,j)=xx(kp2,j+1)+h1*xx(ki2,j)/cd(j) ! ohmic drop is not accurate
          xx(kp1,j)=xx(kp1,j)+xx(kp2,j)+xx(kj,j)*fc*ranode
          else
          xx(ki2,j)=cur*dble(j-1)/dble(n1+1)
          xx(kj,j)=cur/fc/h1/area1/dble(n1+1)
          xx(kp2,j)=xx(kp2,j+1)+h1*xx(ki2,j)/cd(j) ! ohmic drop is not accurate
          xx(kp1,j)=Ua+xx(kp2,j)+xx(kj,j)*fc*ranode
          endif
   79 continue
c
c initialization for variables kS1,kj1,kS2,kj3 ! side reactions, not impedance
      do j=1,nj
      if (nside.ge.1) then
      xx(kS1,j)=c1init
      xx(kS2,j)=c2init
      endif
      xx(kj1,j)=0.0d0
      xx(kj2,j)=0.0d0
      xx(kj3,j)=0.0d0
      enddo

      do j=1,nj
      do i=1,n
      xt(i,j,1)=xx(i,j)
      enddo
      enddo

c Assumed temperature dependence for solid-state diffusion
      dfs1=dfs1save*dexp((EbarS1)*(t-298.15d0)/(t*298.15d0))
      dfs3=dfs3save*dexp((EbarS3)*(t-298.15d0)/(t*298.15d0))

c Initialize solid-phase concentrations. PSA
      if (mvdc1.eq.1.or.mvdc3.eq.1) then
      mpa=1 ! only one size for variable solid-phase diffusion coefficient
        do jj=1,nj
          do jjj=1,nnj
            css(jj,jjj)=xx(2+mpa,jj)
            cssold(jj,jjj)=xx(2+mpa,jj)
c For a true variable solid-phase diffusion coefficient put in the
c functional dependence here and in subroutine vardc
            if (jj.le.n1+1) ds(jj,jjj)=dfs1
            if (jj.ge.n1+n2) ds(jj,jjj)=dfs3
          enddo
        enddo
      endif

      mcL=1
      if (imp.eq.0) call comp(n,lim,k,rr,0,nflag,0,jcount)
	if(rr.eq.0.d0) print *, jcount, ' iterations at zero time step'
      if(nconv.eq.1) write (2,*) 'comp 1 not converged' ! initial
         if(nconv.eq.1) go to 533
      vv=xx(kp1,nj)-xx(kp1,1)
      Rint=-(vv-vvold)/(cur-curold)+RG ! internal resistance
      OCP=vv+(Rint-RG)*cur
      pwrmax=OCP**2/Rint/4.d0
      if (pwrmax.le.cu(L) .and. L.eq.-2) then
      write (2,*) 'WARNING: Power is above estimated power maximum'
      write (2,*) 'Requested power= ',cu(L)
      write (2,*) 'Power maximum= ',pwrmax
      write (2,*) 'Lower the power or change the initial current guess.'
      endif

      cuL=cu(L)
      cur=cu(L) ! good for mc(L)=1, 2, and -1.
c     First guess for current for constant load.
      if(mc(L).eq.-3) cur=OCP/(cu(L)+Rint)
      if(mc(L).eq.0) cur=(OCP-cu(L))/Rint ! constant potential
      if(mc(L).eq.-2) then ! constant power
        cur=cu(L)/OCP
        do i=1,6 ! refine estimate of current for given power
          cur=cur-(cur*OCP-Rint*cur**2-cu(L))/(OCP-2.d0*Rint*cur)
        enddo
      endif
      
	avgeta=0.d0
      if(mc(L).eq.0) avgeta=(OCP-cuL)/2.d0
      if(mc(L).le.-2) avgeta=cur*(Rint-RG)/2.d0

c     Must activate lpow=1 in data file if you want peak powers:
      kkflag=0
      nflag=0

c     calculate mass (kg/m2) of the cell
      call mass(re,rs3,rs1,rf,rpl,rc,rcn,rcp)
      dens=tw/thk

c     Impedance section:
      if (imp.eq.1) then
c     initialize steady state to make impedance runs.
      imp=0
      call comp(n,lim,k,rr,0,nflag,0,jcount)

      imp=1
      ji=1
c     omega is the array of perturbation frequencies in Rad/s
      do 61 iii=1, 110

      omega(ji)=10.d6/(10.d0**0.1d0)**iii
c
      call comp(n,lim,k,rr,0,nflag,0,jcount)
c
   61 continue
      write (2,*)''
      write (2,*)'  Omega            Z real        Z imag  '
      write (2,*)' (Rad/s)          (ohmcm2)      (ohmcm2) '
      do 62 i=1,110
   62 write (2,57)omega(i),vreal(i)*10000,vimag(i)*10000
      write (2,*)''
      write (2,*)'  Omega           Z reala       Z imaga  '
      write (2,*)' (Rad/s)          (ohmcm2)      (ohmcm2) '
      do 63 i=1,110
   63 write (2,57)omega(i),vreala(i)*10000,vimaga(i)*10000
      write (2,*)''
      write (2,*)'  Omega           Z realc       Z imagc  '
      write (2,*)' (Rad/s)          (ohmcm2)      (ohmcm2) '
      do 64 i=1,110
   64 write (2,57)omega(i),vrealc(i)*10000,vimagc(i)*10000
      write (2,*)''
      write (2,*)'  Omega             Z mag         phi  '
      write (2,*)' (Rad/s)           (ohmcm2)      (deg) '
      do 65 i=1,110
   65 write (2,57)omega(i),zmag(i)*10000,phi(i)
      write (2,*)''
      write (2,*)'  Omega            Z maga        phia  '
      write (2,*)' (Rad/s)           (ohmcm2)      (deg) '
      do 66 i=1,110
   66 write (2,57)omega(i),zmaga(i)*10000,phia(i)
      write (2,*)''
      write (2,*)'  Omega             Z magc       phic  '
      write (2,*)' (Rad/s)           (ohmcm2)      (deg) '
      do 67 i=1,110
   67 write (2,57)omega(i),zmagc(i)*10000,phic(i)
      write (2,*)''
      go to 999
      endif ! end of impedance section (general)

c     Obtain better initial guess by running at constant current
c     which is close to the current for the chosen constraint.
      mcL=1
        call comp(n,lim,k,rr,0,nflag,0,jcount)
        if(nconv.eq.1) write (2,*) 'comp 2 not converged' ! refine
c     write (2,*) 'cellpot 1'
        call cellpot(k,vv,0,0,lflag)

      mcL=mc(L)
      cuL=cu(L)
      curold=0.d0
      vvold=OCP

      if(mc(L).eq.-2) then ! constant power
      LS=2 ! Bypass this for speed.
      if(LS.ne.0) go to 56
      mcL=mc(L)
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 4 not converged' ! power
      curold=xx(ki2,n1+1)
      vvold=xx(kp1,nj)-xx(kp1,1)
      cur=curold-1.d0
      mcL=1
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 5 not converged' ! get Rint
      vv=xx(kp1,nj)-xx(kp1,1)
      Rint=-(vv-vvold)/(cur-curold)+RG ! internal resistance
      OCP=vv+(Rint-RG)*cur
      pwrmax=OCP**2/Rint/4.d0
   56 continue

      mcL=mc(L)
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 6 not converged' ! repeat?
      call cellpot(k,vv,1,0,lflag)
c     write (2,*) 'cellpot 3'
      endif ! treatment of constant power
c
c     rr is the size of a time step.
      if(cap1.eq.0.d0 .and. cap3.eq.0.d0) then
      rr=0.02d0
c     tmmax !for constant time steps
      else
      rr=1.5d-13 !  initial time step is 1.5d-13 second
      endif

      ed=0.d0
      heat=0.d0
      qlosstot=0.d0
      iflag=0
      L=0
   53 L=L+1 ! new leg
      mcL=mc(L)
      cuL=cu(L)
      if(L.gt.1) rr=0.d0
  123 k=k+1 ! new time step
      nt=k-1

c     adjust time step to match time of change in current
      time=ts(k-1)+rr
      if(time .ge. tot(L) .and. mc(L).lt.2) then
      rr=tot(L)-ts(k-1)
      time=tot(L)
      iflag=1
      endif
c
  129 ts(k)=ts(k-1)+rr ! new time step
      call calca(k) ! should need only one call at each time step
      dtnow=rr
c
      if(mc(L).eq.-3 .or. mc(L).ge.0) then
c     mc is 1 or 2 so run galvanostatically
c     or mc is 0; run at constant potential
c     or mc is -3; run at constant load
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 7 not converged' ! cur,pot,load
        if (nconv.eq.1) go to 533
c
      else ! iterate on power, mc is -2
c     The strategy is to use the constant-power mode in comp if
c     the possible maximum power is greater than that called for.
c     Otherwise, one wants to report the maximum power.  (This
c     is a little like using a taper at constant potential when
c     a constant-current run reaches a cutoff potential.)
c     The maximum-power mode is abandoned for a while.  2-16-08
      do 610 j=1,nj
      term_s1s(j)=term_s1(j)
  610 terms(j)=term(j)

      LS=3 ! activate this section to get estimate of pwrmax
c     at each time step.
c     This does not get a true pwrmax, since it does not iterate
c     and take the current to the true curmax.  It gives an
c     approximation by calculating Rint from a small change in the
c     current.
      if(LS.ne.0) go to 58
      mcL=1 ! comp is run in constant-current mode
      call comp(n,lim,k,rr,0,nflag,0,jcount) ! at previous current
      if(nconv.eq.1) write (2,*) 'comp 8 not converged' ! power
         if (nconv.eq.1) go to 533
      curold=cur
      vvold=xx(kp1,nj)-xx(kp1,1)
      pwr=cur*vvold
      cur=cur-1.d0 ! decrement current to get internal resistance
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 9 not converged' ! to get Rint
         if (nconv.eq.1) go to 533
      vv=xx(kp1,nj)-xx(kp1,1)
      Rint=-(vv-vvold)/(cur-curold)+RG ! internal resistance
      OCP=vv+(Rint-RG)*cur
      pwrmax=OCP**2/Rint/4.d0
      curmax=OCP/Rint/2.d0 ! estimate of current at max power
   58 continue
      if(pwrmax.gt.cu(L)) then
      mcL=mc(L) ! specified power should be possible
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 10 not converged' ! power
      else

c     iterate to find maximum power ! really don't want to go here
      mcL=-4 ! try to let comp do it
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 11 not converged' ! max power
         if (nconv.eq.1) go to 533
      vv=xx(kp1,nj)-xx(kp1,1)
      endif ! end of treatment of constant power

      endif ! end time-step calculation for set mode

      if(rr.lt.0.9999*dtnow) iflag=0 ! time step decreased in comp
      mcL=mc(L)
      if(mc(L).eq.2) then
        call cellpot(k,vv,0,0,lflag)
      else
        call cellpot(k,vv,1,0,lflag)
      endif
      if(rr.ne.0.d0) coul=coul+rr*xx(ki2,n1+1)
      if (il1.eq.1.and.mod(k,il3).eq.0) call nucamb(k,il2)
      frt=fc/r/t

      sign=1.d0
      vcut = vcutlo(L)
      IF(VV.LT.VCUT.and.cur.eq.0.d0) then
        do i=L+1,Lcurs
          tot(i) = tot(i) + time - tot(L) !Correct time markers
        enddo
        go to 888
      endif
      if(cu(L).lt.0.d0) then
         sign=-1.d0
         vcut = vcuthi(L)
      endif
      IF(sign*VV.LT.sign*VCUT .and. cu(L).ne.0.d0) then
         if(taper(L).eq.1) then !Taper the current by holding VV at VCUT
      if(mc(L).eq.1) then
            mc(L) = 0
            mcL = 0
            cu(L) = vcut
            cuL=vcut
            k=k+1
            ts(k)=ts(k-1)
      call comp(n,lim,k,0.d0,0,nflag,0,jcount)
      call cellpot(k,vv,1,0,lflag)
      if(il1.eq.1) call nucamb(k,il2)
      endif
         else                   !End profile segment
            do i=L+1,Lcurs
               tot(i) = tot(i) + time - tot(L) !Correct time markers
            enddo
            go to 888
         endif
      endif
c     check to see if cutoff potential is exceeded if mc is 2
      if (mc(L).eq.2) then
        IF ((VV.LT.TT(L) .AND. CU(L).GT.0.0) .OR.
     &  (VV.GT.TT(L) .AND. CU(L).LT.0.0)) THEN
          if (dabs(vv-tt(l)) .gt. 1.0d-04) then
            rr=rr/2.0d0
            iflag=1
            go to 129
          else
            time2=time2+rr
c     write (2,*) 'cellpot 5'
            call cellpot(k,vv,1,0,lflag)
            frt=fc/r/t
            iflag=1
          endif
        else
          iflag=0
          time2=time2+rr
c     write (2,*) 'cellpot 6'
          call cellpot(k,vv,1,0,lflag)
          frt=fc/r/t
        end if
      end if
c
c     Increasing time steps: comment out for constant time steps
      rrmax=tmmax
      if(k.le.5) rrmax=0.5d0
      jlim=int(lim/5.0d0)
      if (jlim.lt.5) jlim=5
      if(jcount.lt.jlim.and. k.gt.2 .and. rr.lt.rrmax .and.
     1iflag.eq.0.and.lcount.le.1) then
      if(xx(ki2,n1+5).lt.1200.d0) rr=rr*1.1d0

c      print *,'next time step increased to ', rr,' (s)'
      end if

      if(k.GE.maxt-1) then
        write (2,*) 'kmax=',k,' a larger matrix needed for xt'
      endif
      if (k.GE.501) then ! trim stored solid concentrations
      do 92 kk=3,401,2
      kput=(KK+1)/2
      ts(kput)=ts(kk)
      do 92 j=1,nj
      do 92 i=1,n
   92 xt(i,j,kput)=xt(i,j,kk)
      do 93 kk=402,K
      ts(kk-200)=ts(kk)
      do 93 j=1,nj
      do 93 i=1,n
   93 xt(i,j,kk-200)=xt(i,j,kk)
      k=k-200
      endif

      if (iflag.eq.0 .and. rr.eq.0.d0) then
c     rr is the size of a time step.
      if(cap1.eq.0.d0 .and. cap3.eq.0.d0) then
      rr=0.02d0
c      tmmax !for constant time steps
      else
      rr=1.5d-13 !  initial time step is 1.5d-13 second
      endif
      endif
      if (iflag .eq. 0) go to 123 ! a new time step
  888 continue
      iflag=0
      if (mc(L).eq.2) then
      do 124 m=L,Lcurs
  124 tot(m)=tot(m)+time2
      time2=0.0d0
      end if
      IF(L.EQ.LCURS .AND. LCURS.GE.100) THEN 
c If Lcurs greater than 100 go into an endless loop
      L=0
      tot(1)=TOT(LCURS)
      if (mc(1).lt.2) TOT(1)=TOT(1)+60.0D0*TT(1)
      do 403 i=2,Lcurs
      if (mc(i).lt.2) then
      tot(i)=tot(i-1)+6.0d01*tt(i)
      else
      tot(i)=tot(i-1)
      end if
  403 continue
      ENDIF

  533 if(L.lt.Lcurs) then ! prepare for a new leg
      rr=0.d0
c       calculate zero-time solution for change in current
c     Need some initialization here to get started smoothly on
c     a new leg of the calculation.
      cuL=cu(L+1)
c     First guess for current
      cur=cu(L+1) ! constant current
      if(mc(L+1).eq.-3) cur=OCP/(cu(L+1)+Rint) ! constant load
      if(mc(L+1).eq.0) cur=(OCP-cu(L+1))/Rint ! constant potential

      if(mc(L+1).eq.-2) then ! constant power
      LS=1
      if(LS.eq.1) go to 55 ! bypass maximum-power calculations
c     first calculate power maximum
      curmax=OCP/Rint/2.d0
      pwrmax=OCP**2/Rint/4.d0

      do ii=1,4 ! could use fewer iterations here
      cur=curmax
      mcL=1
      call comp(n,lim,k,rr,0,nflag,0,jcount) 
      if(nconv.eq.1) write (2,*) 'comp 12 not converged'
      curold=cur ! save correctly calculated values at one cur
      vvold=xx(kp1,nj)-xx(kp1,1)
      cur=cur-1.d0 ! decrement current
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 13 not converged'
      vv=xx(kp1,nj)-xx(kp1,1)
      Rint=-(vv-vvold)/(cur-curold)+RG
      OCP=vv+(Rint-RG)*cur
      curmax=OCP/Rint/2.d0
      pwrmax=OCP**2/Rint/4.d0
      enddo ! ii
      if(pwrmax.lt.cuL) then ! can't deliver set power
      cur=curmax
      call comp(n,lim,k,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 14 not converged'
      mcL=mc(L)
      go to 53
      endif

   55 continue
      cur=cu(L+1)/OCP
      do ii=1,6 ! refine estimate of current for given power
c     could solve quadratic
      cur=cur-(cur*OCP-Rint*cur**2-cu(L+1))/(OCP-2.d0*Rint*cur)
      enddo
      endif ! end of power loop

      avgeta=cur*Rint/2.d0
      if(mc(L+1).eq.0) avgeta=(OCP-cuL)/2.d0

      Ua=xx(kp1,1)
      Uc=xx(kp1,nj)
      Ua=0.d0 ! need better open-circuit potentials
      Uc=4.d0 ! need better open-circuit potentials
      do 73 j=1,(n1+1) ! negative electrode
      do mpa=1,npa
      xt(2+mpa,j,k+0)=xx(2+mpa,j)
      enddo ! mpa
      xt(ki2,j,k+0)=cur*dble(j-1)/dble(n1)
      xt(kj,j,k+0)=cur/fc/h1/area1/dble(n1)
   73 xt(kp1,j,k+0)=xx(kp2,j)+Ua+0.60d0*avgeta

      do j=n1+1,n1+n2 ! separator
      xt(ki2,j,k+0)=cur
      enddo
c
      do 75 j=(n2+n1),nj ! positive electrode
      do mpa=1,npa
      xt(2+mpa,j,k+0)=xx(2+mpa,j)
      enddo ! mpa
      xt(ki2,j,k+0)=cur*dble(nj-j)/dble(n3)
      xt(kj,j,k+0)=-cur/fc/h3/area3/dble(n3)
   75 xt(kp1,j,k+0)=xx(kp2,j)+Uc-1.20d0*avgeta

c     straighten out the initialization
      mcL=1
      call comp(n,lim,k+1,rr,0,nflag,0,jcount)
      if(nconv.eq.1) write (2,*) 'comp 15 not converged'
      mcL=mc(L+1)

        if(mc(L+1) .gt. 0) then
      if (nconv.eq.0) k=k+1
          ts(k)=ts(k-1)
          rr=0.0d0 !run a zero time step at the start
c     mcL=1
          call comp(n,lim,k,rr,0,nflag,0,jcount) 
      if(nconv.eq.1) write (2,*) 'comp 16 not converged'
c     write (2,*) 'cellpot 7'
          call cellpot(k,vv,1,0,lflag)
      if (il1.eq.1.and.mod(k,il3).eq.0) call nucamb(k,il2)
        endif
      endif
c      rr = tmmax !uncomment and comment out next 5 lines for constant time steps
      if(cap1.eq.0.d0 .and. cap3.eq.0.d0) then
        rr=0.02d0   !  initial time step is 0.2 seconds
      else
        rr=1.5d-13 !  initial time step is 1.5d-13 second
      endif
      nconv=0
      if(L.le.Lcurs) then
      if (L.eq.1) pow=ed/tw/ts(k)
      if (L.ne.1) pow=ed/tw/(time-timesave)
      timesave=ts(k)
c     write (2,*) ' '
      if (L.eq.Lcurs) write (2,44) tw
      write (2,48) ed/tw/3.6d03
      write (2,49) pow
      write (2,47) heat
      write (2,*) ' '
      ed=0.0d0
      endif
      if(jsol.gt.0 .and. jsol.lt.nj) call sol(k,jsol)
c     Test for seoond leg. WHT 3/10/14
c      RG=RGin(2)
c      RGn=RG/3.d0
c      RGp=RG-RGn
c      RGext=RG/4.d0
      IF(L.LT.LCURS) GO TO 53 ! new leg

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Additional features section:
c
c     peak-power subroutine:
      if(lpow.eq.1) then
      il1=0
      call peak(n,lim,cu(1), vcut)
      endif
c
c     Solid-phase concentration profiles at given time and position
      if(jsol.gt.0 .and. jsol.lt.nj) call sol(k,jsol)

      stop
  999 end
c
c***********************************************************************
      subroutine comp(n,lim,kk,tau,kkflag,nflag,lpow,jcount)
c     Main calculational subroutine, going through the
c     equations that are solved through the cell sandwich.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/power/ ed,Vold,ranodesave,rcathdesave,heat,qlosstot
      common/ssblock/ xp0(16),xx0(16,221),term(221)
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/imp/ vreal(110),vimag(110),phi(110),zmag(110),omega(510),
     1vreala(110),vrealc(110),vimaga(110),vimagc(110),
     1phia(110),phic(110),zmaga(110),zmagc(110)
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/tprop/df(221),cd(221),tm(221),ddo2(221),ddh2(221),
     1ddf(221),dcd(221),dtm(221),dfu(221),d2fu(221),do2(221),dh2(221)
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
      common/side/rksc1,c1init,c2init,rksa1,term_s1(221),vol,
     &rksa2,rksa3,rksc3,rksc2,UsO2,UsH2,cn2,term_s2(221),nside
      common/mat/ b,d
      common/bnd/ a,c,g,x,y
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      common /maxpow/ lcount
      common/RG/ RG,RGn,RGp,RGext
      save iijj
      dimension b(17,17),d(17,35),termn(221)
      dimension termn_s1(221), termn_s2(221)
      dimension a(17,17),c(17,221),g(17),x(17,17),y(17,17)
   99 format (1h ,//5x,'this run just did not converge'//)
  103 format(f13.5,'      ',f12.6,'     ',f12.6)
  104 format (i4,' ',g15.5,' ',g15.5,' ',
     &g15.5,' ',g15.5,' ',g15.5,' ',g15.5)
      ntot=n
      if(imp.eq.1) ntot=2*n
      nx=ntot
      mcLL=mcL
      if(mcLL.eq.-4) mcL=1 ! use const cur for maximum power
      Lcount=0 ! counter for maximum power
      curold=cur ! save correctly calculated values at one cur
      vvold=xx(kp1,nj)-xx(kp1,1)
      rrold=rr
      rr=tau
      if(mcLL.eq.-4) cur=cur-1.d0

      omi=dble(1-imp)
      exbrug=1.5d0
c     The exbrug exponent is being set independently for the separator and
c     positive and negative electrodes as shown in the following three lines...

      exbrug1=1.5d0 !EX for the negative active material.
      exbrug2=1.5d0 !EX for the separator material.
      exbrug3=1.5d0 !EX for the positive active material.	

c     The following statement is for activating the shut-down spearator feature.
      Shutdown=0.0 !If shutdown=0. shut-down separator not active,  
c     if shutdown=1.0   Calcualation performed in comp subroutine. active

c     The following statement may need to be revised due to the above changes.	 
      if(nprop.eq.7 .or. nprop.eq.8) exbrug=3.3d0
c
  666 kadd=0
c     sets first guess to last time-step values
      if(rr.eq.0.d0 .and. lpow.eq.1) kadd=1
      if(kk .eq. 1) kadd = 1
      do 1 j=1,nj
      do 1 i=1,ntot
      if(rr.ne.0.d0) xx(i,j)=xt(i,j,kk-1+kadd)
    1 c(i,j)=xx(i,j)

      do 2 j=1,nj
      do 2 i=1,nnj
    2 css(j,i)=cssold(j,i)
c
c     initialize variables to begin each iteration (jcount is iteration #)
    3 Lcount=Lcount+1
      jcount=0
      do 4 i=1,n
    4 xp(i)=0.0d0

      if(time.eq.0.d0) then
      totLiold=0.d0
      do j=1,nj
      fh=h3*ep3
      if(j.le.n1+n2) fh=h2*ep2
      if(j.le.n1+1) fh=h1*ep1
      if(j.eq.1) fh=0.5*h1*ep1
      if(j.eq.n1+1) fh=0.5d0*(h1*ep1+h2*ep2)
      if(j.eq.n1+n2) fh=0.5d0*(h3*ep3+h2*ep2)
      if(j.eq.nj) fh=0.5d0*h3*ep3
      totLiold=totLiold + fh*c(1,j)
      enddo
      print *, totLiold
      write (3,*) 'totLiold ', totLiold
	endif
      totLio=0.d0
      do j=1,nj
      fh=h3*ep3
      if(j.le.n1+n2) fh=h2*ep2
      if(j.le.n1+1) fh=h1*ep1
      if(j.eq.1) fh=0.5*h1*ep1
      if(j.eq.n1+1) fh=0.5d0*(h1*ep1+h2*ep2)
      if(j.eq.n1+n2) fh=0.5d0*(h3*ep3+h2*ep2)
      if(j.eq.nj) fh=0.5d0*h3*ep3
      totLio=totLio + fh*c(1,j)
      enddo
	print *, totLiold,totLio
      
c
    8 j=0
      jcount=jcount+1
      jcountsum=jcountsum+1
      if(mcL.le.-2 .or. mcL.eq.0) cur=xx(ki2,n1+1)
c     calculate physical properties
      call prop(nj,n2,n1)
c
c     initialize x and y for iteration
      do 9 i=1,ntot
      do 9 k=1,ntot
      x(i,k)=0.0d0
    9 y(i,k)=0.0d0
c
c     store previous iteration of (xp in xp0)  &  (xx in xx0)
      do 6 i=1,n
      xp0(i)=xp(i)
      do 6 jj=1,nj
    6 xx0(i,jj) = xx(i,jj)  

c     for a given iteration, set up governing equations and bc's
c     start at the left interface and move across cell.
c     initialize a,b,d,g arrays for each node
c
   10 j=j+1
      hC=h2
      if(j.le.n1+1)  hC=h1
      if(j.le.n1+1)  area=area1
      if(j.le.n1+1)  capp=capp1
      if(j.le.n1+1)  cap=cap1
      if(j.ge.n1+n2) hC=h3
      if(j.ge.n1+n2) area=area3
      if(j.ge.n1+n2) capp=capp3
      if(j.ge.n1+n2) cap=cap3
      if(j.eq.1 .or. j.eq.n1+1 .or. j.eq.n1+n2 .or. j.eq.nj)
     &area=area/2.d0

c
      do 11 i=1,ntot+1
      g(i)=0.0d0
      xx(i,j)=c(i,j)
      do 11 k=1,ntot+1
      a(i,k)=0.0d0
      b(i,k)=0.0d0
   11 d(i,k)=0.0d0

c
c     Added if-then construct to account for impedance.  While
c     impedance will be taken at zero time, we must solve for
c     concentration profiles. (JPM - 29 Jun 95)
      if(imp.eq.1) then
      rr=1.0d0 ! set rr to 1.0 to avoid divide by 0
      else ! end of impedance section (preliminary)
      if(rr.le.0.0d0) then ! treat as a zero time step
      b(1,1)=1.0d0
      g(1)=xt(1,j,kk-1+kadd)-c(1,j) ! fix electrolyte concentration
      if (nside.ge.1) then
      b(imb1,kS1)=1.0d0
      g(imb1)=xt(kS1,j,kk-1+kadd)-c(kS1,j) !fix oxygen concentration
      if (nside.ge.2) then
      b(imb2,kS2)=1.0d0
      g(imb2)=xt(kS2,j,kk-1+kadd)-c(kS2,j) !fix hydrogen concentration
      endif
      endif
      go to 199
      endif
      endif
c     ________________________________________________
c     Equation 1, mass balance.
c
c     Adsorption behavior is fixed w/ dqn parameter
c     dqn is the dq-/dq derivative equal to either -1 (anion
c     adsorption) or 0 (cation adsorption)
c 111 dq1=0.0d0
c     dq3=0.0d0
  111 dq1=-1.0d0
      dq3=-1.0d0 ! Johnson and Newman, J. Electrochem. Soc., 118 (1971), 510-517.

      termn(j)=0.d0

c Set the stoichiometric coefficient according to the nature of the reacting ion.
c For a cation, set so=1.0d0, for an anion, set so=-1.0d0.

      so=1.0d0
      if (nprop.eq.13) so=-1.d0 ! for an anion in the intercalation reaction

      fac=1.d0
      if(j.eq.n1+2 .and. n1 .gt. 0)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2+1) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
      epn=ep2+epp2
      hn=h2
      acn=0.0d0
      dqn=0.0d0
      if(n1.gt.0 .and. j.le.n1+1) then
        epn=ep1+epp1
        hn=h1
        acn=area1*capp1
        dqn=dq1
      elseif(j.gt.n1+n2) then
        epn=ep3+epp3
        hn=h3
        acn=area3*capp3
        dqn=dq3
      endif

      if(j.gt.1) then ! deal with box to left of point.
        termn(j)=-(df(j)+fac*df(j-1))*(c(1,j)-c(1,j-1))/hn/2.d0
     &  -so*(1.d0-0.5d0*(tm(j)+tm(j-1)))*c(ki2,j-1)/fc
        a(1,1)=epn*hn*0.125d0/rr*omi
     &  -(df(j)+fac*df(j-1))/hn/4.d0+fac*ddf(j-1)*
     &  (c(1,j)-c(1,j-1))/hn/4.d0-so*dtm(j-1)*c(ki2,j-1)/4.d0/fc
        b(1,1)=epn*hn*0.375d0/rr*omi
     &  +(df(j)+fac*df(j-1))/hn/4.d0+ddf(j)*(c(1,j)-c(1,j-1))/hn/4.d0
     &  -so*dtm(j)*c(ki2,j-1)/4.d0/fc
        a(1,ki2)=so*(1.d0-0.5d0*(tm(j)+tm(j-1)))/2.d0/fc
        g(1)=-epn*hn*(0.375d0*(c(1,j)-xt(1,j,kk-1+kadd))
     &  +0.125d0*(c(1,j-1)-xt(1,j-1,kk-1+kadd)))/rr*omi
        if (imp.eq.1) then
          b(1,1+n)=-epn*hn*0.375d0*omega(ji)  /2.0d0
          a(1,1+n)=-epn*hn*0.125d0*omega(ji)  /2.0d0
          b(1+n,1)=epn*hn*0.375d0*omega(ji)   /2.0d0
          a(1+n,1)=epn*hn*0.125d0*omega(ji)   /2.0d0
        endif ! end of impedance section (Equation 1)
      endif

      fac=1.d0
      if(n1.gt.0 .and. j.eq.n1+1)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
      epn=ep2+epp2
      hn=h2
      acn=0.0d0
      dqn=0.0d0
      if(n1.gt.0 .and. j.lt.n1+1) then
        epn=ep1+epp1
        hn=h1
        acn=area1*capp1
        dqn=dq1
      elseif(j.ge.n1+n2) then
        epn=ep3+epp3
        hn=h3
        acn=area3*capp3
        dqn=dq3
      endif
      if(j.ne.nj) then ! deal with box to right of point.
c     At the foil anode, only the box to the right should be used
        termn(j)=termn(j)-(fac*df(j)+df(j+1))*(c(1,j)-c(1,j+1))/hn/2.d0
     &  +so*(1.d0-0.5d0*(tm(j)+tm(j+1)))*c(ki2,j)/fc
        d(1,1)=epn*hn*0.125d0/rr*omi
     &  -(fac*df(j)+df(j+1))/hn/4.d0+ddf(j+1)*(c(1,j)-c(1,j+1))/hn/4.d0
     &  +so*dtm(j+1)*c(ki2,j)/4.d0/fc
        b(1,1)=b(1,1)+epn*hn*0.375d0/rr*omi
     &+(fac*df(j)+df(j+1))/hn/4.d0+fac*ddf(j)*(c(1,j)-c(1,j+1))/hn/4.d0
     &  +so*dtm(j)*c(ki2,j)/4.d0/fc
        b(1,ki2)=b(1,ki2)-so*(1.d0-0.5d0*(tm(j)+tm(j+1)))/2.d0/fc
        g(1)=g(1)-epn*hn*(0.375d0*(c(1,j)-xt(1,j,kk-1+kadd))
     &  +0.125d0*(c(1,j+1)-xt(1,j+1,kk-1+kadd)))/rr*omi
        if (imp.eq.1) then
          b(1,1+n)=b(1,1+n)-epn*hn*0.375d0*omega(ji) /2.0d0
          d(1,1+n)=-epn*hn*0.125d0*omega(ji)       /2.0d0
          b(1+n,1)=b(1+n,1)+epn*hn*0.375d0*omega(ji) /2.0d0
          d(1+n,1)=epn*hn*0.125d0*omega(ji)        /2.0d0
        endif ! end of impedance section (equation 1)
      endif
      g(1)=g(1)+(termn(j)+term(j))/2.d0

      if (imp.eq.1) g(1)=0.0d0

c     ________________________________________________
c Equation ii2, material balance in solid insertion material.
c This material balance includes the option of having side-reaction
c species 2 (hydrogen in the NiMH case) absorb or desorb from
c the solid phase in addition to the electrochemical intercalation
c or deintercalation of hydrogen. PSA 4/22/06

      if(imp.eq.0) then ! impedance is off

        if(j.le.n1+1) then ! in the anode
          if (n1 .eq. 0) then  ! foil anode
            mpa=1 ! only one size with foil anode
            g(ii2) = -c(2+mpa,j)
            b(ii2,2+mpa) = 1.0d0
          else ! nonfoil anode
            if (mvdc1.eq.0) then ! superposition
              if (nside.eq.3) then ! side reaction with hydrogen
                mpa=1 ! only one size with side reactions
c     If side reaction 3 does not involve absorption or desorption,
c     remove c(kj3,j) and modify the appropriate Jacobian.  This feature
c     is included to treat hydrogen absorption/desorption in the NiMH system.
                b(ii2,kj3)=1.0d0/Rad1pa(mpa)
                b(ii2,2+mpa)=ai2(1,mpa)/rr*omi
                b(ii2,kj)=1.0d0/Rad1pa(mpa)
                g(ii2)=ai2(1,mpa)*xt(2+mpa,j,kk-1)/rr*omi-sumpa(mpa,j)-
     &       ai2(1,mpa)/rr*omi*c(2+mpa,j)-(c(kj,j)+c(kj3,j))/Rad1pa(mpa)
              else !no side reaction 3
                b(ii2,kj)=1.d0
                g(ii2)=-c(kj,j)
                do mpa=1,npa ! Calculate g, no side reaction
                  b(ii2,2+mpa)=ai2(1,mpa)/rr*omi*Rad1pa(mpa)
     &            *area1pa(mpa)/area1
                  g(ii2)=g(ii2)+area1pa(mpa)/area1*(ai2(1,mpa)*(xt(2+mpa
     &            ,j,kk-1)-c(2+mpa,j))/rr*omi-sumpa(mpa,j))*Rad1pa(mpa)
                enddo ! mpa

              endif ! end side (and no side) reactions

c Special section for VDC -anode
            elseif (mvdc1.eq.1) then
              mpa=1 ! only one size for variable solid-phase diffusion coefficient
              call vardc(j,jcount,kk)
              b(ii2,2+mpa)=1.0D0
              g(ii2)=css(j,nnj)-xx(2+mpa,j)
            endif ! end calculate b's
          endif ! end nonfoil anode

        elseif (j.ge.n1+n2) then !cathode
          if (mvdc3.eq.0) then ! superposition
            if (nside.eq.3) then ! side reaction with hydrogen
              mpa=1 ! only one size with side reactions
c     If side reaction 3 does not involve absorption or desorption,
c     remove c(kj3,j) and modify the appropriate Jacobian.  This feature
c     is included to treat hydrogen absorption/desorption in the NiMH system.
              b(ii2,2+mpa)=ai(1,mpa)/rr*omi
              b(ii2,kj)=1.0d0/Rad3pa(mpa)
              b(ii2,kj3)=1.0d0/Rad3pa(mpa)
              g(ii2)=ai(1,mpa)*xt(2+mpa,j,kk-1)/rr*omi-sumpa(mpa,j)-
     &      ai(1,mpa)/rr*omi*c(2+mpa,j)-(c(kj,j)+c(kj3,j))/Rad3pa(mpa)
            else ! no side reaction 3
              b(ii2,kj)=1.d0
              g(ii2)=-c(kj,j)
              do mpa=1,npa
                b(ii2,2+mpa)=ai(1,mpa)/rr*omi*Rad3pa(mpa)
     &          *area3pa(mpa)/area3
                g(ii2)=g(ii2)+area3pa(mpa)/area3*(ai(1,mpa)*(xt(2+mpa,j
     &            ,kk-1)-c(2+mpa,j))/rr*omi-sumpa(mpa,j))*Rad3pa(mpa)
              enddo !mpa

            endif !end side reactions

c Special section of VDC - cathode
         elseif (mvdc3.eq.1) then !calculate b's
           mpa=1 ! only one size for variable solid-phase diffusion coefficient
           call vardc(j,jcount,kk)
           b(ii2,2+mpa)=1.0D0
           g(ii2)=css(j,nnj)-xx(2+mpa,j)
         endif ! end calculate b's

        else  !separator
          b(ii2,4+npa)=1.0d0
          g(ii2)=-c(4+npa,j) ! set solid flux = 0 in separator
        endif

      else ! Impedance section

      mpa=1 ! only one size for impedance
        if(j.eq.1 .and. n1.eq.0) then
c No solid concentration for foil electrode:
        b(ii2,2+mpa)=1.0d0
        g(ii2)=-c(2+mpa,j)
        b(ii2+n,2+mpa+n)=1.0d0
        g(ii2+n)=-c(2+mpa+n,j)
        g(ii2)=0.0d0
         elseif(j.gt.n1+1 .and. j.lt.n1+n2) then ! separator
        b(ii2,2+mpa)=1.0d0
        b(ii2+n,2+mpa+n)=1.0d0
        g(ii2)=0.d0
         else ! anode and cathode
        if(j.eq.1) then
        Rad=Rad1
        dfs=dfs1
        elseif(j.eq.n1+n2) then
        Rad=Rad3
        dfs=dfs3
        endif
c Solution of diffusion equation w/Laplace gives:
        b(ii2,kj)=-1.0d0
        b(ii2+n,kj3)=-1.0d0
c New section added April 4, 1998, JSN
        arg=Rad*(2.d0*omega(ji)/dfs)**0.5d0
        e1=expg(-arg)
        e2=e1**2
        if(arg.lt.0.22d0) then
        denom=2.d0*e1*arg**2*(1.d0+arg**4/3.6d2+arg**8/1.8144d6)
        diff2=2.d0*e1*(2.d0*arg+arg**5/6.d1+arg**9/1.99584d7*110.d0)
        diff3=2.d0*e1*(arg**3/3.d0+arg**7/2.52d3+arg**11/1.99584d7)
        Ysr=dfs/Rad*(arg**4/1.8d2+arg**8/4.536d5+arg**12/7.2648576d9)
     %  /1.d0/(1.d0+arg**4/3.6d2+arg**8/1.8144d6)
        else
        denom=1.d0+e2-2.d0*e1*dcos(arg)
        diff2=1.d0-e2+2.d0*e1*dsin(arg)
        diff3=1.d0-e2-2.d0*e1*dsin(arg)
        Ysr=-dfs/Rad*(1.d0-arg/2.d0*diff2/denom)
        endif
        Ysi=dfs/Rad*arg/2.d0*diff3/denom
        real1=-Ysr
        comp1=Ysi
        b(ii2,2+mpa)=real1
        b(ii2+n,2+mpa+n)=real1
        b(ii2,2+mpa+n)=comp1
        b(ii2+n,2+mpa)=-comp1
        g(ii2)=0.0d0
         endif
      g(ii2+n)=0.0d0
      endif ! impedance section (equation ii2)

c     ________________________________________________

c Equation imb2 for variable kS2
c Material balance on the second side-reaction species.
c This is set up for a balance on hydrogen for the NiMH system.

      if (nside.ge.2) then

      sf=1.0d0/2.0d0  !stoichiometric factor for H2 reaction

      termn_s2(j)=0.d0

c aream the area for use in mass balance

      if (j.le.n1+1) aream=area1
      if (j.ge.n1+n2) aream=area3

      fac=1.d0
      if(j.eq.n1+2 .and. n1 .gt. 0)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2+1) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
      epn=ep2+epp2

      hn=h2  !box to left of point.
      epn=epg2
      if (j.le.n1+1) then
      hn=h1
      epn=epg1
      elseif (j.gt.n1+n2) then
      hn=h3
      epn=epg3
      endif

      if (j.eq.n1+n2) sf=0.0d0 !half box w/o generation
      if (j.eq.n1+2) sf=0.0d0 !half box w/0 generation

      if (j.gt.1) then !box to left of point
      termn_s2(j)=-(dh2(j)+fac*dh2(j-1))*(c(kS2,j)-c(kS2,j-1))/hn/2.0d0
     &+sf*hn*aream*vol*(0.375d0*(c(kj3,j)+c(kj2,j))+0.125d0*
     &(c(kj3,j-1)+c(kj2,j-1)))
      g(imb2)=-epn*hn*(0.375d0*(c(kS2,j)-xt(kS2,j,kk-1+kadd))
     &+0.125d0*(c(kS2,j-1)-xt(kS2,j-1,kk-1+kadd)))/rr
      b(imb2,kS2)=epn*hn*0.375d0/rr+(dh2(j)+
     &fac*dh2(j-1))/hn/4.0d0+ddh2(j)*(c(kS2,j)-c(kS2,j-1))/hn/4.0d0
      a(imb2,kS2)=-(dh2(j)+fac*dh2(j-1))/hn/4.0d0+hn*epn*0.125d0/rr
     &+fac*ddh2(j-1)*(c(kS2,j)-c(kS2,j-1))/hn/4.0d0
      b(imb2,kj3)=-0.375d0*hn*sf*aream*vol/2.0d0
      b(imb2,kj2)=-0.375d0*hn*sf*aream*vol/2.0d0
      a(imb2,kj3)=-0.125d0*hn*sf*aream*vol/2.0d0
      a(imb2,kj2)=-0.125d0*hn*sf*aream*vol/2.0d0
      endif

      sf=1.0d0/2.0d0

      fac=1.d0
      if(j.eq.n1+1 .and. n1 .gt. 0)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
      epn=ep2+epp2

      hn=h2 !deal with box to right of point.
      epn=epg2
      if (j.lt.n1+1) then
      hn=h1
      epn=epg1
      elseif (j.ge.n1+n2) then
      hn=h3
      epn=epg3
      endif

      if (j.eq.n1+1) sf=0.0d0 !half box w/o generation
      if (j.eq.n1+n2-1) sf=0.0d0 !half box w/o generation

      if (j.lt.nj) then !box to right of point
      termn_s2(j)=termn_s2(j)-(fac*dh2(j)+dh2(j+1))*
     &(c(kS2,j)-c(kS2,j+1))/hn/2.0d0+
     &sf*hn*aream*vol*(0.375d0*(c(kj3,j)+c(kj2,j))+0.125d0*
     &(c(kj3,j+1)+c(kj2,j+1)))
      g(imb2)=g(imb2)-epn*hn*(0.375d0*(c(kS2,j)-xt(kS2,j,kk-1+kadd))
     &+0.125d0*(c(kS2,j+1)-xt(kS2,j+1,kk-1+kadd)))/rr

      b(imb2,kS2)=b(imb2,kS2)+epn*hn*0.375d0/rr+(fac*dh2(j)+dh2(j+1))
     &/hn/4.0d0+fac*ddh2(j)*(c(kS2,j)-c(kS2,j+1))/hn/4.0d0
      b(imb2,kj3)=b(imb2,kj3)-sf*hn*aream*vol*0.375d0/2.0d0
      b(imb2,kj2)=b(imb2,kj2)-sf*hn*aream*vol*0.375d0/2.0d0
      d(imb2,kS2)=-(fac*dh2(j)+dh2(j+1))/hn/4.0d0+ddh2(j+1)*
     &(c(kS2,j)-c(kS2,j+1))/hn/4.0d0+epn*hn*0.125d0/rr
      d(imb2,kj3)=-sf*hn*aream*vol*0.125d0/2.0d0
      d(imb2,kj2)=-sf*hn*aream*vol*0.125d0/2.0d0
      endif

      g(imb2)=g(imb2)+(termn_s2(j)+term_s2(j))/2.0d0

      endif

c     ________________________________________________

c Equation imb1 for variable kS1
c Material balance on the first side-reaction species.
c This is set up for a balance on oxygen for the NiMH system.

      if (nside.ge.1) then
      sf=-1.0d0/4.0d0  !stoichiometric factor for O2 reaction

	if (nprop.ne.13) then ! Need diff. coeff. for Li-ion side reaction
      do2(j)=1.0d-5 !gas-phase oxygen diffusion coefficient
      ddo2(j)=0.0d0  ! first derivative wrt O2 concentration
	endif

      termn_s1(j)=0.d0

      if (j.le.n1+1) aream=area1
      if (j.ge.n1+n2) aream=area3

      fac=1.d0
      if(j.eq.n1+2 .and. n1 .gt. 0)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2+1) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
      epn=ep2+epp2

      hn=h2  !box to left of point.
c For liquid-phase side reaction see User's Manual.
c     epn=ep2
      epn=epg2
      if (j.le.n1+1) then
      hn=h1
c     epn=ep1
      epn=epg1
      elseif (j.gt.n1+n2) then
      hn=h3
c     epn=ep3
      epn=epg3
      endif

      if (j.eq.n1+n2) sf=0.0d0 !half box w/o generation
      if (j.eq.n1+2) sf=0.0d0 !half box w/0 generation
      if (j.gt.1) then !box to left of point
      termn_s1(j)=-(do2(j)+fac*do2(j-1))*(c(kS1,j)-c(kS1,j-1))/hn/2.0d0
     &+sf*hn*aream*vol*(0.375d0*c(kj1,j)+0.125d0*c(kj1,j-1))
      g(imb1)=-epn*hn*(0.375d0*(c(kS1,j)-xt(kS1,j,kk-1+kadd))
     &+0.125d0*(c(kS1,j-1)-xt(kS1,j-1,kk-1+kadd)))/rr
      b(imb1,kS1)=epn*hn*0.375d0/rr+(do2(j)+fac*do2(j-1))/hn/4.0d0
     &+ddo2(j)*(c(kS1,j)-c(kS1,j-1))/hn/4.0d0
      a(imb1,kS1)=-(do2(j)+fac*do2(j-1))/hn/4.0d0+epn*hn*0.125d0/rr
     &+fac*ddo2(j-1)*(c(kS1,j)-c(kS1,j-1))/hn/4.0d0
      b(imb1,kj1)=-0.375d0*hn*sf*aream*vol/2.0d0
      a(imb1,kj1)=-0.125d0*hn*sf*aream*vol/2.0d0
      endif

      sf=-1.0d0/4.0d0

      fac=1.d0
      if(j.eq.n1+1 .and. n1 .gt. 0)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
      epn=ep2+epp2

      hn=h2 !deal with box to right of point.
c For liquid-phase side reaction see User's Manual.
c     epn=ep2
      epn=epg2
      if (j.lt.n1+1) then
      hn=h1
c     epn=ep1
      epn=epg1
      elseif (j.ge.n1+n2) then
      hn=h3
c     epn=ep3
      epn=epg3
      endif

      if (j.eq.n1+1) sf=0.0d0 !half box w/o generation
      if (j.eq.n1+n2-1) sf=0.0d0 !half box w/o generation

      if (j.lt.nj) then !box to right of point
      termn_s1(j)=termn_s1(j)-(fac*do2(j)+do2(j+1))*(c(kS1,j)
     &-c(kS1,j+1))
     &/hn/2.0d0+sf*hn*aream*vol*(0.375d0*c(kj1,j)+0.125d0*c(kj1,j+1))
      g(imb1)=g(imb1)-epn*hn*(0.375d0*(c(kS1,j)-xt(kS1,j,kk-1+kadd))
     &+0.125d0*(c(kS1,j+1)-xt(kS1,j+1,kk-1+kadd)))/rr

      b(imb1,kS1)=b(imb1,kS1)+epn*hn*0.375d0/rr+(fac*do2(j)+do2(j+1))
     &/hn/4.0d0+fac*ddo2(j)*(c(kS1,j)-c(kS1,j+1))/hn/4.0d0
      b(imb1,kj1)=b(imb1,kj1)-sf*hn*aream*vol*0.375d0/2.0d0
      d(imb1,kS1)=-(fac*do2(j)+do2(j+1))/hn/4.0d0+epn*hn*0.125d0/rr-
     &ddo2(j+1)*(c(kS1,j)-c(kS1,j+1))/hn/4.0d0
      d(imb1,kj1)=-sf*hn*aream*vol*0.125d0/2.0d0
      endif

      g(imb1)=g(imb1)+(termn_s1(j)+term_s1(j))/2.0d0

      endif

  199 if(imp.eq.1) g(12)=0.0d0

c     ________________________________________________
c     Equation ip2, Ohm's law in the liquid phase
c     There are separate versions of this equation for a Li-ion
c     system and for the NiMH system.
      if(j.le.n1) then
      h=h1
      else if(j.lt.n1+n2) then
      h=h2
      else
      h=h3
      endif
      fac=1.0d0
      if(j.eq.n1+1 .and. n1 .gt. 0)
     1fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2)
     1fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)
c
         if(j.eq.nj) then
      b(ip2,kp2)=1.0d0
      g(ip2)=-c(kp2,j) !boundary condition on potential
      if (imp.eq.1) g(ip2)=0.0d0
      go to 12
         endif

      if (nprop.eq.13) then 
c Ohm's law in solution for NiMH different from Lithium
      dcf=(c(1,j+1)-c(1,j))/h 
      r1=(c(1,j+1)+c(1,j))/2.0d0
      r4=c(ki2,j)
      p1=1.0d0-(tm(j)+tm(j+1))/2.0d0 !different from used in mass balance
      p2=(fac*cd(j)+cd(j+1))/2.0d0
      p3=(fac*dcd(j)+dcd(j+1))/2.0d0
      H2OMW=18.016D0
      PKOHMW=56.11D0
      c1=r1/1.0d6
      fu= 1.004D0-36.23D0*c1**0.5+1374.3D0*c1-17850.7d0*c1**1.5d0+
     &55406.0D0*c1**2d0+7.16856D05*c1**2.5d0
      d_fu=-18.115D0*c1**(-0.5)+1374.3D0-2.6776D04*c1**0.5
     &+1.10812D05*c1 + 1.7921D06*c1**1.5
      d2_fu=9.0575D0/c1**1.5-13388.0D0/c1**0.5+1.10812D05+
     & 2.6882D06*c1**0.5
      d_fu=d_fu/1.0d6
      d2_fu=d2_fu/1.0d12
      densw=1.001D0+47.57D0*c1-776.22D0*c1**2
      dens_1d = 47.57D0 - 1552.44D0*c1
      dens_2d= -1552.44D0
      densw=densw*1.0d6
      dens_1d=dens_1d*1.0d0
      dens_2d=dens_2d/1.0d6

      c_solv=(densw-c(1,j)*PKOHMW)/H2OMW

      U1=densw-c(1,j)*PKOHMW
      U2=dens_1d-PKOHMW
      U3=-dcf*((p1*U1+H2OMW*c(1,j))*((fu*d2_fu-d_fu**2.0D0)/
     &(fu**2.0D0)-1.0D0/c(1,j)**2.0D0)+(d_fu/fu+1.0D0/xx(1,j))
     &*(p1*U2+H2OMW))
      U4=-p1*(d_fu*densw/fu +densw/xx(1,j)) -(H2OMW-p1*PKOHMW)*
     &(d_fu*c(1,j)/fu + 1.0D0)

      d(ip2,1)=-frt*c(ki2,j)/p2/p2*p3/2.0d0
      d(ip2,1)=d(ip2,1)+(p1+r1/c_solv)*(d_fu/fu+1/r1)*(1.0d0/h)
      d(ip2,1)=d(ip2,1)+(p1+r1/c_solv)*dcf*(-d_fu/fu/fu*d_fu*0.5d0
     &+0.5d0/fu*d2_fu-0.5d0/r1/r1)
      d(ip2,1)=d(ip2,1)+(d_fu/fu+1/r1)*dcf*(-0.5d0*r1*H2OMW/U1**2.d0*U2
     &+0.5*H2OMW/U1)
      d(ip2,1)=-d(ip2,1)

      b(ip2,1)=frt*c(ki2,j)/p2/p2*p3/2.0d0
      b(ip2,1)=b(ip2,1)+(p1+r1/c_solv)*(d_fu/fu+1/r1)*(1.0d0/h)
      b(ip2,1)=b(ip2,1)+(p1+r1/c_solv)*dcf*(-d_fu/fu/fu*d_fu*0.5d0
     &+0.5d0/fu*d2_fu-0.5d0/r1/r1)
      b(ip2,1)=b(ip2,1)+(d_fu/fu+1/r1)*dcf*(-0.5d0*r1*H2OMW/U1**2.d0*U2
     &+0.5*H2OMW/U1)
      b(ip2,1)=-b(ip2,1)
      d(ip2,kp2)=frt/h
      b(ip2,kp2)=-frt/h
      b(ip2,ki2)=frt/p2
      g(ip2)=-frt*(c(kp2,j+1)-c(kp2,j))/h-frt/p2*c(ki2,j)
     &-(p1+r1/c_solv)*(d_fu*1/fu+1/r1)*dcf

      else ! for Lithium chemistry

      dcf=(xx(1,j+1)-xx(1,j))/h *2.0d0 !factor of 2 added
      r1=(xx(1,j+1)+xx(1,j))/2.0d0
      r4=xx(ki2,j)
      p1=(tm(j)+tm(j+1))/2.0d0
      p2=(fac*cd(j)+cd(j+1))/2.0d0
      p4=(dfu(j)+dfu(j+1))/2.0d0
      d(ip2,1)=(1.0d0-p1)*(1.0d0/r1+p4)/h *2.0d0 !factor of 2 added
      b(ip2,1)=-d(ip2,1)+((1.0d0-p1)*(d2fu(j)-1.0d0/r1/r1)*dcf
     &-(1.0d0/r1+p4)*dcf*dtm(j)+frt*r4*fac*dcd(j)/p2/p2)/2.0d0
      d(ip2,1)=d(ip2,1)+((1.0d0-p1)*(d2fu(j+1)-1.0d0/r1/r1)*dcf
     &-(1.0d0/r1+p4)*dcf*dtm(j+1)+frt*r4*dcd(j+1)/p2/p2)/2.0d0
      d(ip2,kp2)=-frt/h
      b(ip2,kp2)=frt/h
      b(ip2,ki2)=-frt/p2
      g(ip2)=-(1.0d0-p1)*(1.0d0/r1+p4)*dcf+frt*(c(kp2,j+1)-c(kp2,j))/h
     &+frt/p2*c(ki2,j)
      endif

      if (imp.eq.1) g(ip2)=0.0d0
   12 if(imp.eq.1) g(ip2+n)=0.0d0

c     ________________________________________________
c     Equation 3, Butler-Volmer kinetics
c     Note that a zero and a nonzero time step are treated differently.
c     Also note that the equation is written differently depending on
c     whether superposition or a discretized version of the solid-phase
c     material balance is being used.

c Define the individual pore-wall fluxes.  These are necessary for
c multiple particle sizes in order to treat properly a film resistance,
c and for a rigorous energy balance.

      izfl=0
      if (rr.eq.0.0d0) izfl=1 ! need a special flag for zero timestep

      do ii=1,nj ! counter index loop
      do mpa=1,npa
      if (ii.le.n1+1) then
      pwfindiv(mpa,ii)=(ai2(1,mpa)*(xt(2+mpa,ii,kk-1+kadd)-c(2+mpa,ii))
     &/rr*omi-sumpa(mpa,ii))*Rad1pa(mpa)
      if (mvdc1.eq.1) pwfindiv(mpa,ii)=xx(kj,ii)
      pwfindiv_d(mpa,ii)=-Rad1pa(mpa)*ai2(1,mpa)/rr*omi
      elseif (ii.ge.n1+n2) then
      pwfindiv(mpa,ii)=(ai(1,mpa)*(xt(2+mpa,ii,kk-1+kadd)-c(2+mpa,ii))
     &/rr*omi-sumpa(mpa,ii))*Rad3pa(mpa)
      if (mvdc3.eq.1) pwfindiv(mpa,ii)=xx(kj,ii)
      pwfindiv_d(mpa,ii)=-Rad3pa(mpa)*ai(1,mpa)/rr*omi
      else
      pwfindiv(mpa,ii)=0.0d0
      pwfindiv_d(mpa,ii)=0.0d0
      endif
      enddo !mpa
      enddo

      if(rr.eq.0.d0) then ! do a sum of B-V equations to get j
        if(j.gt.n1+1 .and. j.lt.n1+n2) then
          b(ii2,kj)=1.0d0
          g(ii2)=-c(kj,j) ! zero transfer current in separator
        else
          do mpa=1,npa
            call ekin(j,kk,0,0.0d0,mpa,mvdc1,mvdc3)
            if(j.eq.n1+n2 .and. mpa.eq.1) Uc=g0
            if(j.eq.1 .and. mpa.eq.1)     Ua=g0
            if(j.le.n1+1) then
              g(ii2)=g(ii2)+g(2+mpa)*area1pa(mpa)/area1
                do ii=1,n
                  b(ii2,ii)=b(ii2,ii)+b(2+mpa,ii)*area1pa(mpa)/area1
                enddo ! ii
            elseif(j.ge.n1+n2) then
              g(ii2)=g(ii2)+g(2+mpa)*area3pa(mpa)/area3
              do ii=1,n
                b(ii2,ii)=b(ii2,ii)+b(2+mpa,ii)*area3pa(mpa)/area3
              enddo ! ii
            endif
          enddo ! mpa
        endif
        do mpa=1,npa
          do ii=1,n
            b(2+mpa,ii)=0.d0
          enddo
          b(2+mpa,2+mpa)=1.0d0
          g(2+mpa)=xt(2+mpa,j,kk-1+kadd)-c(2+mpa,j) ! fix solid concentration
        enddo ! mpa

      else !rr>0

      if((mvdc1.eq.1.and.j.le.n1+1).or.
     &(mvdc3.eq.1.and.j.ge.n1+n2)) then  ! Variable diffusion coefficient
      mpa=1 ! only one size for variable solid-phase diffusion coefficient
      if(j.gt.n1+1 .and. j.lt.n1+n2) then
      b(2+mpa,kj)=1.0d0
      g(2+mpa)=-c(kj,j) ! zero transfer current in separator
      else
      call ekin(j,kk,0,0.0d0,1,mvdc1,mvdc3)
       if(j.eq.n1+n2) Uc=g0
       if(j.eq.1)     Ua=g0
      endif

      else !using superposition

      do mpa=1,npa
        if(j.gt.n1+1 .and. j.lt.n1+n2) then !separator
          if (imp.eq.0) then
            b(2+mpa,2+mpa)=1.0d0
            g(2+mpa)=-c(2+mpa,j) ! zero solid concentration in separator
          else
            b(2+mpa,kj)=1.0d0
            g(2+mpa)=-c(kj,j)
          endif
        else
          call ekin(j,kk,0,0.0d0,mpa,mvdc1,mvdc3)
          if(j.eq.n1+n2 .and. mpa.eq.1) Uc=g0
          if(j.eq.1 .and. mpa.eq.1)     Ua=g0

          if(j.le.n1+1) then
            g(2+mpa)=g(2+mpa)+xx(kj,j)-Rad1pa(mpa)*(-sumpa(mpa,j)
     &      +ai2(1,mpa)/rr*omi*(xt(2+mpa,j,kk-1+kadd)-c(2+mpa,j)))
            b(2+mpa,2+mpa)=b(2+mpa,2+mpa)-Rad1pa(mpa)*ai2(1,mpa)/rr*omi
            if (imp.eq.0) b(2+mpa,kj)=b(2+mpa,kj)-1.d0
          elseif(j.ge.n1+n2) then
            g(2+mpa)=g(2+mpa)+xx(kj,j)-Rad3pa(mpa)*(-sumpa(mpa,j)
     &      +ai(1,mpa)/rr*omi*(xt(2+mpa,j,kk-1+kadd)-c(2+mpa,j)))
            b(2+mpa,2+mpa)=b(2+mpa,2+mpa)-Rad3pa(mpa)*ai(1,mpa)/rr*omi
            if (imp.eq.0) b(2+mpa,kj)=b(2+mpa,kj)-1.d0
          endif
        endif
      enddo ! mpa
      endif

      endif

      mpa=1 ! only one size for impedance
      if (imp.eq.1) g(2+mpa)=0.0d0
      if(imp.eq.1) g(2+mpa+n)=0.0d0

c ________________________________________________
c    Equation iSk1.  Side reaction 1 kinetics
c    This is set up for oxygen evolution (Tafel) and 
c    recombination (concentration driven) for the NiMH system.
c    We also have a sample side reaction for solvent reduction/oxidation
c    in a Li-ion system.  Tafel kinetics is used at each electrode.
c    For this reaction, we recommend setting rksa1=1.0d-9 and rksac1=1.0d-9
c    in the input file.

      if (nside.ge.1) then

      rksa1_save=rksa1
      rksc1_save=rksc1
      ti0n=dexp((Ebarks1a)*(t-298.15d0)/(t*298.15d0))
      ti0p=dexp((Ebarks1c)*(t-298.15d0)/(t*298.15d0))
      rksa1=rksa1*ti0n
      rksc1=rksc1*ti0p
      

      alphasO2=1.17d0
      UsO2=0.25d0

      if (j.le.n1+1) then  !recombination at the anode
      g(iSk1)=-c(kj1,j)+rksa1*c(kS1,j)
      b(iSk1,kj1)=1.0d0
      b(iSk1,kS1)=-rksa1

c Sample kinetics for side reaction in a Li-ion system: comment out
c the above and uncomment the below.  We use an alpha value of 0.5
c and use a value of 1.0 V for the "equilibrium" potential of
c solvent reduction

c      g(iSk1)=-c(kj1,j)+rksa1*expg(0.5d0*frt*(c(kp1,j)-c(kp2,j)-1.d0))
c      b(iSk1,kj1)=1.0d0
c      b(iSk1,kp1)=-rksa1*0.5d0*frt*expg(0.5d0*frt*(c(kp1,j)-c(kp2,j)
c    &-1.d0))
c      b(iSk1,kp2)=-b(iSk1,kp1)

      elseif (j.ge.n1+n2) then !evolution at the cathode
      g(iSk1)=-c(kj1,j)-rksc1*expg(alphasO2*frt*(c(kp1,j)-c(kp2,j)
     &-UsO2))
      b(iSk1,kj1)=1.0d0
      b(iSk1,kp1)=rksc1*alphasO2*frt*expg(alphasO2*frt*
     &(c(kp1,j)-c(kp2,j)-UsO2))
      b(iSk1,kp2)=-b(iSk1,kp1)

c Sample kinetics for side reaction in a Li-ion system: comment out
c the above and uncomment the below.  We use an alpha value of 0.5
c and use a value of 4.4 V for the "equilibrium" potential of 
c solvent oxidation

c      g(iSk1)=-c(kj1,j)-rksc1*expg(0.5d0*frt*(c(kp1,j)-c(kp2,j)-4.d0))
c      b(iSk1,kj1)=1.0d0
c      b(iSk1,kp1)=rksc1*0.5d0*frt*expg(0.5d0*frt*(c(kp1,j)-
c     1c(kp2,j)-4.d0))
c      b(iSk1,kp2)=-b(iSk1,kp1)

      else !in the separator
      g(iSk1)=c(kj1,j)
      b(iSk1,kj1)=-1.0d0
      endif

      else
      g(iSk1)=0.0d0
      endif

      rksa1=rksa1_save
      rksc1=rksc1_save

c ________________________________________________
c    Equation iSk2.  Side reaction 2 kinetics 
c    This is set up for hydrogen evolution from the negative electrode
c    of the NiMH system.

      if (nside.ge.2) then
      alphasH2=1.0d0
      UsH2=-0.98d0

      rksa2_save=rksa2
      rksc2_save=rksc2
      ti0n=dexp((Ebarks3a)*(t-298.15d0)/(t*298.15d0))
      ti0p=dexp((Ebarks3c)*(t-298.15d0)/(t*298.15d0))
      rksa2=rksa2*ti0n
      rksc2=rksc2*ti0p

      if (j.le.n1+1) then  
c H2 evolution from the negative electrode only.
      g(iSk2)=c(kj2,j)-rksa2*expg(-alphasH2*frt*(c(kp1,j)
     &-c(kp2,j)-UsH2))
      b(iSk2,kj2)=-1.0d0
      b(iSk2,kp1)=-rksa2*alphasH2*frt*
     &expg(-alphasH2*frt*(c(kp1,j)-c(kp2,j)-UsH2))
      b(iSk2,kp2)=-b(iSk2,kp1)
      else !in the separator or the positive
      g(iSk2)=c(kj2,j) 
      b(iSk2,kj2)=-1.0d0
      endif

      else
      if(imp.eq.1) g(3+n)=0.0d0
      endif

      rksa2=rksa2_save
      rksc2=rksc2_save

c ________________________________________________
c    Equation iSk3.  Side reaction 3 kinetics
c    This is set up for hydrogen absorption/desorption in the
c    negative electrode of the NiMH system.
c Positive j is H leaving the solid phase

      if (nside.eq.3) then

      rksa3_save=rksa3
      rksc3_save=rksc3
      ti0e=dexp((Ebarks2a)*(t-298.15d0)/(t*298.15d0))
      ti0r=dexp((Ebarks2c)*(t-298.15d0)/(t*298.15d0))
      rksa3=rksa3*ti0e
      rksc3=rksc3*ti0r

      mpa=1 ! only one size with side reactions
      if (j.le.n1+1) then  
c H2 absorption, using an equilibrium isotherm from the literature

      g(iSk3)=c(kj3,j)+rksa3*(c(kS2,j)*t*8.206d-5-
     &(5.d0+4.9999d0*dtanh((c(2+mpa,j)/ct1-1.04d0)/0.11d0)))
      b(iSk3,kS2)=-rksa3*t*8.2d-5
      b(iSk3,2+mpa)=-rksa3*4.9999d0*
     &((dcosh((c(2+mpa,j)/ct1-1.04d0)/0.11d0))**(-2))*1.0d0/0.11d0/ct1
      b(iSk3,kj3)=-1.0d0

      elseif (j.ge.n1+n2) then !reaction at the positive electrode

      g(iSk3)=c(kj3,j)
      b(iSk3,kj3)=-1.0d0

      else !in the separator
      g(iSk3)=c(kj3,j)
      b(iSk3,kj3)=-1.0d0
      endif

      else
      mpa=1 ! only one size for impedance
      if(imp.eq.1) g(5+n)=0.0d0
      endif

      rksa3=rksa3_save
      rksc3=rksc3_save

c     ________________________________________________
c     Equation ip1, Ohm's law in solid
c
      if(j.eq.1) then
      if (n1 .eq. 0) then ! foil electrode
      g(ip1) = cur - fc*c(kj,j)
      b(ip1,kj) = fc
      if(imp.eq.1) then
      g(ip1) = 1.d0
c     area1 is mult by L here so aL is 1 in double-layer charging term:
      b(ip1,kp2+n)=omega(ji)*capp1
      b(ip1,kp1+n)=-omega(ji)*capp1
      b(ip1+n,kp2)=-omega(ji)*capp1
      b(ip1+n,kp1)=omega(ji)*capp1
      endif ! impedance section (j=1, a)
      else
      b(ip1,kp1)=-1.0d0/h1 ! not a foil electrode
      d(ip1,kp1)= 1.0d0/h1
      b(ip1,ki2)=-1.0d0/sig1
c     next statement requires that cur be the current in the separator
      g(ip1)=-cur/sig1+c(ki2,j)/sig1-(c(kp1,j+1)-c(kp1,j))/h1
      if(imp.eq.1) g(ip1)=-1.0d0/sig1
      endif
      elseif(j.lt.n1+1) then ! anode
      b(ip1,kp1)=-1.0d0/h1
      d(ip1,kp1)= 1.0d0/h1
      b(ip1,ki2)=-1.0d0/sig1
      g(ip1)=-cur/sig1+c(ki2,j)/sig1-(c(kp1,j+1)-c(kp1,j))/h1
      if(imp.eq.1) g(ip1)=-1.0d0/sig1
      elseif(j.eq.n1+1) then
c     This is the current boundary condition.
            if(mcL.eq.-3) then ! constant load
      g(ip1)=(xx(nx+1,j+1)-xx(nx+1,j))/xx(ki2,j)-cuL-RG !WHT 1-11-08
      d(ip1,nx+1)=-1.d0/xx(ki2,j)
      b(ip1,nx+1)=1.d0/xx(ki2,j)
      b(ip1,ki2)=(xx(nx+1,j+1)-xx(nx+1,j))/xx(ki2,j)**2
            elseif(mcL.eq.0) then ! constant potential
      g(ip1)=xx(nx+1,j+1)-xx(nx+1,j)-cuL-RG*xx(ki2,j) !WHT 1-11-08
      d(ip1,nx+1)=-1.d0
      b(ip1,nx+1)=1.d0
      b(ip1,ki2)=RG !WHT 1-11-08
            elseif(mcL.eq.-2) then ! constant power
      g(ip1)=(xx(nx+1,j+1)-xx(nx+1,j))*xx(ki2,j)-cuL-RG*xx(ki2,j)**2 !WHT 1-11-08
      d(ip1,nx+1)=-xx(ki2,j)
      b(ip1,nx+1)=xx(ki2,j)
      b(ip1,ki2)=-(xx(nx+1,j+1)-xx(nx+1,j))+2.d0*RG*xx(ki2,j) !WHT 1-11-08
            elseif(mcL.eq.-4) then ! maximum power
      g(ip1)=xx(nx+1,j+1)-xx(nx+1,j)-xx(ki2,j)*Rint !RG*something
      d(ip1,nx+1)=-1.d0
      b(ip1,nx+1)=1.d0
      b(ip1,ki2)=Rint
      cu0=(xx(nx+1,j+1)-xx(nx+1,j))*xx(ki2,j)
            else ! constant current
      b(ip1,ki2)=1.0d0
      g(ip1)=cur-xx(ki2,j)
            endif

      if(imp.eq.1) g(ip1)=1.0d0
      elseif(j.lt.n1+n2) then
      b(ip1,kp1)=1.0d0 ! separator, set phi1=0
      g(ip1)=-c(kp1,j)
      if(imp.eq.1) g(ip1)=0.0d0
      elseif(j.lt.nj) then ! cathode
      d(ip1,kp1)=1.0d0/h3
      b(ip1,kp1)=-1.0d0/h3
      b(ip1,ki2)=-1.0d0/sig3
      g(ip1)=-cur/sig3+c(ki2,j)/sig3-(c(kp1,j+1)-c(kp1,j))/h3
      if(imp.eq.1) g(ip1)=-1.0d0/sig3
      else ! j=nj
      b(ip1,ki2)=1.0d0
      g(ip1)=-c(ki2,j) !  i2 is no longer used at j=nj
      if(imp.eq.1) g(ip1)=0.d0
      endif
      if(imp.eq.1) g(12)=0.0d0

c     ________________________________________________

c     Equation n+1.  Carry end potentials.  This equation is added to
c     permit handling constant load, etc. without an extra loop.
      if (imp.eq.0) then
      if(j.eq.1) then
      g(nx+1)=xx(kp1,1)-xx(nx+1,1) ! pick up negative potential
      b(nx+1,kp1)=-1.d0
      b(nx+1,nx+1)=1.d0
      elseif(j.eq.nj) then
      g(nx+1)=xx(kp1,nj)-xx(nx+1,nj) ! pick up positive potential
      b(nx+1,kp1)=-1.d0
      b(nx+1,nx+1)=1.d0
      elseif(j.le.n1+1) then
      g(nx+1)=xx(nx+1,j)-xx(nx+1,j-1) ! carry negative potential
      b(nx+1,nx+1)=-1.d0
      a(nx+1,nx+1)=1.d0
      else ! if(j.ge.n1+n2) then
      g(nx+1)=xx(nx+1,j)-xx(nx+1,j+1) ! carry positive potential
      b(nx+1,nx+1)=-1.d0
      d(nx+1,nx+1)=1.d0
      endif
      endif

c     ________________________________________________
c
c     Equation i2div.  Current balance

      if(j.gt.n1+1 .and. j.lt.n1+n2) then ! separator
        g(i2div)=xx(ki2,j-1)-xx(ki2,j)
        if (imp.eq.0) a(i2div,ki2)=-1.d0
        b(i2div,ki2)=1.d0
        if (imp.eq.1) then ! impedance section
c     Set perturbation current = 1 purely real:
          b(i2div,ki2)=1.0d0
          g(i2div)=1.0d0
        endif ! end of impedance section (n1+1<j<n1+n2)
c
      elseif(j.eq.1 .and. n1.eq.0) then ! foil electrode
        g(i2div) = cur - c(ki2,j)
        b(i2div,ki2) = 1.0d0
        if(imp.eq.1) g(i2div)=1.d0 ! set perturbation current = 1 purely real

      else  !no impedance, no foil
        if(cap.eq.0.d0) then
          if (nside.ge.1) then ! Side reactions, both 1,2,3
            b(i2div,ki2)=-1.0d0/hC
            a(i2div,ki2)=1.0d0/hC
            b(i2div,kj)=area*fc
            b(i2div,kj1)=-area*fc
            b(i2div,kj2)=-area*fc
            g(i2div)=c(ki2,j)/hC-area*fc*(c(kj,j)-c(kj1,j)-c(kj2,j)) !not O(h2)
c written only for j=1
            if(j.gt.1) g(i2div)=(c(ki2,j)-c(ki2,j-1))/hC-
     &    area*fc*(c(kj,j)-c(kj1,j)-c(kj2,j))
            if (nside.eq.1) b(i2div,kj2)=0.0d0
          else ! No side reactions
            b(i2div,ki2)=-1.0d0/hC
            a(i2div,ki2)=1.0d0/hC
            b(i2div,kj)=area*fc
            g(i2div)=c(ki2,j)/hC-area*fc*c(kj,j) ! not order h2
c     need to include mpa ! problem
         if(j.gt.1) g(i2div)=(c(ki2,j)-c(ki2,j-1))/hC-area*fc*c(kj,j) !not O(h2)
          endif
        else !nonzero value of capacitances
        if (rr .eq. 0.d0) then ! no instantaneous change across double layer
          g(i2div)=c(kp1,j)-c(kp2,j)-xt(kp1,j,kk-1+kadd)+
     &    xt(kp2,j,kk-1+kadd)
          b(i2div,kp1)=-1.d0
          b(i2div,kp2)= 1.d0
        else !nonzero capacitance and nonzero time step
          b(i2div,ki2)=-1.0d0/hC
          a(i2div,ki2)=1.0d0/hC
          b(i2div,kj)=area*fc
          b(i2div,kj1)=area*fc
          b(i2div,kj2)=area*fc
          b(i2div,kp1)=area*cap/rr*omi*2.0d0
          b(i2div,kp2)=-area*cap/rr*omi*2.0d0
          g(i2div)=(c(ki2,j)+xt(ki2,j,kk-1+kadd))/hC
          if(j.gt.1) g(i2div)=g(i2div)+(-c(ki2,j-1)-xt(ki2,j-1,
     &  kk-1+kadd))/hC
          g(i2div)=g(i2div) -area*fc*(c(kj,j)+xt(kj,j,kk-1+kadd)
     &    -c(kj1,j)-xt(kj1,j,kk-1+kadd)-c(kj2,j)-xt(kj2,j,kk-1+kadd))
     &    -area*cap*(c(kp1,j)-c(kp2,j)-xt(kp1,j,kk-1+kadd)
     &    +xt(kp2,j,kk-1+kadd))/rr*omi*2.0d0
        endif
      endif
c
      if (imp.eq.1) then ! impedance section
c     Add double-layer charging:
      b(i2div,kp2+n)=omega(ji)*area*capp
      b(i2div,kp1+n)=-omega(ji)*area*capp
      b(kj+n,kp2)=-omega(ji)*area*capp
      b(kj+n,kp1)=omega(ji)*area*capp
      g(i2div)=0.0d0
      g(kj+n)=0.0d0
      endif ! end of impedance section (Equation i2div)
      endif

c     ________________________________________________
c
  121 if (imp.eq.1) then ! impedance section
c     Copy real matrices into imag matrices:
      do 122 nct1=1, n
      do 122 nct2=1, n
      if(j.eq.1) x(nct1+kp1,nct2+kp1)=x(nct1,nct2)
      if(j.gt.1) a(nct1+kp1,nct2+kp1)=a(nct1,nct2)
      b(nct1+kp1,nct2+kp1)=b(nct1,nct2)
      if(j.eq.nj) y(nct1+kp1,nct2+6)=y(nct1,nct2)
      if(j.lt.nj) d(nct1+6,nct2+6)=d(nct1,nct2)
  122 continue
      endif
c
      if(imp.eq.0 .and. mcL.ne.1) nx=nx+1 !carry end potentials only if needed

c     double the number of variables (real & imaginary parts)
      if(imp.eq.1) n=2*n
            call band(j)
      if(imp.eq.1) n=n/2

      if(imp.eq.0 .and. mcL.ne.1) nx=nx-1

      if(j.lt.nj) go to 10

      if (imp.eq.0) then
      do 607 jj=1,nj
      do 607 i=1,nx+1
  607 c(i,jj)=xx(i,jj)+c(i,jj)  
      endif ! make corrections for imp.ne.1

c     ________________________________________________
c
c     print impedance results
c
      if (imp.eq.1) then
      vreal(ji)=c(kp1,1)-c(kp1,nj)+RG !WHT 1-11-08
      vimag(ji)=c(kp1+n,1)-c(kp1+n,nj)
c     Half-cell impedance calculations:
      vreala(ji)=c(kp1,1)-c(kp2,n1+n2)+RGn !WHT 1-11-08
      vrealc(ji)=c(kp2,n1+n2)-c(kp1,nj)+RGp !WHT 1-11-08
      vimaga(ji)=c(kp1+n,1)-c(kp2+n,n1+n2)
      vimagc(ji)=c(kp2+n,n1+n2)-c(kp1+n,nj)
c     Z magnitude and phase angle calculation:
      zmag(ji)=(vreal(ji)*vreal(ji)+vimag(ji)*vimag(ji))**0.5d0
      zmaga(ji)=(vreala(ji)*vreala(ji)+vimaga(ji)*vimaga(ji))**0.5d0
      zmagc(ji)=(vrealc(ji)*vrealc(ji)+vimagc(ji)*vimagc(ji))**0.5d0
      phi(ji)=180.0d0/pi*datan(- vimag(ji)/vreal(ji))
      phia(ji)=180.0d0/pi*datan(- vimaga(ji)/vreala(ji))
      phic(ji)=180.0d0/pi*datan(- vimagc(ji)/vrealc(ji))

      write (2,103)omega(ji),(10**4)*vreal(ji),(10**4)*vimag(ji)  !Jun29

      rr=0.01 !to compensate for rr=1 set to avoid div by 0
      ji=ji+1
      return
      endif ! end of impedance section (printing)
c
      do 56 i=1,n
   56 xp(i)=(4.0d0*c(i,2)-3.0d0*c(i,1)-c(i,3))/2.0d0/h1
      nerr=0
      nconcflag = 0
      do 25 j=1,nj
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     shoe horns:

      if(c(1,j).lt.xx(1,j)/100.d0) c(1,j)=xx(1,j)/100.d0
c     Concentration of lithium in solution in the positive is kept from going to zero by the following.
c     The minimum value can be adjusted depending to the specific scenario being tested.

      if(c(1,j).gt.xx(1,j)*100.d0) c(1,j)=xx(1,j)*100.0d0
      if (c(kp2,j).lt.(xx(kp2,j)-0.02d0)) c(kp2,j)=xx(kp2,j)-0.02d0
      if (c(kp2,j).gt.(xx(kp2,j)+0.02d0)) c(kp2,j)=xx(kp2,j)+0.02d0
      if (c(kp1,j).lt.(xx(kp1,j)-0.02d0)) c(kp1,j)=xx(kp1,j)-0.02d0
      if (c(kp1,j).gt.(xx(kp1,j)+0.02d0)) c(kp1,j)=xx(kp1,j)+0.02d0
      if (c(kp2,j).gt. 9.9d0) c(kp2,j)= 9.9d0
      if (c(kp2,j).lt.-9.9d0) c(kp2,j)=-9.9d0
      if (c(kp1,j).gt. 9.9d0) c(kp1,j)= 9.9d0
      if (c(kp1,j).lt.-9.9d0) c(kp1,j)=-9.9d0
c
      do mpa=1,npa
      if (j .ge. n1+n2) then
      if(c(2+mpa,j).lt.xx(2+mpa,j)/100.d0) then
      nerr=nerr+1
      c(2+mpa,j)=xx(2+mpa,j)/100.d0 ! use cs min
      endif
      if(c(2+mpa,j).gt.xx(2+mpa,j)*100.d0) c(2+mpa,j)=xx(2+mpa,j)*100.d0 
      if(ct3-c(2+mpa,j).le.(ct3-xx(2+mpa,j))/100.d0) then
      nerr=nerr+1
      c(2+mpa,j)=ct3-(ct3-xx(2+mpa,j))/100.d0
      endif
      if(c(2+mpa,j).ge.ct3) c(2+mpa,j)=0.999999d0*ct3
c
      else if (j .le. n1+1 .and. n1 .gt. 0) then

      if(c(2+mpa,j).lt.xx(2+mpa,j)/100.d0) nerr=nerr+1
      if(c(2+mpa,j).lt.xx(2+mpa,j)/100.d0) c(2+mpa,j)=xx(2+mpa,j)/100.d0 ! use cs min
      if(c(2+mpa,j).gt.xx(2+mpa,j)*100.d0) c(2+mpa,j)=xx(2+mpa,j)*100.d0 
      if(ct1-c(2+mpa,j).le.(ct1-xx(2+mpa,j))/100.d0) nerr=nerr+1
      if(ct1-c(2+mpa,j).le.(ct1-xx(2+mpa,j))/100.d0)
     & c(2+mpa,j)=ct1-(ct1-xx(2+mpa,j))/100.d0
      if(c(2+mpa,j).ge.ct1) c(2+mpa,j)=0.999999d0*ct1
      endif
c
c     Not shoe-horns, but trips to force smaller time steps
      if (j. ge .n1+n2 .or. (j .le. n1+1 .and. n1 .gt. 0)) then
      if(c(2+mpa,j).lt.1.0d-10) c(2+mpa,j) = 1.d-10
      if(c(2+mpa,j).lt.1.0d-3) nconcflag = 1
      endif
      enddo ! mpa
c
      if(c(1,j).lt.1.0d-10) then
      c(1,j) = 1.d-10
c      c(kj,j) = 0.0d0
      endif
      if(c(1,j).lt.1.0d-03) nconcflag = 2
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      c(nx+1,j)=c(ip1,nj)
      if(j.le.n1+1) c(nx+1,j)=c(ip1,1)
      do 25 i=1,n+1
   25 xx(i,j)=c(i,j) 

      if (jcount .gt. 3*lim .and. rr.eq.0.0d0) then
      write (2,99)

      stop
      endif
c     ________________________________________________
c
c     Decreasing time steps:
      if (jcount .gt. lim .and. rr.gt.1.0d-16) then
        tau=tau/2.0d0
        rr=tau
        Lcount=0
        if(mcLL.eq.-4) cur=curold
        ts(kk)=ts(kk-1)+tau
c     print *,'time step reduced to ', tau, ts(kk),jcount
        if (tau.lt.1.0d-4) then
      nconv=1  !trigger the flag for no convergence
          if (lpow.eq.1) then  !peak power activated
            nflag=1
            go to 66
          endif
          nt=nt-1
          pow=ed/tw/ts(nt+1)
      write (2,*) 'mass is ',tw
      write (2,*) 'energy is ',ed/tw/3.6d03
      write (2,*) 'power is ',pow
          write (2,*) kk-1,' this time step did not converge'
          if (nconcflag .eq. 1) then
            write (2,*), 'Solid concentration driven to zero.'
            write (2,*),'Consider changing vcut or the current density.'
            if (nconcflag.eq.3) then
              write (2,*) 'The solid utilization is above 1.'
            endif
          endif
          if (nconcflag .eq. 2) then
            write (2,*), 'Salt concentration driven to zero.'
            write (2,*),'Consider changing vcut or the current density.'
      endif
      write (2,*) ' '
      if (kk.gt.5) then
      write (2,*) 'dV/dt at point of failure is (V/s)', 
     &    ((xt(kp1,nj,kk-1)-xt(kp1,1,kk-1))-(xt(kp1,nj,kk-2)
     &    -xt(kp1,1,kk-2)))/(ts(kk-1)-ts(kk-2))
      write (2,*) 'Moving on to next leg if not the last'
      write (2,*) ' '
      endif
      return
          else
            iflag=0
            call calca(kk)
            go to 666
          end if
c
      else
c     ________________________________________________
c
c check for convergence
c normalize concentrations
c     calculate total amount of lithium in solution      
      totLinew=0.d0
      do j=1,nj
      fh=h3*ep3
      if(j.le.n1+n2) fh=h2*ep2
      if(j.le.n1+1) fh=h1*ep1
      if(j.eq.1) fh=0.5*h1*ep1
      if(j.eq.n1+1) fh=0.5d0*(h1*ep1+h2*ep2)
      if(j.eq.n1+n2) fh=0.5d0*(h3*ep3+h2*ep2)
      if(j.eq.nj) fh=0.5d0*h3*ep3
      totLinew=totLinew + fh*c(1,j)
      enddo
      factor=totLiold/totLinew
      do j=1,nj
      c(1,j)=factor*c(1,j)
      enddo
c     if(dabs(factor-1.d0) .ge. 1.d-4) nerr=nerr+1
c     if(jcount.ge.4) print *, totLinew,jcount
      
c
      if(nerr.ne.0) go to 8
      do jj=1,nj
      do ki=1,n
      errlim=1.d-10
      if(ki.eq.kj) errlim=1.d-16 !change to -14 if problems with convergence
      dxx=dabs( xx(ki,jj) - xx0(ki,jj) )
      enddo
      enddo
      do 55 ki=1,n
      errlim=1.d-10
      if(ki.eq.kj) errlim=1.d-16 !change to -14 if problems with convergence
      dxp=dabs( xp(ki)-xp0(ki) )
      n1hold = 1
      if (n1 .ge. 11) n1hold = n1-10
      dxx=dabs( xx(ki,n1hold) - xx0(ki,n1hold) )
      dxx2=dabs( xx(ki,n1+n2+10)-xx0(ki,n1+n2+10) )
      dxx3=dabs( xx(ki,2)-xx0(ki,2) )
      if(dxx .gt. 1.d-9*dabs(xx(ki,n1hold)).and.dxx.gt.errlim) then
      go to 8
      endif
      if(dxx2.gt.1.d-9*dabs(xx(ki,n1+n2+10)).and.dxx2.gt.errlim) then
      go to 8
      endif
      if(dxx3.gt.1.d-9*dabs(xx(ki,2)).and.dxx3.gt.errlim) then
      go to 8
      endif
   55 continue


      if(mcLL.eq.-4) then ! maximum power
      vvold=vv
      vv=xx(kp1,nj)-xx(kp1,1)
      pwr=cur*vv


      if(Lcount.eq.1) then
        curlow=0.d0
        curhigh=0.d0
        curmaxx=cur
        pwrmaxx=cur*vv
      else

      if(cur*vv.gt.pwrmaxx) then
        if(curmaxx.lt.cur) then
          curlow=curmaxx
        else
          if(Lcount.gt.2) curhigh=curmaxx
        endif
        curmaxx=cur
        pwrmaxx=cur*vv

      else ! no new power maximum
        if(cur.lt.curmaxx) then
          curlow=cur
        elseif(cur.lt.curhigh) then
          curhigh=cur
        endif
      endif

c     May have a problem with Rint if the difference between
c     cur and curold becomes too small.  However, the bounds
c     from curlow and curhigh should still give valid answers.
      Rint=-(vv-vvold)/(cur-curold)+RG ! internal resistance
      OCP=vv+(Rint-RG)*cur
      curmax=OCP/Rint/2.d0 ! estimate of current at max power
      pwrmax=OCP**2/Rint/4.d0 ! estimate of max power
      if(Lcount.eq.2) curhigh=2.d0*curmax
      endif
      curold=cur

      cur=cur-1.d0
      if(Lcount.gt.1) then
      if(curmax.gt.1.2d0*cur) curmax=1.2d0*cur
      cur=curmax
      if(cur.lt.curlow .or. cur.gt.curhigh) cur=0.5d0*(curlow+curhigh)
      endif
      if(Lcount.lt.18 .and. dabs(cur-curold).gt.1.d-5*dabs(cur)) go to 3
      if(Lcount.lt.5) go to 3
      cur=curold
      endif

      if ((xx(1,nj)/xt(1,nj,kk-1)).lt.0.75d0.and.xx(1,nj).lt.100.0d0)
     1rr=rr/2.0d0

c     ________________________________________________
c
c     if(lpow.ne.1) write (2,*) jcount,' iterations required'
c
      do 60 ll=1, nj  ! save present time results in xt()
      do 60 lk=1,n
   60 xt(lk,ll,kk)=xx(lk,ll)

      do 61 ll=1,nj
      do 61 lk=1,nnj
   61 cssold(ll,lk)=css(ll,lk)


c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(rr.ne.0.0d0) then
      do 58 j=1,nj  ! nonzero time step
      term_s2(j)=termn_s2(j)
      term_s1(j)=termn_s1(j)
   58 term(j)=termn(j)
      else ! zero time step
      do 65 j=1,nj
      term(j)=0.d0
      term_s1(j)=0.0d0
      term_s2(j)=0.0d0

      fac=1.d0
      if(j.eq.n1+2 .and. n1 .gt. 0)
     &fac=((ep2+epp2)*exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2+1) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)

      hn=h2 
      if(n1.gt.0 .and. j.le.n1+1) hn=h1
      if(n1.gt.0 .and. j.le.n1+1) area=area1
      if(j.gt.n1+n2) hn=h3
      if(j.gt.n1+n2) area=area3

      if(j.gt.1) term(j)=
     &-(df(j)+fac*df(j-1))*(c(1,j)-c(1,j-1))/hn/2.d0
     &-(1.d0-0.5d0*(tm(j)+tm(j-1)))*c(ki2,j-1)/fc
      if (j.gt.1) term_s1(j)=
     &-(do2(j)+fac*do2(j-1))*(c(kS1,j)-c(kS1,j-1))/hn/2.0d0
     &+hn*area*vol*(0.375d0*c(kj1,j)+0.125d0*c(kj1,j-1))
      if (j.gt.1) term_s2(j)=
     &-(dh2(j)+fac*dh2(j-1))*(c(kS2,j)-c(kS2,j-1))/hn/2.0d0
     &+hn*area*vol*(0.375d0*(c(kj3,j)+c(kj2,j))+0.125d0*
     &(c(kj3,j-1)+c(kj2,j-1)))

      fac=1.d0
      if(j.eq.n1+1 .and. n1 .gt. 0)
     &fac=((ep2+epp2)**exbrug2)/((ep1+epp1)**exbrug1)
      if(j.eq.n1+n2) fac=((ep3+epp3)**exbrug3)/((ep2+epp2)**exbrug2)

      hn=h2
      if(n1.gt.0 .and. j.lt.n1+1) hn=h1
      if(n1.gt.0 .and. j.lt.n1+1) area=area1
      if(j.ge.n1+n2) hn=h3
      if(j.ge.n1+n2) area=area3

      if(j.lt.nj) term(j)=term(j)
     &-(fac*df(j)+df(j+1))*(c(1,j)-c(1,j+1))/hn/2.d0
     &+(1.d0-0.5d0*(tm(j)+tm(j+1)))*c(ki2,j)/fc
      if (j.lt.nj) term_s2(j)=term_s2(j)
     &-(fac*dh2(j)+dh2(j+1))*(c(kS2,j)-c(kS2,j+1))/hn/2.0d0
     &+hn*area*vol*(0.375d0*(c(kj3,j)+c(kj2,j))+0.125d0*
     &(c(kj3,j+1)+c(kj2,j+1)))
   65 if (j.lt.nj) term_s1(j)=term_s1(j)
     &-(fac*do2(j)+do2(j+1))*(c(kS1,j)-c(kS1,j+1))/hn/2.0d0
     &+hn*area*vol*(0.375d0*c(kj1,j)+0.125d0*c(kj1,j+1))
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      end if
c
   66 continue

      return
      end
c
c***********************************************************************
      subroutine calca(kk)
c     subroutine to calculate diffusion in solid particles
c     of active material, based on a superposition integral.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/gas/ epg1,epg2,epg3
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      dimension ar(4,maxt),bz(6)
c
      do mpa=1,npa ! different particle sizes
      do 319 l=1,nt
      ai2(l,mpa)=0.0d0
  319 ai(l,mpa)=0.0d0
c
      do 70 i=1,kk-1
      ar(1,i)=dfs3*(ts(kk)-ts(i))/Rad3pa(mpa)/Rad3pa(mpa)
      ar(2,i)=dfs3*(ts(kk)-ts(i+1))/Rad3pa(mpa)/Rad3pa(mpa)
      ar(3,i)=dfs1*(ts(kk)-ts(i))/Rad1pa(mpa)/Rad1pa(mpa)
      ar(4,i)=dfs1*(ts(kk)-ts(i+1))/Rad1pa(mpa)/Rad1pa(mpa)
c
      do 69 m=1,2
      t1=ar(m,i)
      t2=ar((m+2),i)
c
      a1=0.0d0
      a12=0.0d0
c
      s=1.644934066848d0
c
c     Bessel's function zeros:
      bz(1)=2.4048255577d0
      bz(2)=5.5200781103d0
      bz(3)=8.6537281103d0
      bz(4)=11.7915344391d0
      bz(5)=14.9309177086d0
c
      if (shape3.gt.2.0d0) then
c     spherical particles:
      if (t1 .gt. 0.06d0) then
c
      do 59 j=1,5
      y1=dble(j*j)*pi*pi*t1
   59 if (y1 .le. 1.5d02) a1=a1+(expg(-y1))/dble(j*j)
      a1=2.0d0*(s-a1)/pi/pi
c
      else
c
      if (t1.LE.0.0d0) then
      a1=0.0d0
      else
      do 60 j=1,3
      z=dble(j)/dsqrt(t1)
      call erfcg(z,e)
      y2=dble(j*j)/t1
      if(y2 .ge. 1.5d02) then
      da=-dble(j)*dsqrt(pi/t1)*e
      else
      da=expg(-y2)-dble(j)*dsqrt(pi/t1)*e
      end if
   60 a1=a1+da
      a1=-t1 + 2.0d0*dsqrt(t1/pi)*(1.0d0+2.0d0*a1)
      end if
c
      end if
      else
c
      if (shape3.lt.2.0d0) then
c     planar particles:
      if(t1 .gt. 0.06d0) then
c
      do 61 j=1,5
      da=((-1.0d0)**j)*(1.0d0 - expg(-dble(2*j+1)**2
     1*pi*pi*t1))/dble(2*j+1)**2
   61 a1=a1+da
      a1=4.0d0*a1/pi/pi
c
      else
c
      do 62 j=1,3
      z=dble(j)/2.0d0/dsqrt(t1)
      call erfcg(z,e)
      da=((-1.0d0)**j)*(expg(-dble(j*j)/4.0d0/t1)
     &-dble(j)/2.0d0*dsqrt(pi/t1)*e)
   62 a1=a1+da
      a1=2.0d0*dsqrt(t1/pi)*(1.0d0+2.0d0*a1)
c
      end if
      else ! shape = 2
c     cylindrical particles:
      if (t1.gt.0.06d0) then
c
      do 63 j=1,5
      da=(1.0d0-expg(-bz(j)*bz(j)*t1))/bz(j)/bz(j)
   63 a1=a1+da
      a1=2.0d0*a1
c
      else
c
      a1=2.0d0*dsqrt(t1/pi)-t1/4.0d0-5.0d0*(t1**1.5d0)/96.0d0/dsqrt(pi)
     1-31.0d0*t1*t1/2048.0d0
c
      end if
      end if
      end if
c
      if (n1 .eq. 0) go to 36
c     (skip calculations of Li diffusion in the solid
c     if the negative electrode is metal foil)
c
      if (shape1.gt.2.0d0) then
c     spherical particles:
      if(t2 .gt. 0.06d0) then
c
      do 64 j=1,5
      y2=dble(j*j)*pi*pi*t2
   64 if(y2 .le. 1.5d02) a12=a12+(expg(-y2))/dble(j*j)
      a12=2.0d0*(s-a12)/pi/pi
c
      else
c
      if (t2.eq.0.0d0) then
      a12=0.0d0
      else
      do 65 j=1,3
      z=dble(j)/dsqrt(t2)
      call erfcg(z,e)
      y2=dble(j*j)/t2
      if(y2 .gt. 1.5d02) then
      da=-dble(j)*dsqrt(pi/t2)*e
      else
      da=expg(-y2)-dble(j)*dsqrt(pi/t2)*e
      end if
   65 a12=a12+da
      a12=-t2 + 2.0d0*dsqrt(t2/pi)*(1.0d0+2.0d0*a12)
      end if
      end if
c
      else
      if (shape1.lt.2.0d0) then
c     planar particles:
      if(t2 .gt. 0.06d0) then
c
      do 66 j=1,5
      da=((-1.0d0)**j)*(1.0d0 - expg(-dble(2*j+1)**2
     1*pi*pi*t2))/dble(2*j+1)**2
   66 a12=a12+da
      a12=4.0d0*a12/pi/pi
c
      else
c
      do 67 j=1,3
      z=dble(j)/2.0d0/dsqrt(t2)
      call erfcg(z,e)
      da=((-1.0d0)**j)*(expg(-dble(j*j)/4.0d0/t2)
     &-dble(j)/2.0d0*dsqrt(pi/t2)*e)
   67 a12=a12+da
      a12=2.0d0*dsqrt(t2/pi)*(1.0d0+2.0d0*a12)
c
      end if
      else
c     cylindrical particles:
      if (t2.gt.0.06d0) then
c
      do 68 j=1,5
      da=(1.0d0-expg(-bz(j)*bz(j)*t2))/bz(j)/bz(j)
   68 a12=a12+da
      a12=2.0d0*a12
c
      else
c
      a12=2.0d0*dsqrt(t2/pi)-t2/4.0d0-5.0d0*(t2**1.5d0)/96.0d0
     &/dsqrt(pi)-31.0d0*t2*t2/2048.0d0
c
      end if
      end if
      end if
c
   36 continue
c
      ar(m,i)=a1
   69 ar((m+2),i)=a12
c
      ai(kk-i,mpa)=ar(1,i)-ar(2,i)
   70 ai2(kk-i,mpa)=ar(3,i)-ar(4,i)
c
      do j=1,nj
      sumpa(mpa,j)=0.d0
      enddo
      do j=1,n1+1
      if (kk.gt.2) then  
      do  54 i=1, kk-2
   54 if(ts(i+1)-ts(i).ne.0.0d0) sumpa(mpa,j)=sumpa(mpa,j)
     &+(xt(2+mpa,j,i+1)- xt(2+mpa,j,i))*ai2(kk-i,mpa)/(ts(i+1)-ts(i))
      endif
      enddo
      do j=n1+n2,nj
      if (kk.gt.2) then  
      do 95 i=1, kk-2
   95 if(ts(i+1)-ts(i).ne.0.0d0) sumpa(mpa,j)=sumpa(mpa,j) ! problem
     &+(xt(2+mpa,j,i+1)- xt(2+mpa,j,i))*ai(kk-i,mpa)/(ts(i+1)-ts(i))
      endif
      enddo
      enddo ! mpa
c
      return
      end
c***********************************************************************
      subroutine erfcg(z,e)
c   an error function complement subroutine.
      implicit real*8(a-h,o-z)
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
c
      a1=0.254829592d0
      a2=-0.284496736d0
      a3=1.421413741d0
      a4=-1.453152027d0
      a5=1.061405429d0
      if(z .lt. 2.747192d0) then
      t3=1.0d0/(1.0d0+0.3275911d0*z)
      e=(a1*t3+a2*t3*t3+a3*(t3**3)+a4*(t3**4)
     1+a5*(t3**5))*expg(-z*z)
      else
c
      if(z .gt. 25.0d0) then
      e=0.0d0
      else
c
      sum=0.0d0
      max=z*z + 0.5d0
      fac=-0.5d0/z/z
      sum=fac
      tl=fac
      n=1
   10 n=n+1
      if(n .gt. max) go to 15
      tn=tl*dble(2*n-1)*fac
      sum=sum + tn
      if(tn .lt. 1.0d-06) go to 15
      tl=tn
      go to 10
   15 e=(expg(-z*z))*(1.0d0+sum)/dsqrt(pi)/z
      end if
      end if
c
      return
      end
c***********************************************************************
      subroutine band(j)
c     subroutine to solve coupled, linear differnece equations
      implicit real*8(a-h,o-z)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/mat/ b,d
      common/bnd/ a,c,g,x,y
      dimension b(17,17),d(17,35)
      dimension a(17,17),c(17,221),g(17),x(17,17),y(17,17)
      dimension e(17,18,221)
      save e, np1
  101 format (15h determ=0 at j=,i4)
      n=nx
      if (j-2)  1,6,8
    1 np1 = n + 1
      do 2 i=1,n
      d(i,2*n+1)= g(i)
      do 2 l=1,n
      lpn= l + n
    2 d(i,lpn)= x(i,l)
      call matinv(n,2*n+1,determ)
      if (determ)  4,3,4
    3 print 101, j
    4 do 5 k=1,n
      e(k,np1,1)= d(k,2*n+1)
      do 5 l=1,n
      e(k,l,1)= - d(k,l)
      lpn= l + n
    5 x(k,l)= - d(k,lpn)
      return
    6 do 7 i=1,n
      do 7 k=1,n
      do 7 l=1,n
    7 d(i,k)= d(i,k) + a(i,l)*x(l,k)
    8 if (j-nj)  11,9,9
    9 do 10 i=1,n
      do 10 l=1,n
      g(i)= g(i) - y(i,l)*e(l,np1,j-2)
      do 10 m=1,n
   10 a(i,l)= a(i,l) + y(i,m)*e(m,l,j-2)
   11 do 12 i=1,n
      d(i,np1)= - g(i)
      do 12 l=1,n
      d(i,np1)= d(i,np1) + a(i,l)*e(l,np1,j-1)
      do 12 k=1,n
   12 b(i,k)= b(i,k) + a(i,l)*e(l,k,j-1)
      call matinv(n,np1,determ)
      if (determ)  14,13,14
   13 print 101,j
   14 do 15 k=1,n
      do 15 m=1,np1
   15 e(k,m,j)= - d(k,m)
      if (j-nj)  20,16,16
   16 do 17 k=1,n
   17 c(k,j)= e(k,np1,j)
      do 18 jj=2,nj
      m= nj - jj + 1
      do 18 k=1,n
      c(k,m)= e(k,np1,m)
      do 18 l=1,n
   18 c(k,m)= c(k,m) + e(k,l,m)*c(l,m+1)
      do 19 l=1,n
      do 19 k=1,n
   19 c(k,1)= c(k,1) + x(k,l)*c(l,3)
   20 return
      end
c***********************************************************************
      subroutine matinv(n,m,determ)
c     matrix inversion program used by BAND(j).
      implicit real*8(a-h,o-z)
      common/mat/ b,d
      dimension b(17,17),d(17,35)
      dimension id(17)
      determ=1.0d0
      do 1 i=1,n
    1 id(i)=0
      do 18 nn=1,n
      bmax=1.1d0
      do 6 i=1,n
      if(id(i).ne.0) go to 6
      bnext=0.0d0
      btry=0.0d0
      do 5 j=1,n
      if(id(j).ne.0) go to 5
      if(dabs(b(i,j)).le.bnext) go to 5
      bnext=dabs(b(i,j))
      if(bnext.le.btry) go to 5
      bnext=btry
      btry=dabs(b(i,j))
      jc=j
    5 continue
      if(bnext.ge.bmax*btry) go to 6
      bmax=bnext/btry
      irow=i
      jcol=jc
    6 continue
      if(id(jc).eq.0) go to 8
      determ=0.0d0
      return
    8 id(jcol)=1
      if(jcol.eq.irow) go to 12
      do 10 j=1,n
      save=b(irow,j)
      b(irow,j)=b(jcol,j)
   10 b(jcol,j)=save
      do 11 k=1,m
      save=d(irow,k)
      d(irow,k)=d(jcol,k)
   11 d(jcol,k)=save
   12 f=1.0d0/b(jcol,jcol)
      do 13 j=1,n
   13 b(jcol,j)=b(jcol,j)*f
      do 14 k=1,m
   14 d(jcol,k)=d(jcol,k)*f
      do 18 i=1,n
      if(i.eq.jcol) go to 18
      f=b(i,jcol)
      do 16 j=1,n
   16 b(i,j)=b(i,j)-f*b(jcol,j)
      do 17 k=1,m
   17 d(i,k)=d(i,k)-f*d(jcol,k)
   18 continue
      return
      end
c***********************************************************************
      subroutine nucamb(kk,il2)
c     subroutine to print detailed profiles as go along
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/tprop/df(221),cd(221),tm(221),ddo2(221),ddh2(221),
     1ddf(221),dcd(221),dtm(221),dfu(221),d2fu(221),do2(221),dh2(221)
      common /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      dimension zz(221)
  109 format(f6.1,', ',f11.2,', ',f8.4,', ',g10.3,', ',f8.3,', ',g10.4
     1,', ',g10.4,', ',g10.4,', ',g10.4,', ',g10.4)
  310 format('Distance   C Elec    C Sol Surf  Liq Pot    Solid Pot  ',
     &'Liq Cur   j main     j side 1    j side 2    j side 3')
  311 format('(microns)    (mol/m3)   x or y     (V)         (V) ',
     1'    (A/m2)     (A/m2)       (A/m2)      (A/m2)      (A/m2)')
   44 format(' t = ',1pe18.6,' min','  temp = ',0pf11.5,' deg C')
c
      do 5 i=1,n1+1
      w=dble(i-1)-1.d-8 ! helps rounding when diff'ing files from diff comps
    5 zz(i) = w*h1*1.0d06
      do 71 i=n1+2,n2+n1
      w=dble(i-(n1+1))-1.d-8
   71 zz(i)=zz(n1+1)+w*h2*1.0d06
      do 72 i=n1+n2+1,nj
      w=dble(i-(n1+n2))-1.d-8
   72 zz(i)=zz(n2+n1)+w*h3*1.0d06
c
      write (4,*) ' '
      write (4,310) 
      write (4,311)

      write (4,44) ts(kk)/60.0d0,t-273.15d0
c            print 44, ts(kk)/60.0d0,t-273.15d0
c
      do 10 j=1,nj,il2
      if (j .le. n1+1) then
      csol=ct1
      else
      csol=ct3
      end if
      if(j.le.n1+1) then
      curden_main=area1*fc*xx(kj,j)
      curden_side=area1*fc*xx(kj3,j)
      else if(j.ge.n1+n2) then
      curden_main=area3*fc*xx(kj,j) !corrected mixed-up area variables 1/15/99
      curden_side=area3*fc*xx(kj3,j)
      else
      curden_main=0.0d0
      curden_side=0.0d0
      endif
      mpa=1 ! only one size for detailed profiles

   10 write (4,109) zz(j),xx(1,j),xx(2+mpa,j)/csol,xx(kp2,j),xx(kp1,j),
     &xx(ki2,j),xx(kj,j)*fc,xx(kj1,j)*fc,xx(kj2,j)*fc,xx(kj3,j)*fc
      if (kk.eq.1) write (4,*) 'that was the zero-time solution!'
      write (4,*) ' '
      write (4,*) ' '

      return
      end
c***********************************************************************
      subroutine peak(n,lim,curr, vcut)
c     subroutine to calculate the peak-power capability at
c     a particular time in the discharge.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
  311 format(f8.5,', ',f7.3,',  ',f8.3,', ',f8.5)
c
c     Peak power current ramp section:
      write (2,*) ' '
      write (2,*) '    PEAK POWER '
      write (2,*) ' '
      write (2,*) 'avg. cell pot, current, power, min pot,'
      write (2,*) '  (V)  ,  (A/m2),  (W/m2),     (V)'
c
c     Duration of current pulse is 30 seconds.
c
      nend = 1
      curpmax = 0.0d0
      curtry = 0.0d0
      pwrmax = 0.0d0
      vvpmax = 0.0d0
      curmin=0.d0
      pwrpmax=0.d0
      rrmax=30.0d0
      cur=curr
      curmax=0.d0
  127 kcount=0
      fact=10.0d0
      vfmax=0.d0
      k=nt+2
      do 126 j=1,nj
      do 126 i=1,n
  126 xt(i,j,k)=xt(i,j,k-1)
      ppow=0.0d0
      ii=0
c     Ramp current:
c
  128 continue
      if (ii.gt.100) return
c
      energ=0.0d0
  130 ii=ii+1
      cur=cur+fact  ! start a new current density
         if(curmax.ne.0.d0) then
         cur=0.5d0*(curpmax+curmax)
         if(mod(ii,2).eq.1) cur=0.5d0*(curpmax+curmin)
         if(cur.eq.0.d0) cur=0.5d0*(curmax+curmin)
         if(vfmax.gt.0.d0 .and. vfmax.lt.vcut) then
         curtry=curpmax+(vcut-vvpmax)/(vfmax-vvpmax)*(curmax-curpmax)
     &+0.01d0*(curmax-curpmax)*dble(mod(ii,3)-1)
         if(curtry.lt.curmax .and. curtry.gt.curpmax) cur=curtry
         endif
      if(vfmax.gt.vcut .and. curmin.gt.0.0d0) then
      v2=pwrmax/curmax
      cur2=curmax
      v1=pwrmin/curmin
      cur1=curmin
      vm=pwrpmax/curpmax
      curm=curpmax
      resis=-((pwrpmax-pwrmin)/(curpmax-curmin)
     &-(pwrpmax-pwrmax)/(curpmax-curmax))/(curmin-curmax)
      Uop=(pwrpmax-pwrmin)/(curpmax-curmin)+resis*(curpmax+curmin)
      curtry=Uop/2.d0/resis+0.1d0*(curmax-curmin)*dble(mod(ii,3)-1)
      if(curtry.lt.curmax .and. curtry.gt.curmin) cur=curtry
      endif
         endif
      kkflag=0
      iflag=0
      nflag=0
      k=nt+2
      timpk=0.0d0
      rr=0.0d0
      ts(k)=ts(k-1)
      call comp(n,lim,k,rr,kkflag,nflag,1,jcount)
      call cellpot(k,vv,0,1,lflag)
      vlast=vv
      rr=0.02d0
  129 kkflag=kkflag+1
      k=k+1
      ts(k)=ts(k-1)+rr
      call calca(k)
c
      call comp(n,lim,k,rr,kkflag,nflag,1,jcount)
c
      if (nflag.eq.1.and.kcount.lt.20) then
      write (2,*) 'Peak current decreased',kcount,fact
      if(cur.lt.curpmax) then
      write (2,*) 'Convergence on power failed; already converged at a
     & higher current'
      write (2,*) 'Best results obtained are:'
      write (2,311) pwrpmax/curpmax,curpmax,pwrpmax,vvpmax
      return
      endif
      curmax=cur
      vfmax=0.d0
      go to 128
      endif
      if (kcount.ge.10) return
c
      call cellpot(k,vv,0,1,lflag)
      energ=energ+(vlast+vv)*(ts(k)-ts(k-1))*cur/2.0d0
c
      timpk=timpk+rr
      if (dabs(timpk-30.0d0).gt.0.1d0) then
c
      if (timpk.lt.30.0d0) then
      vlast=vv
c     Increasing time steps:
      if(jcount.lt.12 .and. kkflag.gt.5 .and. (2.0d0*rr
     1+timpk).lt.30.0d0 .and. iflag.eq.0) then
         rr=rr*2.0d0
      print *,'should not be here'
        print *, 'next time step increased to ', rr,'(s)'
      end if
      if(timpk+rr.gt.30.0d0) iflag=1
      if(timpk+rr.gt.30.0d0) rr=30.0d0-timpk
      go to 129
      end if
c
      end if
      ppow=energ/30.0d0
      write (2,311) ppow/cur,cur,ppow,vv
         if(ppow.gt.pwrpmax .and. vv.gt.vcut) then
         if(curpmax.lt.cur) then
            curmin=curpmax
            pwrmin=pwrpmax
            vfmin=vvpmax
         else
            curmax=curpmax
            pwrmax=pwrpmax
            vfmax=vvpmax
         endif
         curpmax=cur
         pwrpmax=ppow
         vvpmax=vv
         endif
         if(vv.lt.vcut) then
         if(curmax.eq.0.d0 .or. cur.lt.curmax) then
            curmax=cur
            pwrmax=ppow
            vfmax=vv
         endif
         else
         if(cur.gt.curmin .and. cur.lt.curpmax) then
            curmin=cur
            pwrmin=ppow
            vfmin=vv
         endif
         if(cur.gt.curpmax) then
            curmax=cur
            pwrmax=ppow
            vfmax=vv
         endif
         endif
      if(curmax.eq.0.0d0) go to 128
      if(curmin. lt. 0.999d0*curmax) go to 128
      write (2,311) pwrpmax/curpmax,curpmax,pwrpmax,vvpmax
      if(nend .eq. 0) return
      nend = 0
      go to 127
c
      end
c***********************************************************************
      subroutine cellpot(kk,v,li,lpow,lflag)
c     subroutine to calculate the cell potential at the end
c     of a time step and to print results.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/power/ ed,Vold,ranodesave,rcathdesave,heat,qlosstot
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
      common/side/rksc1,c1init,c2init,rksa1,term_s1(221),vol,
     &rksa2,rksa3,rksc3,rksc2,UsO2,UsH2,cn2,term_s2(221),nside
      common/RG/ RG,RGn,RGp,RGext
      save ut1, ut3
  306 format (f8.5,',   ',f6.2,',   ',g12.5,',   ',g12.5,',   '
     &,g12.5,',   ',g12.5)
  307 format(f9.3,', ',f6.4,', ',f6.4,', ',f8.5,', ',f8.5,', ',f8.2,
     &', ', f6.2, ', ', f9.2, ', ', f6.2, ', ', f6.2, ', ', g10.3)
  309 format(f6.3,', ',f8.5,', ',f5.2,', ',f12.3,',',f9.5,', ',f6.2,
     & ',',g10.3, ',', g10.3)

      if (kk.le.1) then
      ut3=utz(1,nj)
      if (n1.ne.0) ut1=utz(1,1)
      endif

c Define the individual pore-wall fluxes.  These are necessary for
c multiple particle sizes in order to treat properly a film resistance,
c and for a rigorous energy balance

      do j=1,nj
c     This makes the fluxes the same for different particle sizes
c     when a new leg of the driving profile is encountered.
      if (rr.eq.0.0d0) then ! time step check
      do mpa=1,npa
      if (j.le.n1+1) then
      pwfindiv(mpa,j)=xx(kj,j)*area1pa(mpa)/area1
      elseif (j.ge.n1+n2) then
      pwfindiv(mpa,j)=xx(kj,j)*area3pa(mpa)/area3
      else
      pwfindiv(mpa,j)=0.0d0
      pwfindiv_d(mpa,j)=0.0d0
      endif
      enddo !mpa

      else !rr>0

        if (li.eq.1) then
      kadd=0
      do mpa=1,npa
      if (j.le.n1+1) then
      pwfindiv(mpa,j)=(ai2(1,mpa)*(xt(2+mpa,j,kk-1+kadd)-xx(2+mpa,j))
     &/rr-sumpa(mpa,j))*Rad1pa(mpa)
      if (mvdc1.eq.1) pwfindiv(mpa,j)=xx(kj,j)
      pwfindiv_d(mpa,j)=-Rad1pa(mpa)*ai2(1,mpa)/rr
      utz(mpa,j)=utz(mpa,j)-pwfindiv(mpa,j)*rr*area1pa(mpa)
     &/(1.0d0-ep1-epp1-epf1-epg1)/ct1/vf1(mpa)
      elseif (j.ge.n1+n2) then
      pwfindiv(mpa,j)=(ai(1,mpa)*(xt(2+mpa,j,kk-1+kadd)-xx(2+mpa,j))
     &/rr-sumpa(mpa,j))*Rad3pa(mpa)
      if (mvdc3.eq.1) pwfindiv(mpa,j)=xx(kj,j)
      pwfindiv_d(mpa,j)=-Rad3pa(mpa)*ai(1,mpa)/rr
      utz(mpa,j)=utz(mpa,j)-pwfindiv(mpa,j)*rr*area3pa(mpa)/
     &(1.0d0-ep3-epp3-epf3-epg3)/ct3/vf3(mpa)
      else
      pwfindiv(mpa,j)=0.0d0
      pwfindiv_d(mpa,j)=0.0d0
      endif
      enddo ! mpa
        endif ! li
      endif ! time step check
      enddo ! j

c     Material and current balance tests:
c Current balance - integrate the current balance equation
      sum=0.0d0 !mat balance on electrolyte
      sum1=0.0d0 !current neg elec, main
      sum1_s1=0.0d0 !side current neg elec, side reaction 1
      sum3=0.0d0 !main current pos elec
      sum3_s1=0.0d0 !side current pos elec, side reaction 1
      sum1_s2=0.0d0 !side current neg elec, side reaction 2
      sum3_s2=0.0d0 !side current pos elec, side reaction 2
      sum1_s3=0.0d0 !side current neg elec, side reaction 3
      sum3_s3=0.0d0 !side current pos elec, side reaction 3
      sum4=0.0d0 !mat balance oxygen 
      sum5=0.0d0 !mat balance hydrogen

      if (n1 .gt. 2) then
      do 85 j=2,n1
      sum1_s1=sum1_s1+area1*fc*h1*xx(kj1,j)
      sum1_s2=sum1_s2+area1*fc*h1*xx(kj2,j)
      sum1_s3=sum1_s3+area1*h1*xx(kj3,j)
      sum1=sum1+area1*fc*h1*xx(kj,j) !j main
      sum4=sum4+xx(kS1,j)*h1*epg1 !conc O2
      sum5=sum5+xx(kS2,j)*h1*epg1 !conc H2
   85 sum=sum+xx(1,j)*(ep1+epp1)*h1 !conc elec
      endif
      sum=sum+(xx(1,1)+xx(1,n1+1))*(ep1+epp1)*h1/2.0d0
      sum1=sum1+area1*fc*h1*(xx(kj,1)+xx(kj,n1+1))/2.0d0
      sum1_s1=sum1_s1+area1*fc*h1*(xx(kj1,1)+xx(kj1,n1+1))
     & /2.0d0
      sum1_s2=sum1_s2+area1*fc*h1*(xx(kj2,1)+xx(kj2,n1+1))
     & /2.0d0
      sum1_s3=sum1_s3+area1*h1*(xx(kj3,1)+xx(kj3,n1+1))/2.0d0
      sum4=sum4+h1*(xx(kS1,1)+xx(kS1,n1+1))/2.0d0*epg1
      sum5=sum5+h1*(xx(kS2,1)+xx(kS2,n1+1))/2.0d0*epg1

      do 86 j=n1+2,n2+n1-1
      sum4=sum4+h2*xx(kS1,j)*epg2
      sum5=sum5+h2*xx(kS2,j)*epg2
   86 sum=sum+xx(1,j)*(ep2+epp2)*h2
      sum=sum+(xx(1,n1+1)+xx(1,n2+n1))*(ep2+epp2)*h2/2.0d0
      sum4=sum4+h2*(xx(kS1,n1+1)+xx(kS1,n1+n2))/2.0d0*epg2
      sum5=sum5+h2*(xx(kS2,n1+1)+xx(kS2,n1+n2))/2.0d0*epg2

      do 87 j=n2+n1+1,nj-1
      sum3=sum3+area3*fc*xx(kj,j)*h3
      sum3_s1=sum3_s1+area3*fc*xx(kj1,j)*h3
      sum3_s2=sum3_s2+area3*fc*xx(kj2,j)*h3
      sum3_s3=sum3_s3+area3*fc*xx(kj3,j)*h3
      sum4=sum4+h3*xx(kS1,j)*epg3
      sum5=sum5+h3*xx(kS2,j)*epg3
   87 sum=sum+xx(1,j)*(ep3+epp3)*h3
      sum3=sum3+area3*fc*(xx(kj,n1+n2)+xx(kj,nj))
     & *h3/2.0d0
      sum3_s1=sum3_s1+area3*fc*(xx(kj1,n1+n2)+xx(kj1,nj))
     & *h3/2.0d0
      sum3_s2=sum3_s2+area3*fc*(xx(kj2,n1+n2)+xx(kj2,nj))
     & *h3/2.0d0
      sum3_s3=sum3_s3+area3*fc*(xx(kj3,n1+n2)+xx(kj3,nj))
     & *h3/2.0d0
      sum=sum+(xx(1,n1+n2)+xx(1,nj))*h3*(ep3+epp3)/2.0d0
      sum4=sum4+h3*(xx(kS1,n1+n2)+xx(kS1,nj))/2.0d0*epg3
      sum5=sum5+h3*(xx(kS2,n1+n2)+xx(kS2,nj))/2.0d0*epg3
      sum4=sum4/thk !adjusts to average concentration
      sum5=sum5/thk
c     calculate total salt in cell from initial profile:
      w=xt(1,n1+2,1)*(dble(n2-1)*(ep2+epp2)*h2+dble(n1)*(ep1+epp1)*h1
     1+dble(n3)*(ep3+epp3)*h3)
      if(lflag.eq.0) w=w-(xt(1,n1+2,1)-xt(1,1,1))*(dble(n1)*
     1(ep1+epp1)*h1+dble(n3)*(ep3+epp3)*h3)
      if(lflag.eq.0) w=w-(xt(1,n1+2,1)-xt(1,1,1))*(ep2+epp2)*h2
c     material-balance parameter should be ca=1.00
      ca=sum/w

c     Calculate cell potential from dif of solid-phase potentials:
      v=xt(kp1,nj,kk)-xt(kp1,1,kk)-RG*cur !WHT 1-11-08


c Calculate the resistances in the cell, R=V/I
      mpa=1 ! only one size for resistance calculation
      xxsave=xx(2+mpa,n1+1)
      call ekin(n1+1,kk,0,0.0d0,mpa,mvdc1,mvdc3)
      xx(2+mpa,n1+1)=xxsave
      Ua_end=g0
      xxsave=xx(2+mpa,n1+n2)
      call ekin(n1+n2,kk,0,0.0d0,mpa,mvdc1,mvdc3)
      xx(2+mpa,n1+n2)=xxsave
      Uc_end=g0
      R_anode=(xt(kp1,1,kk)-xt(kp2,n1+1,kk)-Ua_end)/cur+RGn !WHT 1-11-08
      R_sep=(xt(kp2,n1+2,kk)-xt(kp2,n1+n2,kk))/cur
      R_cathode=RGp-(xt(kp1,nj,kk)-xt(kp2,n1+n2,kk)-Uc_end)/cur !WHT 1-11-08
      R_total=R_anode+R_sep+R_cathode
c
      if(rr.eq.0.d0 .and. kk.eq.1) then
      call ekin(1,kk,1,ut1,mpa,mvdc1,mvdc3)
      Ua=g0
      xx(2+mpa,1)=xxsave
      xxsave=xx(2+mpa,n1+n2)
      call ekin(n1+n2,kk,1,ut3,mpa,mvdc1,mvdc3)
      Uc=g0
      xx(2+mpa,n1+n2)=xxsave
      Uoc = Uc-Ua 
      call temperature(kk,v,ut3,ut1,Uoc)
      tprint=t-273.15d0
	endif

c     Calculate utilization of two electrodes based on coulombs passed:
      if(li.eq.1) then

c     Calculate energy density by running sum of currentxvoltage:
        ed=ed+((Vold-RG*xx(ki2,n1+1))*xx(ki2,n1+1)
     &+(V-RG*cur)*cur)*rr/2.0d0
        Vold=v

c Here the current is always positive, the total current through cell
c Calculate the new utilizations from integrating the intercalation
c current throughout the cell

        ut3=ut3-sum3*rr/fc/(1.0d0-ep3-epf3-epp3-epg3)/dble(n3)/h3/ct3
        if (n1 .gt. 0)   !Fix how utilization is calculated for a foil anode
c Multiply sum1_s3 by fc for the case of absorption
     &  ut1=ut1-(sum1+sum1_s3*fc)*rr/fc/(1.0d0-ep1-epf1-epp1-epg1)/
     &dble(n1)/h1/(ct1)

      th=ts(kk)/6.0d01


c
      xxsave=xx(2+mpa,1)
      call ekin(1,kk,1,ut1,mpa,mvdc1,mvdc3)
      Ua=g0
      xx(2+mpa,1)=xxsave
      xxsave=xx(2+mpa,n1+n2)
      call ekin(n1+n2,kk,1,ut3,mpa,mvdc1,mvdc3)
      Uc=g0
      xx(2+mpa,n1+n2)=xxsave
      Uoc = Uc-Ua 

      call temperature(kk,v,ut3,ut1,Uoc)
        tprint=t-273.15d0


c assumed temperature dependence for solid-state diffusion
      dfs1=dfs1save*dexp((EbarS1)*(t-298.15d0)/(t*298.15d0))
      dfs3=dfs3save*dexp((EbarS3)*(t-298.15d0)/(t*298.15d0))

      pressn2=28.0d0*t*8.206d-5 !calculate press in atm, 32.33 for N2
      presso2=sum4*t*8.206d-5 !oxygen pressure
      pressh2=sum5*t*8.206d-5 !hydrogen pressure

      if(lpow.ne.0) then
c     ! isothermal peak-power output:
        write (2,309) v,ca,cur,v*cur
      else ! not peak power
        if (nside.ne.0.and.nneg.eq.7) then !NiMH with side reaction
          write (2,307) th,ut1,ut3,v,Uoc,
     1    cur,pressh2*0.9869d0,presso2*0.9869d0,
     1    (pressh2+pressn2+presso2)*0.9869d0,tprint,qq
        else
          write (2,307) th,ut1,ut3,v,Uoc,
     1    cur,tprint,qq
          write (8,"(e15.9)") tprint,dfs1,dfs3
        endif
      endif

      if (cur.ne.0.0d0) write (6,306) th,cur,R_anode,R_sep,
     &R_cathode,R_total
      endif ! li

      jref=(n1+1+n1+n2)/2 ! halfcells
      write (3,310) th,xt(kp1,1,kk)-xt(kp2,jref,kk)+cur*RG/2.d0
     %,xt(kp1,nj,kk)-xt(kp2,jref,kk)-cur*RG/2.d0,cur,t-273.15d0,qq/thk
  310 format (f11.6,2f12.6,f12.4,f14.6,1pe16.7)

      print 410, th*60.d0,cur,t-273.15d0
c     print 410, th*60.d0,cur,ts(kk)-273.15d0
          write (2,307) th,ut1,ut3,v,Uoc,
     1    cur,tprint,qq
      
  410 format (f9.3,', ',f12.5,', ',f9.5,', ',f6.3,',',f9.3,', ',f6.2,
     & ',',g10.3,',',f8.5)
      if (li.eq.1) then
c     Calculate total heat generated as running sum
      heat=heat+qq*rr/3600.d0 !Wh/m2
      qlosstot=qlosstot+qloss*rr/3600.d0 !Wh/m2
      endif
c
      return
      end


c***********************************************************************
      subroutine sol(nmax,jj)
      implicit real*8(a-h,o-z)
c     This subroutine calculates the solid-phase concentration profiles,
c     at one position (jj) and at one time.
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      dimension cs(50)
  100   format(f8.4,'  ',f8.4)
c
c     set initial value of solid concentration
      mpa=1 ! only one size for solid profiles
      do 88 i=1, 50
   88 cs(i)=xt(2+mpa,jj,1)
c
      if (jj .le. n1) then
      dfs = dfs1
      Rad = Rad1
      else if (jj .ge. n1+n2) then
      dfs = dfs3
      Rad = Rad3
      else
      print *, 'jsol selected for node in separator'
      return
      endif
c
c     complete calculations for 50 points along radius of particle
      nmax=nmax-1  ! added
      do 10 i=1,50
      y2=0.02d0*dble(i)
c
      sum1=0.0d0
      do 20 kk=1,nmax
      k=nmax+1-kk
c
      t1=(ts(nmax+1)-ts(k))*dfs/Rad/Rad
      sum2=sum1
c
c     calculate c bar (r,t1)
      sum1=0.0d0
      r1=1.0d0
c
      do 89 j=1,15
      r1=-r1
      y1=dble(j*j)*pi*pi*t1
      y3=dble(j)*pi*y2
      if (y1 .gt. 1.50d02) then
      da=0.0d0
      else
      da=expg(-y1)
      end if
   89 sum1=sum1-2.0d0*r1*da*dsin(y3)/dble(j)/pi/y2
      sum1=1.0d0-sum1
c
c     perform superposition
c
      cs(i)=cs(i)+(xt(2+mpa,jj,k+1)+xt(2+mpa,jj,k)-2.0d0*xt(2+mpa,jj,1)
     1)*(sum1-sum2)/2.0d0
   20 continue
c
   10 continue
      nmax=nmax+1  ! added
c
      write (7,*) ' '
      write (7,*) 'time is ',ts(nmax)/60.0d0,' minutes'
      write (7,*) ' '
      do 90 i=1, 50, 1
   90 write (7,100) .02d0*dble(i),cs(i)/ct1
c
      return
      end
c***********************************************************************
      subroutine mass(re,rs3,rs1,rf,rpl,rc,rcn,rcp)
c     calculate mass (kg/m2) of the cell based on 
c     densities and volume fractions of components.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
c
c     mass of positive electrode
      c1=h3*dble(n3)*(re*ep3+rpl*epp3+rs3*(1.0d0-ep3-epf3-epp3-epg3)+
     &rf*epf3)
c     mass of separator
      s=(re*ep2+rpl*epp2+rc*(1.0d0-ep2-epp2-epg2))*h2*dble(n2-1)
c
c     mass of negative electrode
      if (n1.gt.0) a1=h1*dble(n1)*(re*ep1+rpl*epp1+rs1*(1.0d0-ep1-
     1epf1-epp1-epg1)+rf*epf1)
      if (n1.eq.0) a1=h1*rs1
c
c     mass of current collectors
      cc1=rcn*hcn+rcp*hcp
c
      tw=c1+s+a1+cc1+residm
c
      return
      end
c***********************************************************************
      subroutine temperature(kk,v,ut,ut2,Uoc)
c     subroutine to recompute the cell temperature based 
c     on heat generation, heat capacity, and heat losses.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
      common/side/rksc1,c1init,c2init,rksa1,term_s1(221),vol,
     &rksa2,rksa3,rksc3,rksc2,UsO2,UsH2,cn2,term_s2(221),nside
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      common/RG/ RG,RGn,RGp,RGext
c
c     The entropy and open-circuit potential for each electrode
c     should be given in ekin with respect to a Li reference electrode.
c     Heat generation is positive if exothermic.
c
c     The energy balance is now done according to the method given
c     in Rao and Newman, J. Electrochem. Soc., 144 (1997), 2697-2704.
c     Cp(dT/dt)-Q=-Integral(sum over reactions (a*in,i*Uh,i))dx-IV

c     The heat-transfer coefficient is for heat transferred out of
c     one side of the cell; it is defined based on cell area.
c     htcc is a per-cell heat-transfer coefficient.
c
      htcc=htc/dble(Ncell)

      sum_main1=0.0d0
      sum_main3=0.0d0
      sum_ox1=0.0d0
      sum_ox3=0.0d0
      sum_hye1=0.0d0
      sum_hyr1=0.0d0

      dudtO2=-0.00168d0
      dudtO2=0.0d0
      dudtH2=0.0d0
      dHread=100.0d0 !units here J/mole of H2 readsorbed. all reversible.

      do j=2,n1
      do mpa=1,npa
      utzs=utz(mpa,j)
      xxsave=xx(2+mpa,j)
      call ekin(j,kk,1,utzs,mpa,mvdc1,mvdc3)
      xx(2+mpa,j)=xxsave
      Umain1=g0
      dUdTmain1=dudt
      sum_main1=sum_main1+area1pa(mpa)*fc*pwfindiv(mpa,j)*h1
     &*(Umain1-t*dUdTmain1)
      if (nside.ge.1 .and. nprop.eq.13) then
      sum_ox1=sum_ox1+area1*fc*xx(kj1,j)*h1*(UsO2+r*t/fc*dlog(xx(kS1,j)*
     &8.206d-5*t)-t*dudtO2)
      sum_hye1=sum_hye1+area1*fc*xx(kj3,j)*h1*(UsH2-r*t/fc/2.0d0*
     &dlog(xx(kS2,j)*8.206d-5*t)-t*dudtH2)
      sum_hyr1=sum_hyr1+area1*xx(kj2,j)*h1*dHread
      endif
      enddo !mpa
      enddo

      do mpa=1,npa
      utzs=utz(mpa,1)
      xxsave=xx(2+mpa,1)
      call ekin(1,kk,1,utzs,mpa,mvdc1,mvdc3)
      xx(2+mpa,1)=xxsave
      Umain1=g0
      dUdTmain1=dudt
      sum_main1=sum_main1+area1pa(mpa)*fc*pwfindiv(mpa,1)*h1*
     &(Umain1-t*dUdTmain1)/2.d0
      if (nside.ge.1 .and. nprop.eq.13) then
      sum_ox1=sum_ox1+area1*fc*xx(kj1,1)*h1*(UsO2+r*t/fc*dlog(xx(kS1,1)*
     &8.206d-5*t)-t*dudtO2)
      sum_hye1=sum_hye1+area1*fc*xx(kj3,1)*h1*(UsH2-r*t/fc/2.0d0*
     &dlog(xx(kS2,1)*8.206d-5*t)-t*dudtH2)
      sum_hyr1=sum_hyr1+area1*xx(kj2,1)*h1*dHread
      endif
      enddo !mpa

      do mpa=1,npa
      utzs=utz(mpa,n1+1)
      xxsave=xx(2+mpa,n1+1)
      call ekin(n1+1,kk,1,utzs,mpa,mvdc1,mvdc3)
      xx(2+mpa,n1+1)=xxsave
      Umain1=g0
      dUdTmain1=dudt
      sum_main1=sum_main1+area1pa(mpa)*fc*pwfindiv(mpa,n1+1)*h1*
     &(Umain1-t*dUdTmain1)/2.0d0
      if (nside.ge.1 .and. nprop.eq.13) then
      sum_ox1=sum_ox1+area1*fc*xx(kj1,n1+1)*h1*(UsO2+r*t/fc*
     &dlog(xx(kS1,n1+1)*8.206d-5*t)-t*dudtO2)
      sum_hye1=sum_hye1+area1*fc*xx(kj3,n1+1)*h1*(UsH2-r*t/fc/2.0d0*
     &dlog(xx(kS2,n1+1)*8.206d-5*t)-t*dudtH2)
      sum_hyr1=sum_hyr1+area1*xx(kj2,n1+1)*h1*dHread
      endif
      enddo !mpa

      do j=n1+n2+1,nj-1
      do mpa=1,npa
      utzs=utz(mpa,j)
      xxsave=xx(2+mpa,j)
      call ekin(j,kk,1,utzs,mpa,mvdc1,mvdc3)
      xx(2+mpa,j)=xxsave
      Umain3=g0
      dUdTmain3=dudt
      sum_main3=sum_main3+area3pa(mpa)*fc*pwfindiv(mpa,j)
     &*h3*(Umain3-t*dUdTmain3)
      if (nside.ge.1 .and. nprop.eq.13) then
      sum_ox3=sum_ox3+area3*fc*xx(kj1,j)*h1*(UsO2+r*t/fc*
     &dlog(xx(kS1,j)*8.206d-5*t)-t*dudtO2)
      endif
      enddo !mpa
      enddo

      do mpa=1,npa
      utzs=utz(mpa,n1+n2)
      xxsave=xx(2+mpa,n1+n2)
      call ekin(n1+n2,kk,1,utzs,mpa,mvdc1,mvdc3)
      xx(2+mpa,n1+n2)=xxsave
      Umain3=g0
      dUdTmain3=dudt
      sum_main3=sum_main3+area3pa(mpa)*fc*pwfindiv(mpa,n1+n2)*h3*
     &(Umain3-t*dUdTmain3)/2.0d0
      if (nside.ge.1 .and. nprop.eq.13) then
      sum_ox3=sum_ox3+area3*fc*xx(kj1,n1+n2)*h1*(UsO2+r*t/fc*
     &dlog(xx(kS1,n1+n2)*8.206d-5*t)-t*dudtO2)
      endif
      enddo !mpa

      do mpa=1,npa
      utzs=utz(mpa,nj)
      xxsave=xx(2+mpa,nj)
      call ekin(nj,kk,1,utzs,mpa,mvdc1,mvdc3)
      xx(2+mpa,nj)=xxsave
      Umain3=g0
      dUdTmain3=dudt
      sum_main3=sum_main3+area3pa(mpa)*fc*pwfindiv(mpa,nj)*h3*
     &(Umain3-t*dUdTmain3)/2.0d0
      if (nside.ge.1 .and. nprop.eq.13) then
      sum_ox3=sum_ox3+area3*fc*xx(kj1,nj)*h1*(UsO2+r*t/fc*
     &dlog(xx(kS1,nj)*8.206d-5*t)-t*dudtO2)
      endif
      enddo !mpa

      sumtotal=sum_main1+sum_main3+sum_ox1+sum_ox3+sum_hye1+sum_hyr1

      if (kk.eq.1) then
      qq=cur*(Uoc-v)-cur**2*RGext !WHT 1-11-08
      else
      qq=-cur*v-sumtotal-cur**2*RGext !WHT 1-11-08
      endif

c     This energy balance includes a residual mass.  Input near top.
      if(lht.eq.0)
     &t=t+(rr/(Cp*(dens*thk+residm)))*(htcc*(tam-t)+qq)
      
		if(rr.eq.0.d0) then
c	print *, t, ' iterations on zero time step'
c	stop
	endif

c This is the original energy balance.
c     t=t+(rr/(dens*Cp*thk))*(htcc*(tam-t)+cur*(Uoc-v-t*Soc)) ! fixed 6-30-01

      if (lht.eq.1) then
      print *,'lht.eq.1 is NOT WORKING'
c
c     Calculate htc instead of temperature: the heat-transfer coefficient
c     required to keep the temperature constant is
c     calculated as a function of time.  The heat-transfer coef.
c     is calculated for heat transferred out of one side of the
c     cell stack.  Htcc is defined as a per-cell heat-transfer
c     coefficient.
c
         if (t.ne.tam) then
c           htc=dble(Ncell)*cur*(Uoc-v-t*Soc)/(t-tam) ! fixed 6-30-01
         else
           htc=0.0d0
         endif
         htcc=htc/dble(Ncell)
      endif
c
      return
      end
c***********************************************************************
      double precision function expg(x)
      implicit real*8 (a-h,o-z)
      expg=0.d0
      if(x.gt.-700.d0) expg=dexp(x)
      return
      end

c***********************************************************************
      double precision function cosh(x)
      implicit real*8 (a-h,o-z)
      cosh=1.d0
      cosh=(expg(x)+expg(-x))/2.d0
      return
      end

c***********************************************************************
      subroutine ekin(j,kk,lag,utz,mpa,mvdc1,mvdc3)
      implicit real*8(a-h,o-z)
c     This subroutine evaluates the Butler-Volmer equations.
c     It also provides a library of data for various positive
c     and negative active materials.
      parameter(maxt=900)
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/power/ ed,Vold,ranodesave,rcathdesave,heat,qlosstot
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
      common/mat/ b,d
      common/bnd/ a,c,g,x,y
      common/side/rksc1,c1init,c2init,rksa1,term_s1(221),vol,
     &rksa2,rksa3,rksc3,rksc2,UsO2,UsH2,cn2,term_s2(221),nside
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/RG/ RG,RGn,RGp,RGext
      dimension b(17,17),d(17,35)
      dimension a(17,17),c(17,221),g(17),x(17,17),y(17,17)
c
c     Calculate average open-circuit potential in either
c     electrode if lag=1, otherwise lag=0

c Section for Paxton kinetic data
c NOTE: always set rate constants rka1, rka3 to 1.0 if using NiMH

      if (nprop.eq.13) then
      RTF=r*t/fc ! change danger
      ATCN=0.65d0
      CTCN=0.35d0
      ATCP=0.35d0
      CTCP=0.65d0
      H2OMW=18.016D0
      PKOHMW=56.11D0
      Pcon=xx(1,j)/1.0d6
      Pcs=xx(2+mpa,j)/1.0d6
      CNGMAX=ct1/1.0d6
      CPSMAX=ct3/1.0d6

      EXCDN= 7.85D-04 /(0.012644D0**CTCN*0.046814D0**ATCN*
     1        (0.5D0*CNGMAX)**(CTCN+ATCN))
      EXCDP= 1.04D-04 /(0.012644D0**CTCP*0.046814D0**ATCP*
     1         (0.5D0*CPSMAX)**(CTCP+ATCP))     

      DN = 1.001D0 + 47.57D0*Pcon - 776.22D0*Pcon**2
      DN1D = 47.57D0 - 1552.44D0*Pcon
      DN2D = -1552.44D0

      AC = 1.004D0 - 36.23D0*Pcon**0.5 + 1374.3D0*Pcon
     1 - 17850.7*Pcon**1.5 + 55406.0D0*Pcon**2
     1 +7.16856D05*Pcon**2.5
      AC1D=-18.115D0*Pcon**(-0.5)+1374.3D0-2.6776D04*Pcon**0.5
     1 +1.10812D05*Pcon + 1.7921D06*Pcon**1.5
      AC2D = 9.0575D0/Pcon**1.5 - 13388.0D0/Pcon**0.5
     1 + 1.10812D05 + 2.6882D06*Pcon**0.5
      WAC = 1.0002D0 - 21.238D0*Pcon - 4.1312D03*Pcon**2.0D0
      WAC1D= -21.238D0 -2.0d0*4.1312D03*Pcon    
      endif


c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     OPEN-CIRCUIT POTENTIAL FUNCTIONS:
c
c     g0 is the open-circuit potential in terms of the solid
c     concentration, xx(2+mpa,j), with respect to a lithium metal
c     electrode
c     g1 is the derivative of the open-circuit potential wrt
c     the solid concentration
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     FOR THE NEGATIVE ELECTRODE
c
      if (j .le. n1+1) then
c
      if (lag.eq.1) xx(2+mpa,j)=utz*ct1
c

c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
      go to (51,52,53,54,55,56,57,58),nneg
      write (2,*) 'Please enter data for negative electrode
     & #8 in subroutine ekin'
      stop
   51 go to 111  ! Li foil
   52 go to 112  ! Carbon (petroleum coke)
   53 go to 113  ! MCMB 2528 graphite (Bellcore)
   54 go to 114  ! TiS2
   55 go to 115  ! Tungsten oxide (LixWO3 with 0<x<0.67)
   56 go to 116  ! Lonza KS6 graphite (Bellcore)
   57 go to 117  ! Metal Hydride (from Paxton)
   58 go to 118  ! add your own

c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Li foil (works only if n1 = 0)
c
  111 g0 = 0.0d0
      g1 = 0.0d0
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      dudt = 0.0d0
      go to 97
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Carbon (petroleum coke):
c
  112 c1=-0.132056d0
      c2=1.40854d0
      c3=-3.52312d0
      g0=c1+c2*expg(c3*xx(2+mpa,j)/ct1)
      g1=c2*c3*expg(c3*xx(2+mpa,j)/ct1)/ct1
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      dudt = 0.0d0
      go to 97
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
  113 continue
c   MCMB 2510 carbon (Bellcore)
c      c1=-0.160d0
c      c2=1.32d0
c      c3=-3.0d0
c      g0=c1+c2*expg(c3*xx(2+mpa,j)/ct1)
c      g0=g0+10.d0*expg(-2000.d0*xx(2+mpa,j)/ct1)
c      g1=c2*c3*expg(c3*xx(2+mpa,j)/ct1)/ct1
c      g1=g1-10.d0*2000.d0/ct1*expg(-2000.d0*xx(2+mpa,j)/ct1)
c     MCMB 2528 graphite measured by Chris Bogatu 2000, Telcordia and PolyStor materials
c     for 0.01 < x < 0.9
      sto = xx(2+mpa,j)/ct1
      g0 = 0.194d0+1.5d0*expg(-120.0d0*sto)
     & +0.0351d0*dtanh((sto-0.286d0)/0.083d0)
     & - 0.0045d0*dtanh((sto-0.849d0)/0.119d0)
     & - 0.035d0*dtanh((sto-0.9233d0)/0.05d0)
     & - 0.0147d0*dtanh((sto-0.5d0)/0.034d0)
     & - 0.102d0*dtanh((sto-0.194d0)/0.142d0)
     & - 0.022d0*dtanh((sto-0.9d0)/0.0164d0)
     & - 0.011d0*dtanh((sto-0.124d0)/0.0226d0)
     & + 0.0155d0*dtanh((sto-0.105d0)/0.029d0)
      g1 = -1.5d0*(120.0d0/ct1)*expg(-120.0d0*sto)
     &+(0.0351d0/(0.083d0*ct1))*((dcosh((sto-0.286d0)/0.083d0))**(-2))
     &-(0.0045d0/(ct1*0.119d0))*((dcosh((sto-0.849d0)/0.119d0))**(-2))
     &-(0.035d0/(ct1*0.05d0))*((dcosh((sto-0.9233d0)/0.05d0))**(-2))
     &-(0.0147d0/(ct1*0.034d0))*((dcosh((sto-0.5d0)/0.034d0))**(-2))
     &-(0.102d0/(ct1*0.142d0))*((dcosh((sto-0.194d0)/0.142d0))**(-2))
     &-(0.022d0/(ct1*0.0164d0))*((dcosh((sto-0.9d0)/0.0164d0))**(-2))
     &-(0.011d0/(ct1*0.0226d0))*((dcosh((sto-0.124d0)/0.0226d0))**(-2))
     &+(0.0155d0/(ct1*0.029d0))*((dcosh((sto-0.105d0)/0.029d0))**(-2))
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      dudt = 1.0d-3*(0.00527 + 3.29927*sto - 91.79326*sto**2
     & + 1004.91101*sto**3 - 5812.27813*sto**4 + 19329.75490*sto**5
     & - 37147.89470*sto**6 + 38379.18127*sto**7 - 16515.05308*sto**8)
     & /(1 - 48.09287*sto + 1017.23480*sto**2 - 10481.80419*sto**3
     & + 59431.30001*sto**4 - 195881.64880*sto**5 + 374577.31520*sto**6
     & - 385821.16070*sto**7 + 165705.85970*sto**8)

      if(nneg.ne.0) go to 97 ! bypass new form
c     MCMB 2528 graphite measured by Chris Bogatu 2000, 
c     Telcordia and PolyStor materials.
c     Modified May 2003 to match data from Joongpyo Shim
c     for 0.01 < x < 0.99
      sto = xx(2+mpa,j)/ct1
      g0 = 0.124d0+1.5d0*expg(-150.0d0*sto)
     & +0.0351d0*dtanh((sto-0.286d0)/0.083d0)
     & - 0.0045d0*dtanh((sto-0.9d0)/0.119d0)
     & - 0.035d0*dtanh((sto-0.99d0)/0.05d0)
     & - 0.0147d0*dtanh((sto-0.5d0)/0.034d0)
     & - 0.102d0*dtanh((sto-0.194d0)/0.142d0)
     & - 0.022d0*dtanh((sto-0.98d0)/0.0164d0)
     & - 0.011d0*dtanh((sto-0.124d0)/0.0226d0)
     & + 0.0155d0*dtanh((sto-0.105d0)/0.029d0)
      g1 = -1.5d0*(150.0d0/ct1)*expg(-150.0d0*sto)
     &+(0.0351d0/(0.083d0*ct1))*((dcosh((sto-0.286d0)/0.083d0))**(-2))
     &-(0.0045d0/(ct1*0.119d0))*((dcosh((sto-0.9d0)/0.119d0))**(-2))
     &-(0.035d0/(ct1*0.05d0))*((dcosh((sto-0.99d0)/0.05d0))**(-2))
     &-(0.0147d0/(ct1*0.034d0))*((dcosh((sto-0.5d0)/0.034d0))**(-2))
     &-(0.102d0/(ct1*0.142d0))*((dcosh((sto-0.194d0)/0.142d0))**(-2))
     &-(0.022d0/(ct1*0.0164d0))*((dcosh((sto-0.98d0)/0.0164d0))**(-2))
     &-(0.011d0/(ct1*0.0226d0))*((dcosh((sto-0.124d0)/0.0226d0))**(-2))
     &+(0.0155d0/(ct1*0.029d0))*((dcosh((sto-0.105d0)/0.029d0))**(-2))
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      dudt = 0.0d0
      go to 97
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     TiS2
c
  114 delt=-5.58d-04
      zeta=8.1d0
c     ct1=2.9d04
      g0=2.17d0+(dlog((ct1-xx(2+mpa,j))/xx(2+mpa,j))+delt*xx(2+mpa,j)
     &+zeta)/frt
      g1=(delt-ct1/xx(2+mpa,j)/(ct1-xx(2+mpa,j)))/frt
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      go to 97
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Tungsten oxide (LixWO3 with 0<x<0.67):
c     literature data from Whittingham
c
  115 c1=2.8767d0
      c2=-0.9046d0
      c3=0.76679d0
      c4=-0.15975d0
      c5=0.671d0
      g0=c1+c2*xx(2+mpa,j)/ct1+c3*xx(2+mpa,j)*xx(2+mpa,j)/ct1/ct1+
     1c4*expg(100.0d0*(xx(2+mpa,j)/ct1-0.671d0))
      g1=c2/ct1+2.0d0*c3*xx(2+mpa,j)/ct1/ct1+
     1c4*100.0d0*expg(100.0d0*(xx(2+mpa,j)/ct1-0.671d0))/ct1
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      go to 97
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Bellcore graphite (Lonza KS6)
c
  116 c1=0.7222d0
      c2=0.13868d0
      c3=0.028952d0
      c4=0.017189d0
      c5=0.0019144d0
      c6=0.28082d0
      c7=0.79844d0
      c8=0.44649d0
      xtem=xx(2+mpa,j)/ct1
      g0=c1+c2*xtem+c3*(xtem**0.5d0)-c4*(xtem**(-1.0d0))
     1+c5*(xtem**(-1.5d0))+c6*expg(15.0d0*(0.06d0-xtem))-c7*
     1expg(c8*(xtem-0.92d0))
      g1=c2/ct1+c3*0.5d0*(ct1**(-0.5d0))*(xx(2+mpa,j)**(-0.5d0))+
     1c4*ct1*(xx(2+mpa,j)**(-2.0d0))-c5*1.5d0*(ct1**1.5d0)*(xx(2+mpa,j)
     &**(-2.5d0))-c6*15.0d0/ct1*expg(15.0d0*(0.06d0-xtem))
     1-c7*c8/ct1*expg(c8*(xtem-0.92d0))
c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      go to 97
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  117 continue
c     Metal hydride (based on Blaine Paxton's masters thesis)

      CNORM=xx(2+mpa,j)/ct1
      c1=40.0d0
      c2=50.0d0
      c3=5000.0d0
      c4=0.8d0
      c5=0.2d0
      c6=-0.91d0
      g0=(expg(c2*(c5-cnorm))-expg(c1*(cnorm-c4)))/c3+c6
      g1=-1.0d0/c3*(c2*expg(c2*(c5-cnorm))+c1*expg(c1*(cnorm-c4)))
      g1=g1/ct1


      dudt=0.0d0 
c kinetic expressions for MH here

      if (lag.eq.1) go to 99

      RTF=r*t/fc
      ATCN=0.58d0
      CTCN=0.35d0
      ATCP=0.36d0
      CTCP=0.58d0
      H2OMW=18.016D0
      PKOHMW=56.11D0
      Pcon=xx(1,j)/1.0d6
      Pcs=xx(2+mpa,j)/1.0d6
      CNGMAX=ct1/1.0d6
      CPSMAX=ct3/1.0d6

      EXCDN= 7.85D-04 /(0.012644D0**CTCN*0.046814D0**ATCN*
     1         (0.5D0*CNGMAX)**(CTCN+ATCN))
      EXCDP= 1.04D-04 /(0.012644D0**CTCP*0.046814D0**ATCP*
     1         (0.5D0*CPSMAX)**(CTCP+ATCP))

      DN = 1.001D0 + 47.57D0*Pcon - 776.22D0*Pcon**2
      DN1D = 47.57D0 - 1552.44D0*Pcon
      DN2D = -1552.44D0

        AC = 1.004D0 - 36.23D0*Pcon**0.5 + 1374.3D0*Pcon
     1 - 17850.7*Pcon**1.5 + 55406.0D0*Pcon**2
     1 +7.16856D05*Pcon**2.5
        AC1D=-18.115D0*Pcon**(-0.5)+1374.3D0-2.6776D04*Pcon**0.5
     1 +1.10812D05*Pcon + 1.7921D06*Pcon**1.5
      AC2D = 9.0575D0/Pcon**1.5 - 13388.0D0/Pcon**0.5
     1 + 1.10812D05 + 2.6882D06*Pcon**0.5
      WAC = 1.0002D0 - 21.238D0*Pcon - 4.1312D03*Pcon**2.0D0
      WAC1D= -21.238D0 -2.0d0*4.1312D03*Pcon  

c assumed temperature dependence for exchange current
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      if(cngmax.le.pcs) print *, 'large H ',cngmax,pcs

      U0=rka1*EXCDN/H2OMW**ATCN
      U1=DN-Pcon*PKOHMW
      U2=DN1D-PKOHMW
      U3= (AC*Pcon)**CTCN*(WAC*U1)**ATCN
      U4= PCS**CTCN*(CNGMAX-Pcs)**ATCN
      U6= (AC*Pcon)**CTCN*ATCN*(WAC*U1)**(ATCN-1.0D0)*
     1(WAC*U2+WAC1D*U1)+(WAC*U1)**ATCN*CTCN*
     1(AC*PCON)**(CTCN-1.0D0)*(AC+Pcon*AC1D)
      U7= CTCN*(CNGMAX-Pcs)**ATCN*PCS**(CTCN-1.0D0)-
     1ATCN*Pcs**CTCN*(CNGMAX-Pcs)**(ATCN-1.0D0)

      h0=U0*U3*U4/fc*1.0d4
c     print *,'negative exchange current', h0*fc
c     pause
      h1=U0*U3*U7/fc/1.0d2
      h2=U0*U6*U4/fc/1.0d2
      r1a=ATCN*frt
      r1c=CTCN*frt
      r2a=r1a*(xx(kp1,j)-xx(kp2,j)-g0)
      r2c=r1c*(xx(kp1,j)-xx(kp2,j)-g0)

      de=expg(-r2c)-expg(r2a)
      pe=expg(-r2c)+expg(r2a)

      b(2+mpa,1)=h2*de
      b(2+mpa,kp2)=h0*(r1c*expg(-r2c)+r1a*expg(r2a))
      b(2+mpa,kp1)=-b(2+mpa,kp2)
      b(2+mpa,2+mpa)=h1*de+h0*(r1c*g1*expg(-r2c)+r1a*g1*expg(r2a))
      b(2+mpa,kj)=1.0d0
      g(2+mpa)=-h0*de-xx(kj,j)

      go to 99 ! go to 99 not 97 to avoid other kinetics!
c have completely finished here

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c 118 continue

c     write (2,*) 'Please enter data for negative electrode
c    & #8 in subroutine ekin'
c     stop
  118 c1=-0.160d0
      c2=1.32d0
      c3=-3.0d0
      g0=c1+c2*expg(c3*xx(2+mpa,j)/ct1)
      g0=g0+10.d0*expg(-2000.d0*xx(2+mpa,j)/ct1)
      g1=c2*c3*expg(c3*xx(2+mpa,j)/ct1)/ct1
      g1=g1-10.d0*2000.d0/ct1*expg(-2000.d0*xx(2+mpa,j)/ct1)
      ti0n=dexp((Ebarka)*(t-298.15d0)/(t*298.15d0))
      rka1=rka1save*ti0n
      dudt = 0.0d0
      go to 97

c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     KINETIC EXPRESSIONS FOR NEGATIVE ELECTRODE
c
c     h0 is the exchange current density (A/m2)
c     h1 is the derivative of io wrt solid concentration, xx(2+mpa,j)
c     h2 is the derivative of io wrt electrolyte concen., xx(1,j)
c
c Use full Paxton kinetics
c
c     NONAQUEOUS LIQUIDS
c
   97 if (lag.eq.1) go to 99
      alpha=0.5d0
      alphc=0.5d0
      if (n1 .eq. 0) then
      h0=rka1*dsqrt(xx(1,j))
      h1 = 0.0d0
      h2 = rka1/dsqrt(xx(1,j))/2.0d0
      else
      h0=rka1*dsqrt(xx(1,j))*dsqrt(ct1-xx(2+mpa,j))*dsqrt(xx(2+mpa,j))
      h1=rka1*dsqrt(xx(1,j))*dsqrt(ct1-xx(2+mpa,j))*dsqrt(xx(2+mpa,j))
     &*ct1/(ct1-xx(2+mpa,j))/xx(2+mpa,j)/2.0d0
      h2=rka1*dsqrt(ct1-xx(2+mpa,j))*dsqrt(xx(2+mpa,j))
     &/dsqrt(xx(1,j))/2.0d0
      endif
 
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     POLYMER
c
c     cmax=3.92d03
c     alpha=0.5d0
c     alphc=0.5d0
c     h0=rka1*dsqrt(xx(1,j))*dsqrt(cmax-xx(1,j))*dsqrt(ct1-xx(2+mpa,j))
c    1*dsqrt(xx(2+mpa,j))
c     h1=-rka1*dsqrt(xx(1,j))*dsqrt(cmax-xx(1,j))*dsqrt(ct1-xx(2+mpa,j))
c    1*dsqrt(xx(2+mpa,j))*(1.0d0/(ct1-xx(2+mpa,j))-1.0d0/xx(2+mpa,j))/2.0d0
c     h2=-rka1*dsqrt(xx(2+mpa,j))*dsqrt(ct1-xx(2+mpa,j))*dsqrt(cmax-xx(1,j))
c    1*dsqrt(xx(1,j))*(1.0d0/(cmax-xx(1,j))-1.0d0/xx(1,j))/2.0d0
c

      end if
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     FOR THE POSITIVE ELECTRODE
c
      if (j .ge. n1+n2) then
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
      if (lag.eq.1) xx(2+mpa,j)=utz*ct3
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13),npos
      write (2,*) 'Please enter data for positive electrode
     & #13 in subroutine ekin'
      stop
    1 go to 201  ! TiS2
    2 go to 202  ! Spinel Mn2O4 (lower plateau)
    3 go to 203  ! NaCoO2:  Sodium cobalt oxide
    4 go to 204  ! Spinel Mn2O4 (upper plateau)
    5 go to 205  ! Tungsten oxide (LixWO3 with 0<x<0.67)
    6 go to 206  ! CoO2 (Cobalt dioxide)
    7 go to 207  ! V2O5 (Vanadium oxide)
    8 go to 208  ! NiO2 (Nickel dioxide)
    9 go to 209  ! Spinel Mn2O4 (Bellcore)
   10 go to 210  ! V6O13 (Vanadium oxide)
   11 go to 211  ! LiAl0.2Mn1.8O4F0.2 spinel (Bellcore)
   12 go to 212  ! NiOOHHy from Albertus and Newman
   13 go to 213  ! add your own
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     TiS2
c
  201 delt=-5.58d-04
      zeta=8.1d0
      g0=2.17d0+(dlog((ct3-xx(2+mpa,j))/xx(2+mpa,j))+delt*xx(2+mpa,j)
     &+zeta)/frt
      g1=(delt-ct3/xx(2+mpa,j)/(ct3-xx(2+mpa,j)))/frt
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Spinel Mn2O4 (lower plateau)
c
  202 c1=2.06307d0
      c2=-0.869705d0
      c3=8.65375d0
      c4=0.981258d0
      a1=c3*(xx(2+mpa,j)/ct3-c4)
      g0=c1+c2*(dtanh(a1))
      g1=c2*c3/ct3/(dcosh(a1))/(dcosh(a1))
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     NaCoO2:  Sodium Cobalt Oxide (P2 phase, 0.3<y<0.92)
c
  203 c1=4.4108d0
      c2=-2.086d0
      c3=0.10465d0
      c4=133.42d0
      c5=89.825d0
      c6=0.16284d0
      c7=145.01d0
      c8=71.92d0
      c9=0.01d0
      c10=200.0d0
      c11=0.3d0
      c12=0.885d0
      a1=xx(2+mpa,j)/ct3
      g0=c1+c2*a1+c3*dtanh(-c4*a1+c5)+c6*dtanh(-c7*a1+c8)+c9*
     1expg(-c10*(a1-c11))-c9*expg(c10*(a1-c12))
      g1=c2/ct3-c3*c4/dcosh(-c4*a1+c5)/dcosh(-c4*a1+c5)/ct3-c6*c7
     1/dcosh(-c7*a1+c8)/dcosh(-c7*a1+c8)/ct3
     &-c9*c10*expg(-c10*(a1-c11))/ct3
     1-c9*c10*expg(c10*(a1-c12))/ct3
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Spinel Mn2O4 (upper plateau 0.2<y<1.0)
c     Literature version
c
  204 a1=4.06279d0
      a2=0.0677504d0
      a3=21.8502d0
      a4=12.8268d0
      a5=0.105734d0
      a6=1.00167d0
      a7=-0.379571d0
      a8=1.575994d0
      a9=0.045d0
      a10=71.69d0
      a11=0.01d0
      a12=200.0d0
      a13=0.19d0
      g0=a1+a2*dtanh(-a3*xx(2+mpa,j)/ct3+a4)-a5*((a6-xx(2+mpa,j)
     &/ct3)**a7-a8)-a9*expg(-a10*((xx(2+mpa,j)/ct3)**8))+a11
     1*expg(-a12*(xx(2+mpa,j)/ct3-a13))
      g1=(1.0d0/ct3)*(-a2*a3/dcosh(-a3*xx(2+mpa,j)/ct3+a4)/dcosh(-a3
     1*xx(2+mpa,j)/ct3+a4)+a5*a7*(a6-xx(2+mpa,j)/ct3)**(-1.0d0+a7)+
     1a9*a10*8.0d0*((xx(2+mpa,j)/ct3)**7)*expg(-a10*
     1(xx(2+mpa,j)/ct3)**8))-a11*a12/ct3*expg(-a12*(xx(2+mpa,j)
     &/ct3-a13))
c
      if (g0.gt.4.5) then
      g0=4.5d0
      g1=0.0d0
c     write (2,*) 'U theta overflow - positive'
      else if (g0.lt.3.0) then
      g0=3.0d0
      g1=0.0d0
c     write (2,*) 'U theta underflow - positive'
      end if
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Tungsten oxide (LixWO3 with 0<x<0.67)
c     literature data from Whittingham
c
  205 c1=2.8767d0
      c2=-0.9046d0
      c3=0.76679d0
      c4=-0.15975d0
      c5=0.671d0
      g0=c1+c2*xx(2+mpa,j)/ct3+c3*xx(2+mpa,j)*xx(2+mpa,j)/ct3/ct3+
     1c4*expg(100.0d0*(xx(2+mpa,j)/ct3-0.671d0))
      g1=c2/ct3+2.0d0*c3*xx(2+mpa,j)/ct3/ct3+
     1c4*100.0d0*expg(100.0d0*(xx(2+mpa,j)/ct3-0.671d0))/ct3
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     CoO2 (Cobalt dioxide)
c
  206 continue
c     Marc Doyle's fit
c     r1=4.825510d0
c     r2=0.950237d0
c     r3=0.913511d0
c     r4=0.600492d0
c     g0=r1-r2*expg(-((xx(2+mpa,j)/ct3-r3)/r4)**2)
c     g1=2.0d0*r2*(xx(2+mpa,j)/ct3-r3)*expg(-((xx(2+mpa,j)/ct3
c    1-r3)/r4)**2)/r4/r4/ct3
c
c     Measured by Oscar Garcia 2001 using Quallion electrodes for
c     0.5 < y < 0.99.  Fit revised by Karen Thomas in May 2003 to
c     match Doyle's fit for y < 0.4 and Garcia's data at larger y.
c     Valid for 0 < y < 0.99. Note that capacity fade is found to
c     occur experimentally if y goes below 0.5; this is not included
c     in the model.
      stretch=1.062d0
      sto = xx(2+mpa,j)*stretch/ct3
      g0 = 2.16216d0+0.07645d0*dtanh(30.834d0-54.4806d0*sto)
     & + 2.1581d0*dtanh(52.294d0-50.294d0*sto)
     & - 0.14169d0*dtanh(11.0923d0-19.8543d0*sto)
     & + 0.2051d0*dtanh(1.4684d0-5.4888d0*sto)
     & + 0.2531d0*dtanh((-sto+0.56478d0)/0.1316d0)
     & - 0.02167d0*dtanh((sto-0.525d0)/0.006d0)
      g1 = 0.07645d0*(-54.4806d0/ct3)*
     &((1.0d0/dcosh(30.834d0-54.4806d0*sto))**2)
     &+2.1581d0*(-50.294d0/ct3)*((dcosh(52.294d0-50.294d0*sto))**(-2))
     &+0.14169d0*(19.854d0/ct3)*((dcosh(11.0923d0-19.8543d0*sto))**(-2))
     &-0.2051d0*(5.4888d0/ct3)*((dcosh(1.4684d0-5.4888d0*sto))**(-2))
     &-0.2531d0/0.1316d0/ct3*((dcosh((-sto+0.56478d0)/0.1316d0))**(-2))
     & - 0.02167d0/0.006d0/ct3*((dcosh((sto-0.525d0)/0.006d0))**(-2))
      dudt = (-0.19952+0.92837*sto-1.36455*sto**2+0.61154*sto**3)/
     &(1-5.66148*sto+11.47636*sto**2-9.82431*sto**3+3.04876*sto**4)
      dudt = dudt/1000.0d0 !V/K
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     V2O5 (Vanadium oxide 0<y<0.95)
c
  207 r1=3.3059d0
      r2=0.092769d0
      r3=14.362d0
      r4=6.6874d0
      r5=0.034252d0
      r6=100.0d0
      r7=0.96d0
      r8=0.00724d0
      r9=80.0d0
      r10=0.01d0
      a2=xx(2+mpa,j)/ct3
      g0=r1+r2*dtanh(-r3*a2+r4)-r5*expg(r6*(a2-r7))+r8*expg(r9*
     1(r10-a2))
      g1=-r2*r3/dcosh(-r3*a2+r4)/dcosh(-r3*a2-r4)/ct3-r5*r6*expg(r6*
     1(a2-r7))/ct3-r8*r9*expg(r9*(r10-a2))/ct3
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     NiO2 (Nickel dioxide 0.45<y<1.0)
c
  208 continue
c      From Marc Doyle  
c      r1=6.515d0
c      r2=2.3192d0
c      r3=5.3342d0
c      r4=0.41082d0
c      r5=200.0d0
c      r6=0.44d0
c      r7=0.24247d0
c      r8=60.0d0
c      r9=0.99d0
c      a3=xx(2+mpa,j)/ct3
c      g0=r1+r2*a3-r3*a3**0.5d0+r4*expg(r5*(r6-a3))-r7*expg(r8*(a3-r9))
c      g1=r2/ct3-0.5d0*r3*(xx(2+mpa,j)**(-0.5d0))/ct3**0.5d0
c     1-r4*r5*expg(r5*(r6-a3))/ct3-r7*r8*expg(r8*(a3-r9))/ct3
c      LiNi0.8Co0.2O2 (Gen 1) from PolyStor, measured by Chris Bogatu 2000
c      for 0.4 < y < 0.99
      sto = xx(2+mpa,j)/ct3
      g0 = 3.836d0-0.128d0*dtanh((sto-0.929d0)/0.123d0)
     & - 0.177d0*dtanh((sto-0.634d0)/0.123d0)
     & - 0.076d0*dtanh((sto-0.475d0)/0.0498d0)
      g1 = (-0.128d0/(ct3*0.123d0))*
     &((1.0d0/dcosh((sto-0.929d0)/0.123d0))**2)
     & - (0.177d0/(ct3*0.123d0))*
     &((1.0d0/dcosh((sto-0.634d0)/0.123d0))**2)
     & - (0.076d0/(ct3*0.0498d0))*
     &((1.0d0/dcosh((sto-0.475d0)/0.0498d0))**2)
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Spinel Mn2O4 (Bellcore 0.17<y<1.0)
c
  209 a1=4.19829d0
      a2=0.0565661d0
      a3=14.5546d0
      a4=8.60942d0
      a5=0.0275479d0
      a6=0.998432d0  ! would prefer this to be >=1
      a7=-0.492465d0
      a8=1.901110d0
      a9=0.157123d0
      a10=0.04738d0
      a11=0.810239d0
      a12=40.0d0
      a13=0.133875d0
cccc  if(xx(2+mpa,j).gt.a6*ct3) write (2,*) '#109 in ekin, j=',j
c     g0=a1+a2*dtanh(-a3*xx(2+mpa,j)/ct3+a4)-a5*((a6-xx(2+mpa,j)
c    &/ct3)**a7-a8)-a9*expg(-a10*((xx(2+mpa,j)/ct3)**8))+a11
c    1*expg(-a12*(xx(2+mpa,j)/ct3-a13))
      g0=a1+a2*dtanh(-a3*xx(2+mpa,j)/ct3+a4)
     1-a9*expg(-a10*((xx(2+mpa,j)/ct3)**8))+a11
     1*expg(-a12*(xx(2+mpa,j)/ct3-a13))+a5*a8
      if(xx(2+mpa,j).lt.a6*ct3) g0=g0-a5*((a6-xx(2+mpa,j)/ct3)**a7)

c     g1=(1.0d0/ct3)*(-a2*a3/dcosh(-a3*xx(2+mpa,j)/ct3+a4)/dcosh(-a3
c    1*xx(2+mpa,j)/ct3+a4)+a5*a7*(a6-xx(2+mpa,j)/ct3)**(-1.0d0+a7)+
c    1a9*a10*8.0d0*((xx(2+mpa,j)/ct3)**7)*expg(-a10*
c    1(xx(2+mpa,j)/ct3)**8))-a11*a12/ct3*expg(-a12*(xx(2+mpa,j)
c    &/ct3-a13))
      g1=(1.0d0/ct3)*(-a2*a3/dcosh(-a3*xx(2+mpa,j)/ct3+a4)/dcosh(-a3
     1*xx(2+mpa,j)/ct3+a4)
     1+a9*a10*8.0d0*((xx(2+mpa,j)/ct3)**7)*expg(-a10*
     1(xx(2+mpa,j)/ct3)**8))-a11*a12/ct3*expg(-a12*(xx(2+mpa,j)
     &/ct3-a13))
      if(xx(2+mpa,j).lt.a6*ct3)g1=g1+a5*a7*(a6-xx(2+mpa,j)/ct3)
     &**(-1.d0+a7)/ct3

cccc  if(xx(2+mpa,j).gt.a6*ct3) write (2,*) 'did it'
c
      if (g0.gt.6.0) then
      g0=6.0d0
      g1=0.0d0
c     write (2,*) 'U theta overflow - positive'
      else if (g0.lt.3.0) then
      g0=3.0d0
      g1=0.0d0
c     write (2,*) 'U theta underflow - positive'
      end if
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      go to 98
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Nonstoichiometric Vanadium oxide (V6O13) added by Karen Thomas March 14, 1999,
c based on data from West, Zachau-Christiansen, and Jacobsen,
c Electrochim. Acta vol 28, p. 1829, 1983.
c valid for 0.1 < x < 8.25 in LixV6O13.  Enter csx according to
c LiyVO2.167, where 0.05 <y <1, and cot3 is based on 8 Li inserted.
c Fit for electrical conductivity based on data from same
c paper, corrected for porosity.  Electrical conductivity
c of V6O13- carbon filler composite based on model of Meredith and
c Tobias in Advances in Electrochemistry and Electrochemical Engineering
c vol. 2, 1962.

 210  continue
c     convert solid concentration in mol/m3 to x in mol/mol V6O13
      sto = xx(2+mpa,j)*8.0d0/ct3
      a1 = -10.0d0/6.0d0
      a2 = 8.0d0/6.0d0/ct3
      g0 = 1.9d0 + 0.13d0*dtanh(a1*sto + 1.7d0)
     & + 0.2d0*dtanh(sto*a1 + 6.2d0)
     & + 0.56d0*dtanh(2.0d0*a1*sto + 29.5d0)
      g1=-(1.3d0*a2)/(dcosh(a1*sto+1.7d0)*dcosh(a1*sto+1.7d0))
     & -(2.0d0*a2)/(dcosh(a1*sto+6.2d0)*dcosh(a1*sto+6.2d0))
     & -(11.2d0*a2)/(dcosh(2.0d0*a1*sto+29.5d0)
     &*dcosh(2.0d0*a1*sto+29.5d0))
c     conductivity varies with state of charge for this material
      vanox=120.0d0*expg(-1.5d0*sto)-9.0d0*dtanh(1.5d0*sto-6.5d0)+9.0d0
      sigcarb = 100.0d0/0.0038d0 ! S/m
      boo = sigcarb/vanox
      sig3=vanox*(2.0d0*(boo+2.0d0)+2.0d0*(boo - 1.0d0)*epf3)
     & *((2.0d0-epf3)*(boo + 2.0d0)+2.0d0*(boo-1.0d0)*epf3)/
     & (2.0d0*(boo+2.0d0)-(boo-1.0d0)*epf3)/
     & ((2.0d0-epf3)*(boo+2.0d0)-(boo-1.0d0)*epf3)
       sig3 = sig3*(1.0d0 - ep3 - epp3)**1.5d0
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
       go to 98
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c LiyAl0.2Mn1.8O4F0.2 spinel from Bellcore.  Open circuit potential
c and dudt as function of y measured by Karen Thomas.
c Fit for cell 073 at 25 C, May 2000, for 0.21 < y < 1.0
c Stretch feature added May 28, 2013, to keep fractional charge < 1.d0.
 211  continue
       a1=4.918d0
       a2=0.032d0
       a3=12.7d0
       a4=7.4d0
       a5=0.0505d0
       a6=1.006d0
       a7=-0.49d0
       a8=0.085d0
       a9=0.685d0
       a10=0.24d0
       a11=0.2d0
       a12=26.0d0
       a13=0.24d0
       a14=26.0d0
      stretch=1.006d0 ! amount to shrink xx(3,j) by to keep the 
c     fractional state of charge < 1.0.
c     stretch=1.d0 ! no stretch
c
      xx(3,j)=xx(3,j)*stretch ! to give a value which can be > ct3.
      g0=a1+a2*dtanh(-a3*xx(2+mpa,j)/ct3+a4)
     1-a9*expg(-a10*((xx(2+mpa,j)/ct3)**a14))+a11
     1*expg(-a12*(xx(2+mpa,j)/ct3-a13))-(xx(2+mpa,j)/ct3)*a8
c     if(xx(2+mpa,j).lt.a6*ct3)
c    & g0=g0-a5*((a6-xx(2+mpa,j)/ct3)**a7)
      g0=g0-a5*((a6-xx(2+mpa,j)/ct3)**a7)
c
      g1=(1.0d0/ct3)*(-a2*a3)/dcosh(-a3*xx(2+mpa,j)/ct3+a4)/
     1 dcosh(-a3*xx(2+mpa,j)/ct3+a4)
     1 +a9/ct3*a10*a14*((xx(2+mpa,j)/ct3)**(a14-1.0d0))
     1 *expg(-a10*((xx(2+mpa,j)/ct3)**a14))
     1 -a11*a12/ct3*expg(-a12*(xx(2+mpa,j)/ct3-a13))
     1 -a8/ct3
c     if(xx(2+mpa,j).lt.a6*ct3)
c    & g1=g1+a5*a7*(a6-xx(2+mpa,j)/ct3)**(-1.0d0+a7)/ct3
      g1=g1+a5*a7*(a6-xx(2+mpa,j)/ct3)**(-1.0d0+a7)/ct3
c
      if (g0.gt.6.0d0) then
      g0=6.0d0
      g1=0.0d0
c     write (2,*) 'U theta overflow - positive'
c     else if (g0.lt.3.0d0) then
c     g0=3.0d0
c     g1=0.0d0
c     write (2,*) 'U theta underflow - positive'
      end if
c
      a1=-0.009d0
      alph1=12.0d0
      bet1=0.252d0
      a2 = -0.031d0
      alph2 = 9.0d0
      bet2 = 0.536d0
      a3 = 0.11d0
      a4 = -0.38d0
      a5 = 9.6d0
      a6 = 0.964d0
      a7 = -0.025d0
      alph3 = 8.5d0
      bet3 = 0.785d0
      a8 = -0.008d0
      alph4 = 8.0d0
      bet4 = 0.40d0
      sto = xx(2+mpa,j)/ct3
c
      dudt = a1*(alph1/(bet1**alph1))*(sto**(alph1-1.0d0))
     & *expg(-(sto/bet1)**alph1) +
     & a2*(alph2/(bet2**alph2))*(sto**(alph2-1.0d0))
     & *expg(-(sto/bet2)**alph2)
     & +a3 +a4*expg(a5*(sto-a6)) +
     & a7*(alph3/(bet3**alph3))*(sto**(alph3-1.0d0))
     & *expg(-(sto/bet3)**alph3) +
     & a8*(alph4/(bet4**alph4))*(sto**(alph4-1.0d0))
     & *expg(-(sto/bet4)**alph4)
      dudt = dudt/1000.0d0 !V/K
c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p
      xx(3,j)=xx(3,j)/stretch ! to restore xx to its correct value
      go to 98

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  212 continue
c NiOOHHy Nickel Positive electrode
c Fit by Paul Albertus, JES

c These numbers are for discharging
      c1=0.41d0
      c2=0.1d0
      c3=20.0d0
      c4=50.0d0
      c5=2.0d0
      c6=-0.016d0
      c7=2.0d0
      c8=5.0d0
      c9=0.015d0
      c10=1.07d0
      c11=-0.01d0
      c12=.01d0
      c13=2.7d0
      c14=0.9d0
      c15=0.01d0
      c16=0.65d0
      c17=50.0d0

c These numbers are for charging
      ca1=0.41d0
      ca2=2.0d0
      ca3=11.0d0
      ca4=20.0d0
      ca5=2.0d0
      ca6=-.008d0
      ca7=2.0d0
      ca8=0.1d0
      ca9=0.065d0
      ca10=1.00d0
      ca11=-0.3d0
      ca12=.001d0
      ca13=10.0d0
      ca14=0.7d0

      CNORM=xx(2+mpa,j)/ct3
      if (cur.ge.0.0d0) then !discharge
      g0=c1+c2*expg(c3*(-cnorm+c11))-c8*expg(c4*(cnorm-c10))
     &+c6*dlog(dabs((cnorm/(1.0d0-cnorm))))-c9-c12*expg(c13*cnorm-c14)
     &+c15*expg(-c17*((cnorm-c16)**2))
      g1=-c3*c2*expg(c3*(-cnorm+c11))-c8*c4*expg(c4*(cnorm-c10))
     &+c6/cnorm/(1.0d0-cnorm)-c12*c13*expg(c13*cnorm-c14)
     &+c15*expg(-c17*((cnorm-c16)**2))*(-c17)*2.0d0*(cnorm-c16)
      else !charge
      g0=ca1+ca2*expg(ca3*(-cnorm+ca11))-ca8*expg(ca4*(cnorm-ca10))
     &+ca6*dlog(dabs((cnorm/(1.0d0-cnorm))))+ca9-
     &ca12*expg(ca13*(cnorm-ca14))
      g1=-ca3*ca2*expg(ca3*(-cnorm+ca11))-ca8*ca4*expg(ca4*(cnorm-ca10))
     &+ca6/cnorm/(1.0d0-cnorm)-ca12*ca13*expg(ca13*(cnorm-ca14))
      endif
      g1=g1/ct3

      dudt=0.0d0

c kinetic expressions here for positive

      if (lag.eq.1) go to 99

      RTF=r*t/fc
      ATCN=0.58d0
      CTCN=0.35d0
      ATCP=0.36d0
      CTCP=0.58d0
      H2OMW=18.016D0
      PKOHMW=56.11D0
      Pcon=xx(1,j)/1.0d6
      Pcs=xx(2+mpa,j)/1.0d6
      CNGMAX=ct1/1.0d6
      CPSMAX=ct3/1.0d6

      EXCDN= 7.85D-04 /(0.012644D0**CTCN*0.046814D0**ATCN*
     1         (0.5D0*CNGMAX)**(CTCN+ATCN))
      EXCDP= 1.04D-04 /(0.012644D0**CTCP*0.046814D0**ATCP*
     1         (0.5D0*CPSMAX)**(CTCP+ATCP))     

      DN = 1.001D0 + 47.57D0*Pcon - 776.22D0*Pcon**2
      DN1D = 47.57D0 - 1552.44D0*Pcon
      DN2D = -1552.44D0

        AC = 1.004D0 - 36.23D0*Pcon**0.5 + 1374.3D0*Pcon
     1 - 17850.7*Pcon**1.5 + 55406.0D0*Pcon**2
     1 +7.16856D05*Pcon**2.5
        AC1D=-18.115D0*Pcon**(-0.5)+1374.3D0-2.6776D04*Pcon**0.5
     1 +1.10812D05*Pcon + 1.7921D06*Pcon**1.5
      AC2D = 9.0575D0/Pcon**1.5 - 13388.0D0/Pcon**0.5
     1 + 1.10812D05 + 2.6882D06*Pcon**0.5
      WAC = 1.0002D0 - 21.238D0*Pcon - 4.1312D03*Pcon**2.0D0
      WAC1D= -21.238D0 -2.0d0*4.1312D03*Pcon  

c assumed temperature dependence for exchange current
      ti0p=dexp((Ebarkc)*(t-298.15d0)/(t*298.15d0))
      rka3=rka3save*ti0p

      U0=rka3*EXCDP/H2OMW**ATCP
      U1=DN-Pcon*PKOHMW
      U2=DN1D-PKOHMW
      U3= (AC*Pcon)**CTCP*(WAC*U1)**ATCP
      U4= Pcs**CTCP*(CPSMAX-Pcs)**ATCP
      U6= (AC*Pcon)**CTCP*ATCP*(WAC*U1)**(ATCP-1.0D0)*
     1(WAC*U2+WAC1D*U1)+(WAC*U1)**ATCP*CTCP*
     1(AC*Pcon)**(CTCP-1.0D0)*(AC+Pcon*AC1D)
      U7= CTCP*(CPSMAX-Pcs)**ATCP*Pcs**(CTCP-1.0D0)-
     1ATCP*Pcs**CTCP*(CPSMAX-Pcs)**(ATCP-1.0D0)

      h0=U0*U3*U4/fc*1.0d4
      h1=U0*U3*U7/fc/1.0d2
      h2=U0*U6*U4/fc/1.0d2
      r1a=ATCP*frt
      r1c=CTCP*frt
      r2a=r1a*(xx(kp1,j)-xx(kp2,j)-g0)
      r2c=r1c*(xx(kp1,j)-xx(kp2,j)-g0)

      de=expg(-r2c)-expg(r2a)
      pe=expg(-r2c)+expg(r2a)


      b(2+mpa,1)=h2*de
      b(2+mpa,kp2)=h0*(r1c*expg(-r2c)+r1a*expg(r2a))
      b(2+mpa,kp1)=-b(2+mpa,kp2)
      b(2+mpa,2+mpa)=(h1*de+h0*(r1c*g1*expg(-r2c)+r1a*g1*expg(r2a)))
      b(2+mpa,kj)=1.0d0
      g(2+mpa)=-h0*de-xx(kj,j)

      go to 99 ! go to 99 not 97 to avoid other kinetics

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  213 continue

      write (2,*) 'Please enter data for positive electrode
     & #13 in subroutine ekin'
      stop

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     KINETIC EXPRESSIONS FOR THE POSITIVE ELECTRODE
c
c     h0 is the exchange current density (A/m2)
c     h1 is the derivative of io wrt solid concentration, xx(2+mpa,j)
c     h2 is the derivative of io wrt electrolyte concen., xx(1,j)
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     NONAQUEOUS LIQUIDS
c
   98 if (lag.eq.1) go to 99
      alpha=0.5d0
      alphc=0.5d0
      h0=rka3*dsqrt(xx(1,j))*dsqrt(ct3-xx(2+mpa,j))*dsqrt(xx(2+mpa,j))
      h1=rka3*dsqrt(xx(1,j))*dsqrt(ct3-xx(2+mpa,j))*dsqrt(xx(2+mpa,j))
     &*ct3/(ct3-xx(2+mpa,j))/xx(2+mpa,j)/2.0d0
      h2=rka3*dsqrt(ct3-xx(2+mpa,j))*dsqrt(xx(2+mpa,j))
     &/dsqrt(xx(1,j))/2.0d0
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     POLYMER
c
c     alpha=0.5d0
c     alphc=0.5d0
c     h0=rka3*dsqrt(xx(1,j))*dsqrt(cmax-xx(1,j))*dsqrt(ct3-xx(2+mpa,j))
c    1*dsqrt(xx(2+mpa,j))
c     h1=-rka3*dsqrt(xx(1,j))*dsqrt(cmax-xx(1,j))*dsqrt(ct3-xx(2+mpa,j))
c    1*dsqrt(xx(2+mpa,j))*(1.0d0/(ct3-xx(2+mpa,j))-1.0d0/xx(2+mpa,j))/2.0d0
c     h2=-rka3*dsqrt(xx(2+mpa,j))*dsqrt(ct3-xx(2+mpa,j))*dsqrt(cmax-xx(1,j))
c    1*dsqrt(xx(1,j))*(1.0d0/(cmax-xx(1,j))-1.0d0/xx(1,j))/2.0d0
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
      end if
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (lag.eq.1) go to 99
c
      r1=alpha*frt
      if(j.le.n1+1) then
           an1=1.0d0
           an2=0.0d0
         else
           an1=0.0d0
           an2=1.0d0
      endif
      ranode=ranodesave*dexp((Ebarr1)*(298.15d0-t)/(t*298.15d0))
      rcathde=rcathdesave*dexp((Ebarr3)*(298.15d0-t)/(t*298.15d0))

c Use the individual pore-wall flux for film resistance. PSA
      if (kk.le.2.or.mvdc1.eq.1.or.mvdc3.eq.1.or.izfl.eq.1) then
      r2=r1*(xx(kp1,j)-xx(kp2,j)-g0-fc*xx(kj,j)
     &*(an1*ranode+an2*rcathde))
      else 
      r2=r1*(xx(kp1,j)-xx(kp2,j)-g0-fc*pwfindiv(mpa,j)
     &*(an1*ranode+an2*rcathde))
      endif

      de=-2.d0*r2-r2**3/3.d0-r2**5/60.d0
      if(dabs(r2).gt.200.d0) then
      if(r2.gt.200.d0) de=7.0d86
      if(r2.lt.-200.d0) de=-7.0d86
      pe=7.0d86
      else
      if(dabs(r2).gt.1.0d-7) de=expg(-r2)-expg(r2)
      pe=expg(-r2)+expg(r2)
      endif

c     Try to put in a radial mass-transfer coefficient, so that the denominator
c     becomes 1 + coef/cold(1,j)*exp(-2r)
      coef=0.8d2

	if(t.gt.tam) coef=0.01d2

      denominator=1.d0+coef/xx(1,j)*expg(-r2)

c     if(j.eq.n1+n2+5) print *, 'denom ',denominator,c(1,j)*h0/coef
c    &,c(kj,j)

      b(2+mpa,1)=h2*de/denominator
      b(2+mpa,kp2)=h0*r1*pe/denominator
      b(2+mpa,kp1)=-b(2+mpa,kp2)
      b(2+mpa,2+mpa)=(h1*de+h0*r1*g1*pe)/denominator

c Modify Jacobians for individual pore-wall flux. PSA.

      if (kk.ge.2.and.mvdc1.eq.0.and.mvdc3.eq.0.and.izfl.eq.0) then
      b(2+mpa,2+mpa)=b(2+mpa,2+mpa)+
     &(fc*b(2+mpa,kp2)*(an1*ranode+an2*rcathde)*pwfindiv_d(mpa,j))
     &/denominator
      b(2+mpa,kj)=1.0d0
      else
      b(2+mpa,kj)=1.0d0+fc*b(2+mpa,kp2)*(an1*ranode+an2*rcathde)
     &/denominator
c     if(j.eq.n1+n2+5) print *, 'denom ',denominator,c(1,j)*h0/coef
      endif

      g(2+mpa)=-h0*de/denominator-xx(kj,j)
      b(2+mpa,kp2)=b(2+mpa,kp2)+h0*de/denominator**2
     &*(-coef/xx(1,j)*r1*expg(-r2))
	b(2+mpa,kp1)=-b(2+mpa,kp2)
      b(2+mpa,1)=b(2+mpa,1)+h0*de/denominator**2
     &*(coef/xx(1,j)**2*expg(-r2))

   99 return
      end

c***********************************************************************
      subroutine prop(nj,n2,n1)
c     subroutine to create library of properties of various
c     electrolytes.
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3,shutdown
      common/tprop/df(221),cd(221),tm(221),ddo2(221),ddh2(221),
     1ddf(221),dcd(221),dtm(221),dfu(221),d2fu(221),do2(221),dh2(221)
      common/temp/ thk,htc,dudt,Cp,dens,tam,g0,qq,qloss,residm,ncell,lht
      common/side/rksc1,c1init,c2init,rksa1,term_s1(221),vol,
     &rksa2,rksa3,rksc3,rksc2,UsO2,UsH2,cn2,term_s2(221),nside
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
c
      do 99 j=1,nj
      ee=ep2+epp2
      if(j .lt. n1+2 .and. n1 .gt. 0) ee=ep1+epp1
      if(j .gt. n2+n1) ee=ep3+epp3 

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14),nprop
    1 go to 101  ! AsF6 in methyl acetate
    2 go to 102  ! Perchlorate in PEO
    3 go to 103  ! Sodium Triflate in PEO
    4 go to 104  ! LiPF6 in PC (Sony cell simulation)
    5 go to 105  ! Perchlorate in PC (West's simulation)
    6 go to 106  ! Triflate in PEO
    7 go to 107  ! LiPF6 in EC/DMC and p(VdF-HFP) (Bellcore)
    8 go to 108  ! LiPF6 in EC/DMC and p(VdF-HFP) (Bellcore) cell #2
    9 go to 109  ! Ideal ionomer, t+ = 1.0
   10 go to 110  ! LiTFSI in PEMO (from Steve Sloop, 1999)
   11 go to 111  ! LiPF6 in EC:DMC
   12 go to 112  ! LiTFSI in PEO at 85 C (from Ludvig Edman, 2000)
   13 go to 113  ! 30% KOH in H20 (from Paxton)  
   14 go to 114  ! add your own
c
c  Key to labeling of transport properties:
c  As written salt must be a binary electrolyte with v+=v-=1
c  df(j) - Diffusion coefficient of the salt (m2/s)
c  ddf(j) - First derivative of df(j) wrt electrolyte concentration
c  cd(j) - Conductivity of the salt (S/m)
c  dcd(j) - First derivative of cd(j) wrt electrolyte concentration
c  tm(j) - Transference number. For system in which cation reacts use
c          t+ (e.g. Li systems), for systems in which anion reacts use
c          t- (e.g. NiMH system).
c  dtm(j) - First derivative of tm(j) wrt electrolyte concentration
c  dfu(j) - Activity factor for salt (dlnf/dc) 
c  d2fu(j) - Derivative of dfj(j) wrt electrolyte concentration (d2lnf/dc2)
c 
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     AsF6 in methyl acetate
c
c     diffusion coefficient of the salt (m2/s)
  101 df(j)=(ee**1.5d0)*1.54d-09
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     conductivity of the salt (S/m)
      cd(j)=2.5d0*(ee**1.5d0)
      dcd(j)=0.0d0
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c     transference number of lithium
      tm(j)=0.20d0
      dtm(j)=0.0d0
c     activity factor for the salt
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Perchlorate in PEO
c
c     diffusion coefficient of the salt (m2/s)
  102 df(j)=(ee**1.5d0)*1.78d-12
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     conductivity of the salt (S/m)
      cd(j)=1.6d-02*ee**1.5d0
      dcd(j)=0.0d0
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c     transference number of lithium
      tm(j)=0.10d0
      dtm(j)=0.0d0
c     activity factor for the salt
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Sodium Triflate in PEO
c
  103 r0=1.3041d-07
      r1=4.4978d-07
      r2=-3.1248d-07
      r3=-2.2383d-07
      r4=8.9264d-09
c
c     diffusion coefficient of the salt (m2/s)
      df(j)=0.0001d0*(ee**1.5d0)*(r0+r1*xx(1,j)/1000.d0+
     1r2*((xx(1,j)/1000.d0)**0.5d0) + r3*((xx(1,j)/1000.d0)**1.5d0)
     1+ r4*((xx(1,j)/1000.d0)**3))
      ddf(j)=0.0001d0*(ee**1.5d0)*(r1/1000.d0 +
     10.5d0*r2*(xx(1,j)**(-0.5d0))/(1000.0d0**0.5d0)
     1+ 1.5d0*r3*(xx(1,j)**0.5d0)/1000.d0**1.5d0 +
     13.0d0*r4*(xx(1,j)**2)/1000.d0**3)
      if (xx(1,j).ge.3.0d03) then
      df(j)=(ee**1.5d0)*1.6477d-12
      ddf(j)=0.0d0
      end if
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     conductivity of the salt (S/m)
      r7=4.32d-05
      r8=0.00017d0
      r9=0.000153d0
      r10=3.73d-05
c
      cd(j)=100.d0*(ee**1.5d0)*(r7+r8*xx(1,j)/1000.d0+
     &r9*xx(1,j)*xx(1,j)
     1/1000000.d0+r10*xx(1,j)*xx(1,j)*xx(1,j)/1000000000.d0)
      dcd(j)=100.d0*(ee**1.5d0)*(r8/1000.d0+2.0d0*r9*xx(1,j)/1.d06+
     &3.0d0*r10*xx(1,j)*xx(1,j)/1000000000.d0)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c
c     transference number of lithium
c
      if(xx(1,j).lt.0.3d03) then
      r5=0.32141d0
      r6=2.5768d0
      r11=71.369d0
      r12=643.63d0
      r13=1983.7d0
      r14=2008.d0
      r15=287.46d0
      tm(j)=r5-r6*xx(1,j)/1000.d0+r11*xx(1,j)*xx(1,j)/1000000.d0
     1-r12*((xx(1,j)/1000.d0)**3)+r13*((xx(1,j)/1000.d0)**4)
     1-r14*((xx(1,j)/1000.d0)**5)+r15*((xx(1,j)/1000.d0)**6)
      dtm(j)=-r6/1000.d0+2.0d0*r11*xx(1,j)/1000000.d0-
     13.0d0*r12*(xx(1,j)**2)/(1000.d0**3) +
     14.0d0*r13*(xx(1,j)**3)/(1000.d0**4) -
     15.0d0*r14*(xx(1,j)**4)/(1000.d0**5) +
     16.0d0*r15*(xx(1,j)**5)/(1000.d0**6)
      else
      tm(j)=0.0d0
      dtm(j)=0.0d0
      end if
c
      if(xx(1,j).ge.0.70d03) then
      r5=4.5679d0
      r6=4.506d0
      r11=0.60173d0
      r12=1.0698d0
      tm(j)=-r5+r6*expg(-((xx(1,j)/1000.d0-r11)/r12)**2)
      dtm(j)=-r6*(xx(1,j)/1000.d0-r11)*2.d0
     1*expg(-((xx(1,j)/1000.d0-r11)/r12)**2)/r12/r12/1000.d0
      end if
c
      if(xx(1,j).ge.2.58d03) then
      tm(j)=-4.4204d0
      dtm(j)=0.0d0
      end if
c
c     activity factor for the salt:  (dlnf/dc) and (d2lnf/dc2)
c
      if(xx(1,j).gt.0.45d03) then
      r17=0.98249d0
      r18=1.3527d0
      r19=0.71498d0
      r20=0.16715d0
      r21=0.014511d0
      thermf=r17-r18*xx(1,j)/1000.d0+r19*xx(1,j)*xx(1,j)/1000000.d0-
     1r20*xx(1,j)*xx(1,j)*xx(1,j)/1000000000.d0+r21*xx(1,j)*xx(1,j)
     1*xx(1,j)*xx(1,j)/1000000000000.d0
      dthermf=-r18/1000.d0+2.d0*r19*xx(1,j)/1000000.d0-
     13.d0*r20*xx(1,j)*xx(1,j)/1000000000.d0+4.d0*r21*xx(1,j)*xx(1,j)
     1*xx(1,j)/1000000000000.d0
      end if
c
      if(xx(1,j).le.0.45d03) then
      r23=0.99161d0
      r24=0.17804d0
      r25=55.653d0
      r26=303.57d0
      r27=590.97d0
      r28=400.21d0
      thermf=r23-r24*xx(1,j)/1000.d0-r25*xx(1,j)*xx(1,j)/1000000.d0+
     1r26*xx(1,j)*xx(1,j)*xx(1,j)/1000000000.d0-r27*xx(1,j)*xx(1,j)
     1*xx(1,j)*xx(1,j)/1000000000000.d0+r28*xx(1,j)*xx(1,j)*xx(1,j)
     1*xx(1,j)*xx(1,j)/1000000000000000.d0
      dthermf=-r24/1000.d0-2.d0*r25*xx(1,j)/1000000.d0+
     13.d0*r26*xx(1,j)*xx(1,j)/1000000000.d0-4.d0*r27*xx(1,j)
     1*xx(1,j)*xx(1,j)/1000000000000.d0+5.d0*r28*xx(1,j)*xx(1,j)
     1*xx(1,j)*xx(1,j)/1000000000000000.d0
      end if
c
      dfu(j)=(-1.d0+2.d0*thermf)/xx(1,j)
      d2fu(j)=1.d0/xx(1,j)/xx(1,j)-2.d0*thermf/xx(1,j)/xx(1,j)+
     12.d0*dthermf/xx(1,j)
c
      if(xx(1,j).ge.3.00d03) then
      dfu(j)=-0.9520d0/xx(1,j)
      d2fu(j)=0.9520d0/xx(1,j)/xx(1,j)
      end if
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     LiPF6 in PC (Sony cell simulation used in JES vol. 141, 1994, p. 982 paper)
c
c     this is actually the diff coeff for perchlorate
c     diffusion coefficient of the salt (m2/s)
  104 df(j)=(ee**1.5d0)*2.58d-10
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     conductivity of the salt (S/m) from Barthel et al., Ber. Bunsenges. Phys. Chem. vol 83, 1979, p. 911
c     pmax=0.5409d0 !Marc's fit to Sony cell data
      pmax=0.035d0 !from Barthel et al.
      pu=0.857d0
      aa=1.093d0
      bb=0.04d0
      rho=1.2041d03
      fun=pmax*((1.0d0/rho/pu)**aa)*expg(bb*((xx(1,j)/rho-pu)**2)
     1-(aa/pu)*(xx(1,j)/rho-pu))
      fun2=2.0d0*(bb/rho)*(xx(1,j)/rho-pu)-aa/pu/rho
      cd(j)=0.0001d0+(ee**1.5d0)*((xx(1,j))**aa)*fun
      dcd(j)=(ee**1.5d0)*fun*(aa*(xx(1,j)**(aa-1.0d0))+(xx(1,j)**aa)
     1*fun2)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c     transference number of lithium
      tm(j)=0.20d0
      dtm(j)=0.0d0
c
c     activity factor for the salt (dlnf/dc and d2lnf/dc2)
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Perchlorate in PC (West's simulation)
c
c     diffusion coefficient of the salt (m2/s)
  105 df(j)=(ee**1.5d0)*2.58d-10
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     conductivity of the salt (S/m)
      pmax=0.542d0
      pu=0.6616d0
      aa=0.855d0
      bb=-0.08d0
      rho=1.2041d03
      fun=pmax*((1.0d0/rho/pu)**aa)*expg(bb*((xx(1,j)/rho-pu)**2)
     1-(aa/pu)*(xx(1,j)/rho-pu))
      fun2=2.0d0*(bb/rho)*(xx(1,j)/rho-pu)-aa/pu/rho
      cd(j)=0.0001d0+(ee**1.5d0)*((xx(1,j))**aa)*fun
      dcd(j)=(ee**1.5d0)*fun*(aa*(xx(1,j)**(aa-1.0d0))+(xx(1,j)**aa)
     1*fun2)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c     transference number of lithium
      tm(j)=0.20d0
      dtm(j)=0.0d0
c     activity factor for the salt
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Triflate in PEO
c
  106 r0=-5.0891863844d-05
      r1=8.38645199394d-07
      r2=-5.19747901855d-10
      r3=8.0832709407d-14
c     diffusion coefficient of the salt (m2/s)
      df(j)=(ee**1.5d0)*7.5d-12
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     conductivity of the salt (S/m)
      cd(j)=(ee**1.5d0)*100.0d0*(r0 + r1*xx(1,j)
     &+r2*xx(1,j)*xx(1,j)+r3*xx(1,j)*xx(1,j)*xx(1,j))
      dcd(j)=(ee**1.5d0)*100.0d0*(r1 + 2.0d0*r2*xx(1,j) +
     &3.0d0*r3*xx(1,j)**2)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c     transference number of lithium
c     rough conc. dependence of t+ - highly suspect
      tm(j)=0.0107907d0 + 1.48837d-04*xx(1,j)
      dtm(j)=1.48837d-04
c     activity factor for the salt
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     LiPF6 in EC/DMC and p(VdF-HFP) (Bellcore)
c     This is the 1:2 v/v mixture of EC/DMC (eq. 2 of paper)
c     D and t+ given below were fit from discharge curves
c
c     diffusion coefficient of the salt (m2/s)
c 107 df(j)=(ee**1.5d0)*9.00d-11
  107 df(j)=(ee**3.3d0)*7.50d-11
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c
c     conductivity of the salt (S/m)
c     This is the conductivity of the liquid + salt only (no polymer)
c
c     kappa (c) for EC/DMC 2:1 with LiPF6 at 25 deg C
c     cd(j)=ee**1.5d0*(0.0911d0+1.9101d0*xx(1,j)/1000.d0-1.052d0*xx(1,j)
c    1*xx(1,j)/1000.d0/1000.d0+0.1554d0*(xx(1,j)/1000.d0)**3)
c     derivative of kappa (c) for EC/DMC 2:1 at 25 deg C
c     dcd(j)=(ee**1.5d0)*(1.9101d0/1000.d0-2.0d0*1.052d0*xx(1,j)
c    &/1000.d0/1000.d0+3.0d0*0.1554d0/1000.d0*(xx(1,j)/1000.d0)**2)
c
c     kappa (c) for EC/DMC 1:2 w/ LiPF6 at 25 deg C
c     Note Bruggeman exponent should be adjusted to account for
c     polymer phase - this also affects "fac" parameter in Ohm's
c     law equation number 2
      r1=0.00010793d0
      r2=0.0067461d0
      r3=0.0052245d0
      r4=0.0013605d0
      r5=0.00011724d0
      cd(j)=(ee**3.3d0)*(r1+r2*xx(1,j)/1000.d0
c     cd(j)=(ee**1.5d0)*(r1+r2*xx(1,j)/1000.d0
     1-r3*xx(1,j)*xx(1,j)/1000000.d0
     1+r4*(xx(1,j)/1000.d0)**3-r5*(xx(1,j)/1000.d0)**4)*100.d0
      dcd(j)=(ee**3.3d0)*(r2/10.d0-r3*2.0d0*xx(1,j)/10000.d0
c     dcd(j)=(ee**1.5d0)*(r2/10.d0-r3*2.0d0*xx(1,j)/10000.d0
     1+3.0d0*r4*xx(1,j)*xx(1,j)/10000000.d0
     1-0.4d0*r5*(xx(1,j)/1000.d0)**3)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c
c     transference number of lithium
      tm(j)=0.363d0
      dtm(j)=0.0d0
c
c     activity factor for the salt (dlnf/dc and d2lnf/dc2)
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     LiPF6 in EC/DMC and p(VdF-HFP) (Bellcore) cell #2
c     This is the 2:1 v/v mixture of EC/DMC (eq. 1 of paper)
c     D and t+ given below were fit from discharge curves
c
c     diffusion coefficient of the salt (m2/s)
c 108 df(j)=(ee**1.5d0)*9.00d-11
  108 df(j)=(ee**3.3d0)*7.50d-11
      ddf(j)=0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c
c     conductivity of the salt (S/m)
c
c     kappa (c) for EC/DMC 2:1 w/ LiPF6 at 25 deg C
c     Note Bruggeman exponent should be adjusted to account for
c     polymer phase - this also affects "fac" parameter in Ohm's
c     law equation number 2
      r1=0.00041253d0
      r2=0.005007d0
      r3=0.0047212d0
      r4=0.0015094d0
c     r5=0.0016018d0
      r5=0.00016018d0
      cd(j)=(ee**3.3d0)*(r1+r2*xx(1,j)/1000.d0
c     cd(j)=(ee**1.5d0)*(r1+r2*xx(1,j)/1000.d0
     1-r3*xx(1,j)*xx(1,j)/1000000.d0
     1+r4*(xx(1,j)/1000.d0)**3-r5*(xx(1,j)/1000.d0)**4)*100.d0
      dcd(j)=(ee**3.3d0)*(r2/10.d0-r3*2.0d0*xx(1,j)/10000.d0
c     dcd(j)=(ee**1.5d0)*(r2/10.d0-r3*2.0d0*xx(1,j)/10000.d0
     1+3.0d0*r4*xx(1,j)*xx(1,j)/10000000.d0
     1-0.4d0*r5*(xx(1,j)/1000.d0)**3)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
c
c     transference number of lithium
      tm(j)=0.363d0
      dtm(j)=0.0d0
c
c     activity factor for the salt (dlnf/dc and d2lnf/dc2)
      dfu(j)=0.0d0
      d2fu(j)=0.0d0
      go to 99

c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Ideal Ion Exchange Membrane, t+ = 1.0
c
  109 continue
      df(j) = (ee**1.5d0)*1.0d-11
      ddf(j) = 0.0d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
c     Note: df chosen to be consistent with reported values
c     cd = 0.01d0 S/m (USABC's reported goals)
      cd(j) = (ee**1.5d0)*0.01d0
      dcd(j) = 0.0d0
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
      tm(j) = 1.0d0
      dtm(j) = 0.0d0
      dfu(j) = 0.0d0
      d2fu(j) = 0.0d0
      go to 99

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     USABC Ideal Polymer, t+ = 0.3
c
c      df(j) = (ee**1.5d0)*1.0d-11
c      ddf(j) = 0.0d0
c      cd(j) = (ee**1.5d0)*0.1d0
c      dcd(j) = 0.0d0
c      tm(j) = 0.3d0
c      dtm(j) = 0.0d0
c      dfu(j) = 0.0d0
c      d2fu(j) = 0.0d0
c     LiTFSI in PEMO from Steve Sloop at 40 C
  110 continue
      df(j) = (ee**1.5d0)*((-3.0d-17)*xx(1,j) + 6.0d-13)
      ddf(j) = (ee**1.5d0)*(-3.0d-17)
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
      cd(j) = (ee**1.5d0)*expg((-1.6d-06)*xx(1,j)*xx(1,j)
     & + 3.5d-03*xx(1,j) - 5.9d0)
      dcd(j) = cd(j)*((-3.2d-06)*xx(1,j) + 3.5d-03)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
      tm(j) = 0.6991d0-0.0004d0*xx(1,j)
      dtm(j) = -0.0004d0
      dfu(j) = 0.0021d0
      d2fu(j) = 0.0d0
      go to 99
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c LiPF6 in EC:DMC.
c D and t+ are fit to data from J. Power Sources, vol. 82, 1999, p. 859
c for LiPF6 in EC:EMC.
c cd is from measurements made at
c Bellcore, as reported in Marc Doyle's dissertation.
c
  111 continue
c      exbrug1=1.5d0
c      exbrug2=1.5d0
c      exbrug3=1.5d0
c     Need to set 300 in the following line to the initial temperature + 1.
      if(t.le.tam) Too=tam
      if(t.gt.Too) Too=t
      if(n1.gt.0 .and. j.le.n1+1) exbrug=exbrug1
      if(j.gt.n1 .and. j.le.n1+n2) exbrug=exbrug2
c     The following two lines are for use when a shutdown separator is 
c     being used.  Comment then out when no shut-down is wanted.
      if(shutdown.eq.1.0) then
      if(j.gt.n1 .and. j.le.n1+n2) ee=ee/(1+2.d3*(1+dtanh(((Too*0.00999
     &-3.7)-0.33)/0.010)))
      end if
      if(j.eq.n1+n2+1) exbrug=exbrug3

      df(j) = (ee**exbrug)*5.34d-10*expg(-0.65d0*xx(1,j)/1000.0d0)
      ddf(j) = -0.65d0*df(j)/1000.d0
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
      cd(j) = (ee**exbrug)*(0.0911d0+1.9101d0*xx(1,j)/1000.0d0 -
     & 1.052d0*((xx(1,j)/1000.0d0)**2) +
     & 0.1554d0*((xx(1,j)/1000.0d0)**3))
      dcd(j) = (ee**exbrug)*(1.9101d0/1000.0d0 -
     & 2.0d0*1.052d0*xx(1,j)/1000.0d0/1000.0d0
     & + 0.1554d0*3.0d0*((xx(1,j)/1000.0d0)**2)/1000.0d0)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
      tm(j) =0.4d0
      dtm(j) = 0.0d0
      dfu(j) = 0.0d0
      d2fu(j) = 0.0d0
      go to 99
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c LiTFSI in PEO at 85 C, from Ludvig Edman (LBL) July 2000
  112 continue
      df(j) = (ee**1.5d0)*(-6.73d-19*xx(1,j)*xx(1,j)+
     & 4.92d-16*xx(1,j) +5.5d-12)
      ddf(j) = (ee**1.5d0)*(-13.4d-19*xx(1,j)+4.92d-16)
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd
      cd(j)=(ee**1.5d0)*(6.53d-14*(xx(1,j)**3)-5.73d-10*xx(1,j)*
     & xx(1,j) +1.2d-06*xx(1,j)+4.25d-05)
      dcd(j) = (ee**1.5d0)*(19.59d-14*xx(1,j)*xx(1,j)-11.4d-10*xx(1,j)
     & + 1.2d-06)
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap
      tm(j) = -5.05d-08*xx(1,j)*xx(1,j) + 3.77d-04*xx(1,j)-0.0834d0
      dtm(j) = -10.1d-08*xx(1,j) + 3.77d-04
      dfu(j) = 0.0d0
      d2fu(j) = 0.0d0
      go to 99
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c 30% KOH in H2O.  Data from Blaine Paxton's masters thesis
c note the unit conversions

  113 continue

      do2(j)=(ee**1.5d0)*1.0d-5 !gas-phase oxygen diffusion coefficient
      ddo2(j)=0.0d0  ! first derivative wrt O2 concentration
      dh2(j)=(ee**1.5d0)*1.0d-5
      ddh2(j)=0.0d0  ! first derivative wrt H2 concentration

      c1=xx(1,j)/1.0d6
      DF(J) = (ee**1.5d0)*(2.8509D-5 -2.9659D-4*c1**0.5D0
     1+1.3768D-2*c1-0.14199D0*c1**1.5D0+0.42661D0*c1**2.0D0)
      DF(J)=DF(J)/1.0d4
      DDF(J) = (-7.4148D-05*C1**(-0.5D0)+1.3768D-2 -0.212985D0
     1*C1**0.5D0+0.85322D0*C1)*(ee**1.5d0)
      DDF(J) = DDF(J)/1.0d10
      tdd=dexp((EbarD)*(t-298.15d0)/(t*298.15d0))
      df(j)=df(j)*tdd
      ddf(j)=ddf(j)*tdd

      CD(J) = 2.325D-02 + 210.95D0*C1 -2.2077D04*C1**2.0D0
     1  +6.2907D05*C1**3.0D0
      CD(J) = CD(J)*ee**1.5D0*1.0d2
      DCD(J) = 210.95D0 -4.4154D04*C1+1.8872D06*C1**2.0D0
      DCD(J) = DCD(J)*ee**1.5D0/1.0d4
      tdkap=dexp((Ebarkap)*(t-298.15d0)/(t*298.15d0))
      cd(j)=cd(j)*tdkap
      dcd(j)=dcd(j)*tdkap


c For an anion involved in the electrode reaction put the 
c transference number of the anion here.  
c For KOH solution, t+=0.23 and t-=0.77
      tm(j)=0.77d0  
      dtm(j)=0.0d0

c U1=f, U2=df/dc (in Paxton units), U3=d2f/dc2 (in Paxton units)

      U1 = 1.004D0 - 36.23D0*C1**0.5d0 + 1374.3D0*C1
     &- 17850.7d0*C1**1.5d0 + 55406.0D0*C1**2.0d0
     &+7.16856D5*C1**2.5d0   
      U2=-18.115*C1**(-0.5d0)+1374.4d0-2.6776d4*C1**0.5d0+110812.0d0*C1
     &+1.7921d6*C1**1.5d0
      U3=9.0575*C1**(-1.5d0)-13388.0d0*C1**(-0.5d0)+110812.0d0+2.6881d6
     &*C1**0.5d0
      dfu(j)=1/U1*U2/1.0d6
      d2fu(j)=(1/U1*U3-U2*U2/(U1**2.0d0))/1.0d12

      go to 99

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c add your own
  114 continue
      if (nprop .eq. 14) then
      write (2,*) 'Please enter data for electrolyte
     & #14 in subroutine prop'
      stop
      endif
c      df(j) = (ee**1.5d0)* (diffusion coef. in m^2/s)
c      ddf(j) = d(df)/d(xx(1,j))
c      cd(j) = (ee**1.5d0)* (conductivity in S/m)
c      dcd(j) =  d(cd)/d(xx(1,j))
c      tm(j) = (cation transference number)
c      dtm(j) = (d(tm)/d(xx(1,j))
c      dfu(j) = 0.0d0
c      d2fu(j) = 0.0d0
      go to 99
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   99 continue
c
      return
      end

c***********************************************************************

      subroutine vardc(j,jcount,kk)
      implicit real*8(a-h,o-z)
      parameter(maxt=900)
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      common /calc/ ai(maxt,5),ai2(maxt,5),ts(maxt),h,h1,h2,h3,hcn,
     1hcp,rr,rrmax,cuL,sumpa(5,221),Rad1pa(5),Rad3pa(5),area1pa(5)
     &,area3pa(5),mcL
      common/pindiv/pwfindiv(5,221),pwfindiv_d(5,221),vf1(5),vf3(5),izfl
      common/const/ fc,r,t,frt,cur,ep3,ep2,pi,ep1,epf3,epf1,
     &epp1,epp2,epp3,shape3,shape1,capp1,capp3,nneg,nprop,npos
      common/gas/ epg1,epg2,epg3
      common/var/ xp(16),xx(17,221),xt(16,221,maxt)
     &,exbrug,exbrug1,exbrug2,exbrug3
      common/activ/ EbarD,Ebarkap,Ebarka,Ebarkc,Ebarr1,Ebarr3,
     &Ebars1,Ebars3,Ebarks1a,Ebarks1c,Ebarks2a,Ebarks2c,
     &Ebarks3a,Ebarks3c
      common /n/ tmmax,imp,ji,nx,nt,n1,n2,nj,n3,nconv,npa,iSk3,kj3
     &,ii2,ki2,kj,i2div,ip2,kp2,ip1,kp1,imb2,kS2,imb1,kS1,iSk1,kj1
     &,iSk2,kj2
      common/cprop/ sig3,area3,rka3,rka3save,ct3,dfs3,Rad3,cap1,cap3,
     1sig1,area1,rka1,rka1save,ct1,dfs1,Rad1,tw,dfs1save,dfs3save
      dimension ds1d(150),oflx(221),ocsufs(221)
      save

      errlimit2=1.0d-10 !max error for convergence

C SOLVE FOR THE CONC PROFILE FOR ALL X

C RESET CHANGE VARIABLES
      DO 376 JJ=1,NNJ
      CC(1,JJ)=0.0D0
  376 CONTINUE

C CONVERGENCE LOOP STARTS HERE
      DO 112 MM=1,lims  !go up to lims tries on convergence

C SOLVE FOR DS AS A FUNCTION OF CONC
      DO 114 JJ=1,NNJ
      IF (jcount.eq.1.and.mm.eq.1) then 
      cssold(j,jj)=css(j,jj)
      dsold(j,jj)=ds(j,jj)
      endif
c For a true variable solid-phase diffusion coefficient put in the 
c functional dependence here
      if (j.le.n1+1) then
      ds(j,jj)=dfs1
      ds1d(jj)=0.0d0
      elseif (j.ge.n1+n2) then
      ds(j,jj)=dfs3
      ds1d(jj)=0.0d0
      endif

  114 CONTINUE

      ERR2=0.0D0
      DO 888 JJ=1,NNJ
C ASSIGN VALUES TO BAND2
      if (j.le.n1+1) hh=hha
      if (j.ge.n1+n2) hh=hhc
      DIST=HH*(JJ-1)
      DP=DIST+HH/2.0D0
      DM=DIST-HH/2.0D0
      IF (JJ.EQ.1) THEN  !jj=1 is at r=0
      BB(1,1)=HH**3.0D0/rr/6.0D0+0.25D0*HH*(DS(J,JJ)+
     &DS(J,JJ+1)-(CSS(J,JJ+1)-CSS(J,JJ))*DS1D(JJ))
      DD(1,1)=-0.25D0*HH*(DS(J,JJ)+DS(J,JJ+1)
     &+(CSS(J,JJ+1)-CSS(J,JJ))*DS1D(JJ+1))
      GG(1)=-HH**3.0D0*(CSS(J,JJ)-CSSOLD(J,JJ))/6.0D0
     &/rr+0.25D0*HH*(DS(J,JJ+1)+DS(J,JJ))
     &*(CSS(J,JJ+1)-CSS(J,JJ))+0.25D0*HH*(DSOLD(J,JJ+1)
     &+DSOLD(J,JJ))*(CSSOLD(J,JJ+1)-CSSOLD(J,JJ))
      ELSEIF (JJ.EQ.NNJ) THEN !jj=nnj is at particle surface
      AA(1,1)=0.25D0*DM**2.0D0*(-DS(J,JJ)-DS(J,JJ-1)+
     &(CSS(J,JJ)-CSS(J,JJ-1))*DS1D(JJ-1))/HH
      BB(1,1)=0.5D0*DM**2.0D0*HH/rr+0.25D0*DM**2.0D0*
     &(DS(J,JJ)+DS(J,JJ-1)+(CSS(J,JJ)-CSS(J,JJ-1))*
     &DS1D(JJ))/HH
      GG(1)=-0.5D0*DM**2.0D0*HH*(CSS(J,JJ)-CSSOLD(J,JJ))/
     &rr-0.25D0*DM**2.0D0*(DS(J,JJ)+DS(J,JJ-1))*
     &(CSS(J,JJ)-CSS(J,JJ-1))/HH-0.25D0*DM**2.0D0*
     &(DSOLD(J,JJ)+DSOLD(J,JJ-1))*(CSSOLD(J,JJ)
     &-CSSOLD(J,JJ-1))/HH-((xx(kj,j)+xt(kj,j,kk-1))/2.0d0
     &+(xx(kj2,j)+xt(kj2,j,kk-1))/2.0d0)*DIST**2.0D0
      ELSE
      AA(1,1)=0.25D0*DM**2.0D0*(-DS(J,JJ)-DS(J,JJ-1)
     &+(CSS(J,JJ)-CSS(J,JJ-1))*DS1D(JJ-1))/HH
      BB(1,1)=HH*DIST**2.0D0/rr+0.25D0*DM**2.0D0*
     &(DS(J,JJ)+DS(J,JJ-1)+(CSS(J,JJ)-CSS(J,JJ-1))*
     &DS1D(JJ))/HH+0.25D0*DP**2.0D0*(DS(J,JJ)+DS(J,JJ+1)
     &-(CSS(J,JJ+1)-CSS(J,JJ))*DS1D(JJ))/HH
      DD(1,1)=-0.25D0*DP**2.0D0*(DS(J,JJ)+DS(J,JJ+1)
     &+(CSS(J,JJ+1)-CSS(J,JJ))*DS1D(JJ+1))/HH
      GG(1)=-HH*DIST**2.0D0*(CSS(J,JJ)-CSSOLD(J,JJ))/rr
     &+0.25D0*(-DM**2.0D0*(DS(J,JJ)+DS(J,JJ-1))*(CSS(J,JJ)
     &-CSS(J,JJ-1))+DP**2.0D0*(DS(J,JJ+1)+DS(J,JJ))*
     &(CSS(J,JJ+1)-CSS(J,JJ)))/HH
     &+0.25D0*(-DM**2.0D0*(DSOLD(J,JJ)+DSOLD(J,JJ-1))*
     &(CSSOLD(J,JJ)-CSSOLD(J,JJ-1))+DP**2.0D0*
     &(DSOLD(J,JJ+1)+DSOLD(J,JJ))*(CSSOLD(J,JJ+1)
     &-CSSOLD(J,JJ)))/HH

      ENDIF

      ERR2=ERR2+DABS(GG(1))
      CALL BAND2(JJ)

888   CONTINUE !run until hit nnj

      DO 322 JJ=1,NNJ
      CSS(J,JJ)=CSS(J,JJ)+CC(1,JJ)
  322 CONTINUE

c set a shoe-horn for the solid concentration

      do l=1,nj
      do ll=1,nnj
      if (l.le.n1+1) then
      if (css(l,ll).lt.cssold(l,ll)/1.d2) css(l,ll)=cssold(l,ll)/1.d2
      if (css(l,ll).gt.ct1) css(l,ll)=ct1*0.999999d0
      elseif (l.ge.n1+n2) then
      if (css(l,ll).lt.cssold(l,ll)/1.d2) css(l,ll)=cssold(l,ll)/1.d2
      if (css(l,ll).gt.ct3) css(l,ll)=ct3*0.999999d0
      endif
      enddo
      enddo

C CHECK TO SEE IF CONCENTRATIONS CONVERGED
      ERR2=ERR2/NNJ

      IF (ERR2.LT.errlimit2) GOTO 345

  112 CONTINUE
      PRINT *,'SOLID PHASE DIFF NOT CONVERGE',kk,j

      STOP


  345 CONTINUE

      return 
      end

c***********************************************************************


      subroutine band2(j)
      implicit real*8(a-h,o-z)
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nn,nnj,np,mvdc1,mvdc3,lims
      dimension e(4,5,200)
  101 format (15h determ=0 at j=,i4)
      if (j-2)  1,6,8
    1 np1= nn + 1
      do 2 i=1,nn
      dd(i,2*nn+1)= gg(i)
      do 2 l=1,nn
      lpn= l + nn
    2 dd(i,lpn)= xxx(i,l)
      call matinv2(nn,2*nn+1,determ)
      if (determ)  4,3,4
    3 print 101, j
    4 do 5 k=1,nn
      e(k,np1,1)= dd(k,2*nn+1)
      do 5 l=1,nn
      e(k,l,1)= - dd(k,l)
      lpn= l + nn
    5 xxx(k,l)= - dd(k,lpn)
      return
    6 do 7 i=1,nn
      do 7 k=1,nn
      do 7 l=1,nn
    7 dd(i,k)= dd(i,k) + aa(i,l)*xxx(l,k)
    8 if (j-nnj)  11,9,9
    9 do 10 i=1,nn
      do 10 l=1,nn
      gg(i)= gg(i) - yy(i,l)*e(l,np1,j-2)
      do 10 m=1,nn
   10 aa(i,l)= aa(i,l) + yy(i,m)*e(m,l,j-2)
   11 do 12 i=1,nn
      dd(i,np1)= - gg(i)
      do 12 l=1,nn
      dd(i,np1)= dd(i,np1) + aa(i,l)*e(l,np1,j-1)
      do 12 k=1,nn
   12 bb(i,k)= bb(i,k) + aa(i,l)*e(l,k,j-1)
      call matinv2(nn,np1,determ)
      if (determ)  14,13,14
   13 print 101, j
   14 do 15 k=1,nn
      do 15 m=1,np1
   15 e(k,m,j)= - dd(k,m)
      if (j-nnj)  20,16,16
   16 do 17 k=1,nn
   17 cc(k,j)= e(k,np1,j)
      do 18 jj=2,nnj
      m= nnj - jj + 1
      do 18 k=1,nn
      cc(k,m)= e(k,np1,m)
      do 18 l=1,nn
   18 cc(k,m)= cc(k,m) + e(k,l,m)*cc(l,m+1)
      do 19 l=1,nn
      do 19 k=1,nn
   19 cc(k,1)= cc(k,1) + xxx(k,l)*cc(l,3)
   20 return
      end
      
c***********************************************************************

      subroutine matinv2(nn,m,determ)
      implicit real*8(a-h,o-z)
      COMMON /vdc/ aa(1,1),bb(1,1),cc(1,150),dd(1,3),gg(1),xxx(1,1),
     1yy(1,1),hha,hhc,cssold(220,150),css(220,150),utz(5,221),
     &dsold(220,150),ds(220,150),time,nx,nnj,np,mvdc1,mvdc3,lims
      dimension id(4)
      determ=1.0
      do 1 i=1,nn
   1  id(i)=0
      do 18 nm=1,nn
      bmax=1.1
      do 6 i=1,nn
      if(id(i).ne.0) go to 6
      bnext=0.0
      btry=0.0
      do 5 j=1,nn
      if(id(j).ne.0) go to 5
      if(dabs(bb(i,j)).le.bnext) go to 5
      bnext=dabs(bb(i,j))
      if(bnext.le.btry) go to 5
      bnext=btry
      btry=dabs(bb(i,j))
      jc=j
   5  continue
      if(bnext.ge.bmax*btry) go to 6
      bmax=bnext/btry
      irow=i
      jcol=jc
   6  continue
      if(id(jc).eq.0) go to 8
      determ=0.0
      return
   8  id(jcol)=1
      if(jcol.eq.irow) go to 12
      do 10 j=1,nn
      save=bb(irow,j)
      bb(irow,j)=bb(jcol,j)
  10  bb(jcol,j)=save
      do 11 k=1,m
      save=dd(irow,k)
      dd(irow,k)=dd(jcol,k)
  11  dd(jcol,k)=save
  12  f=1.0/bb(jcol,jcol)
      do 13 j=1,nn
  13  bb(jcol,j)=bb(jcol,j)*f
      do 14 k=1,m
  14  dd(jcol,k)=dd(jcol,k)*f
      do 18 i=1,nn
      if(i.eq.jcol) go to 18
      f=bb(i,jcol)
      do 16 j=1,nn
  16  bb(i,j)=bb(i,j)-f*bb(jcol,j)
      do 17 k=1,m
  17  dd(i,k)=dd(i,k)-f*dd(jcol,k)
  18  continue
      return
      end

c***********************************************************************
c     That's All Folks!
c***********************************************************************

c     1.  Radial mass transfer has been added to the program to allow very high
c         currents to be run.  This addition reflects the reality that lithium ions
c         cannot move at infinite rates to the reaction site.
c     2.  An additional change has been made that corrects the lithium-ion balance if 
c         concentration  becomes negative.  This can occur only at very high currents.  However, with
c         the radial mass transfer added in 1, this change may not be required any longer.
c     3.  A value for the internal cell resistance, RG, has been added to the program.
c         This value is broken down into two major portions.  The first part is for the
c         resistance that is internal to the can but outside the cell stack.  The second
c         part is that portion which is internal to the cell stack.  This quantity is
c         further divided into two parts (one for the negative and one for the positive).
c     3.  EX (exbrug exponent) is now set separately for the negative (exbrug1), the separator 
c         (exbrug2), and the positive (exbrug3). EX has been shown to range from 1.5 to 3. This 
c         insight comes from information in the literature and from of other battery systems.
c     4   Using information from the literature, we have added a function for a special 
c         separator (called a shut-down separator) whose porosity goes to zero above a given 
c         temperature range.
c     5.  We need to add a vapor pressure vs. temperature function for a given solution.
c     6.  We need a function for heat generation caused by the temperature rising above the
c         stability of the anode film.
c     7.  The program can now be run with very low load and RG resistances.  The 
c         resulting current is ~7,200 A/mm2.  We can get these high current
c         results only by running at least two legs, with the first leg giving a better
c         estimate for the time zero calculation.  This approach can be followed 
c         any time the zero-time solution does not converge.
c     8.  We now can use constant current for a short period of time (e.g. 0.00001 s)
c         before switching to the next leg.  This approach can be used with constant load
c         by running the constant current first.  This involves a little bit of trial and
c         error approach.

