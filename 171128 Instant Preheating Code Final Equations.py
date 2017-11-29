# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 15:00:42 2017

@author: Charlotte
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import integrate, special
import cPickle as pickle
from matplotlib import rc
from sympy import pprint
import os

plt.rcParams['xtick.labelsize'] = 30 
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Computer Modern"
rc('text', usetex=True)
#Numerical control
tolerance = 10**(-13)
app = 1
mpl = 100.0 
#------------------------------------------------------------------------------
#----Variables----
#------------------------------------------------------------------------------
phi = sp.symbols('phi')
m = sp.symbols('m')
#N_star = 60.0 #55.0 # 63.49
As = 2.208e-9
#n = 123.0
#alpha = 0.01
#radiation_initial_choice = 0.0001
#Integration set-up IP
start1 = 0.0
end1 = 2.0 * 10**(10) / mpl
numsteps1 = 8000000
T1 = np.linspace(start1,end1,numsteps1)
#Integration set-up Radiation
start2 = 0.0
end2 = 1.0 * 10**(6) / mpl
numsteps2 = 8000000
T2 = np.linspace(start2,end2,numsteps2)

#------------------------------------------------------------------------------
#----Equations----
#------------------------------------------------------------------------------
def Get_Equations(n,N_star):
    M = mpl * sp.sqrt((3.0*np.pi*sp.sqrt(2.0*alpha*As)/(N_star+(sp.sqrt(3.0*alpha)/2.0)))*sp.exp((3.0*alpha)/(4.0*(N_star+(sp.sqrt(3.0*alpha)/2.0)))))
    Vfun = M**4.0*sp.exp(-2.0*n*sp.exp((2.0*phi)/(sp.sqrt(6.0*alpha)*mpl)))
    V = sp.lambdify(phi,Vfun,"numpy") #Evaluatable potential, ie V(0.1) returns a number
    dV = sp.lambdify(phi,sp.diff(Vfun,phi,1),"numpy") #Derivative
    ddV = sp.lambdify(phi,sp.diff(Vfun,phi,2),"numpy")
    ep_minus_1 = sp.lambdify(phi,(((mpl**2.0 / 2.0) * (((sp.diff(Vfun,phi,1))/Vfun)**2.0))-1.0),"numpy")
    phi_star = (sp.sqrt(6.0*alpha)*mpl/2.0)*(sp.ln((3.0*alpha)/(4.0*n*(N_star+(sp.sqrt(3.0*alpha)/2.0)))))
    integral = sp.lambdify(phi,(Vfun/(sp.diff(Vfun,phi,1))),"numpy") 
    return [M,Vfun,V,dV,ddV,ep_minus_1,phi_star,integral]

##############################################
#### NO RADIATION INTEGRATION SET UP
##############################################
#-----------------------------------------------------------------------------
# Initial conditions
#------------------------------------------------------------------------------
def Get_IC_no_rad(phi_star):
    phi0 = float(phi_star)
    phidot0 = 0.0
    H0 =  ((phidot0**2/(6.0*mpl**2)) + (V(phi0)/(3.0*mpl**2)))**0.5
    ICs_no_rad = [phi0,phidot0,H0]
    return ICs_no_rad
#------------------------------------------------------------------------------
# Definition of the background differential equations
#------------------------------------------------------------------------------
def Hdotsolver_no_rad(phidot):
    Hdot = -phidot**2/(2*mpl*mpl)
    return Hdot

def phidotdotsolver_no_rad(phi,phidot,H):
    phidotdot = -dV(phi) - 3*H*phidot
    return phidotdot
       
def Hdotphidotdotsolver_array_no_rad(phi,phidot,H):
    Hdot = np.zeros_like(phi)
    phidotdot = np.zeros_like(phi)
    for i in range(0,np.size(phi)):
        Hdot[i] = Hdotsolver_no_rad(phidot[i])
        phidotdot[i] = phidotdotsolver_no_rad(phi[i],phidot[i],H[i])
    S = [Hdot,phidotdot]
    return S

def Inflation_no_rad(y,t):
    #y[0] = phi
    #y[1] = phidot
    #y[2] = H
    phi = y[0]
    phidot = y[1]
    H = y[2]
    Hdot = Hdotsolver_no_rad(phidot)
    phidotdot = phidotdotsolver_no_rad(phi,phidot,H)
    dy0 = phidot
    dy1 = phidotdot
    dy2 = Hdot
    dy = np.array([dy0,dy1,dy2])
    return dy

#------------------------------------------------------------------------------
# Integration
#------------------------------------------------------------------------------
def Integration_no_rad(ICs_no_rad,N_star):
    phi0 = ICs_no_rad[0]
    print "phi0/mpl is: ", ICs_no_rad[0]/mpl
    print "phidot0 is: ", ICs_no_rad[1]
    print "H0 is: ", ICs_no_rad[2]  
    print "Starting first integration with no radiation"
#    try:
#        f = open('C:\Users\Charlotte\Documents\Lancaster 2\Instant Preheating in Quintessential Inflation\IP_at_zero\Python_Data_Files/QI_IP_no_rad_n_'+str(n)+'_a_'+str(alpha)+'_phi0_'+str(phi0)+'_N_'+str(N_star)+'.lst','r')
#        [phi_no_rad,phidot_no_rad,phidotdot_no_rad,H_no_rad,Hdot_no_rad] = pickle.load(f)
#        print "pickle worked "
#        f.close()
#        return [phi_no_rad,phidot_no_rad,phidotdot_no_rad,H_no_rad,Hdot_no_rad]
#    except: 
    print "no pickle"
    SOL = integrate.odeint(Inflation_no_rad,ICs_no_rad,T1,rtol=tolerance,atol=tolerance)    
    phi_no_rad = SOL[:,0]
    phidot_no_rad = SOL[:,1]
    H_no_rad = SOL[:,2]
    [Hdot_no_rad,phidotdot_no_rad] = Hdotphidotdotsolver_array_no_rad(phi_no_rad,phidot_no_rad,H_no_rad)
#    f = open('C:\Users\Charlotte\Documents\Lancaster 2\Instant Preheating in Quintessential Inflation\IP_at_zero\Python_Data_Files/QI_IP_no_rad_n_'+str(n)+'_a_'+str(alpha)+'_phi0_'+str(phi0)+'_N_'+str(N_star)+'.lst','w')
#    pickle.dump([phi_no_rad,phidot_no_rad,phidotdot_no_rad,H_no_rad,Hdot_no_rad], f)
#    f.close()
    return [phi_no_rad,phidot_no_rad,phidotdot_no_rad,H_no_rad,Hdot_no_rad]

def find_phi_e(ep_minus_1):
    phi_end_comp = fsolve(ep_minus_1,(-0.1*mpl))[0]
    phi_end_analytic = (sp.sqrt(6.0*alpha)*mpl/2.0)*sp.ln((sp.sqrt(3.0*alpha))/(2.0*n))
    #comparing_accuracy(phi_end_comp,phi_end_analytic)    
    return([phi_end_comp,phi_end_analytic])

def get_non_canon_variables(phi_no_rad,phi_end,phidot_no_rad,phidotdot_no_rad):
    phi_non_canon_no_rad = np.sqrt(6*alpha)*mpl*np.tanh(phi_no_rad/(np.sqrt(6*alpha)*mpl))
    phi_end_non_canon = np.sqrt(6*alpha)*mpl*np.tanh(phi_end[0]/(np.sqrt(6*alpha)*mpl))
    phidot_non_canon_no_rad = phidot_no_rad/((np.cosh(phi_no_rad/(np.sqrt(6*alpha)*mpl)))**2)
    phidotdot_non_canon_no_rad = (phidotdot_no_rad/((np.cosh(phi_no_rad/(np.sqrt(6*alpha)*mpl)))**2))- (((2*phidot_no_rad*phidot_no_rad)/(np.sqrt(6*alpha)*mpl))*(1/((np.cosh(phi_no_rad/(np.sqrt(6*alpha)*mpl)))**2))*(np.tanh(phi_no_rad/(np.sqrt(6*alpha)*mpl))))
    V_fun_non_canon = V0*(sp.exp(-n))*(sp.exp(n*(1.0-(phi/(sp.sqrt(6*alpha)*mpl))))-1.0)
    V_non_canon = sp.lambdify(phi,V_fun_non_canon,"numpy")
    dV_non_canon = sp.lambdify(phi,sp.diff(V_fun_non_canon,phi,1),"numpy")
    return [phi_non_canon_no_rad,phi_end_non_canon,phidot_non_canon_no_rad,V_non_canon,dV_non_canon,phidotdot_non_canon_no_rad]

def get_IP_variables(g,nu,phi_vals,phidot_vals,phidotdotvals,H_vals): #IP, will work for either variables
    adi_upper = nu + np.sqrt(phidot_vals/g)
    mod_phi = abs(phi_vals)
    phi_IP = 0
    phidot_IP = 0
    adi_upper_real = 0
    t_IP = 0
    phi_IP_exit = 0
    phidot_IP_exit = 0
    t_IP_exit = 0  
    plt.figure()
    plt.plot(T1,mod_phi,label=r'$\|phi|$')
    plt.plot(T1,adi_upper,"k--")
    plt.legend()
    plt.show()
    for x in range(0,len(phi_vals)):        
            if mod_phi[x] < adi_upper[x] and phi_IP == 0: # first values just within zone
                print "Inside PP region"
                phi_IP = phi_vals[x]
                phidot_IP = phidot_vals[x]
                phidotdot_IP = phidotdotvals[x]
                t_IP = T1[x]
                adi_upper_real = adi_upper[x]
            if t_IP > 0 and mod_phi[x] > adi_upper[x] and T1[x] > t_IP: #values at edge of PP zone for density calcs
                phi_IP_exit = phi_vals[x]
                phidot_IP_exit = phidot_vals[x]
                phidotdotIP_exit = phidotdotvals[x]
                t_IP_exit= T1[x]
                H_IP = H_vals[x]
                print "found exit, index = ", x     
                break
    print "phi_IP/mpl is ", phi_IP/mpl, " the constraint/mpl at this point is ", adi_upper_real/mpl
    return[phi_IP,phidot_IP,t_IP,phi_IP_exit,phidot_IP_exit,t_IP_exit,phidotdotIP_exit,H_IP]

def get_phi_F_from_n_and_alpha(n,alpha):
    phi_F = (((108.0*(np.log(10.0))) - (2.0*n) + (np.log(2.0*n))) *mpl * (np.sqrt(6.0*alpha)))/2.0
    return phi_F

def get_Omega_from_phi_F(phi_F,phi_IP):
    Omega = sp.exp((2.0/3.0)*(1.0 - (((phi_F - phi_IP)*(np.sqrt(3.0/2.0)))/mpl)))
    return Omega
   
def get_g_from_Omega(OMEGA,phi_IP_val,phidot_IP_val):
    kinetic_non_canon_before_ip = phidot_IP_val*phidot_IP_val*0.5
    V_non_canon_ip = V_non_canon(phi_IP_val)
    g = sp.sqrt((8.0*np.pi*np.pi*np.pi*OMEGA*(V_non_canon_ip + kinetic_non_canon_before_ip))/(phidot_IP_val*phidot_IP_val))
    print "g is ", g    
    return g  

def get_Omega_from_g(g_choice,phi_IP_val,phidot_IP_val):
    kinetic_non_canon_before_ip = phidot_IP_val*phidot_IP_val*0.5
    V_non_canon_ip = V_non_canon(phi_IP_val)
    Omega = (g_choice*g_choice*phidot_IP_val*phidot_IP_val) / ((V_non_canon_ip + kinetic_non_canon_before_ip) * 8.0*np.pi*np.pi*np.pi)      
    return Omega
    
def get_accurate_N(T_reh):
    N = 61.93 + sp.ln((V(phi_end[0])**0.25)/mpl)+(1.0/3.0)*(sp.ln((V(phi_end[0])**0.25)/T_reh))
    return N

def get_T_reh_main(rho_chi,rho_phi):
    T_reh = ((30.0/(np.pi*np.pi*106.75))*((rho_chi*rho_chi*rho_chi)/(rho_phi*rho_phi)))**0.25    
    return T_reh  

def get_infl_obvs(N):
    ns = 1 - (2.0/N)
    r = 12.0*alpha/((N + (np.sqrt(3*alpha)/2.0))**2)
    running = 2.0/(N**2 - 2*N)
    return[ns,r,running]

def get_rho_chi(g):
    rho_chi = (g*g*IP_variables[4]*IP_variables[4])/(8.0*np.pi*np.pi*np.pi)     
    return rho_chi
    
def get_rho_phi_end(IP_variables):
    rho_phi = 0.5*IP_variables[4]*IP_variables[4] + V_non_canon(IP_variables[3])
    return rho_phi
    
def get_rho_phi_IP(IP_variables,rho_chi):
    rho_phi = 0.5*IP_variables[4]*IP_variables[4] +V_non_canon(IP_variables[3]) - rho_chi   
    return rho_phi  
    
def backreaction(phi_non_canon,phidot_non_canon,phidotdot_non_canon,dV_non_canon,H,g):
    LHS = phidotdot_non_canon + 3*H*phidot_non_canon + dV_non_canon(phi_non_canon)
    RHS = - (g**(5.0/2.0))*(phidot_non_canon**(3.0/2.0))/(8.0*np.pi*np.pi*np.pi)
    if abs(LHS) > abs((RHS*10.0)) :
        print "\n Back-reaction minimal, LHS is : ", LHS
        print "RHS is : ", RHS
        return False
    else:
        print "Back-reaction important at this g, LHS is: ", LHS
        print "RHS is : ", RHS
        return True

###########################
###      LIVE CODE
###########################

#### LOOPS FOR N ######    
#alpha_vals = [1.0]#(np.sqrt(2.0/3.0))]
#for alpha in alpha_vals:
#    N_guess = 60.0
#    q = 0
#    results = [[],[],[],[],[],[],[],[]]
#    n_min = 122
#    n_max = 125
#    n = (n_min+n_max)/2
#    g_ideal = 0.002
#    g_guess = 1.0
#    while n < n_max and n > n_min:
#        while q < 5:
#        # and n > n_min:
#            print "\n Loop at beginning, n is: ", n
#            print "alpha is: ", alpha
#            print "N_star is: ", N_guess
#            [M,Vfun,V,dV,ddV,ep_minus_1,phi_star,integral] = Get_Equations(n,N_guess)
#            phi_end = find_phi_e(ep_minus_1)
#            print "phi_*/mpl canonical is: ", phi_star/mpl
#            print "phi_end/mpl canonical is at : ", phi_end[0]/mpl
#            print "M/mpl^4 is : ", M/(mpl**4)
#            print "Beginning iteration in canonical variables: "
#            [phi_no_rad,phidot_no_rad,phidotdot_no_rad,H_no_rad,Hdot_no_rad] = Integration_no_rad(Get_IC_no_rad(phi_star),N_guess)
#            
#            V0 = M**4*(sp.exp(-n))
#            print "Switching to non-canonical variables"
#            print "V0 is ", V0
#            [phi_non_canon_no_rad,phi_end_non_canon,phidot_non_canon_no_rad,V_non_canon,dV_non_canon,phidotdot_non_canon_no_rad] = get_non_canon_variables(phi_no_rad,phi_end,phidot_no_rad,phidotdot_no_rad)
#            print "We have the non_canonical variables, we can now look at instant preheating."
#            print "\n Instant Preheating in non-canonical variables: "
#            
#            IP_variables = get_IP_variables(g_guess,0.0,phi_non_canon_no_rad,phidot_non_canon_no_rad,phidotdot_non_canon_no_rad,H_no_rad)
#            print IP_variables
#            print " from n and alpha, we can calculate the distance the field will roll to in the canonical values" 
#            phi_F = get_phi_F_from_n_and_alpha(n,alpha)
#            print "phi_F/mpl = " , phi_F/mpl
#            print "from this we can deduce Omega: ", 
#            Omega_test = get_Omega_from_phi_F(phi_F,IP_variables[3])
#            print "Omega = ", Omega_test
#            print "This omega is from canonical-values, but it should be the same for non-canonical"
#            
#            print 'Entering iteration loop for g'
#            p=0
#            
#            while p < 5:
#                print "\n p = ", p
#                print "For g= ", g_guess, ", instant preheating in non-canonical variables: "
#                IP_variables_guess = get_IP_variables(g_guess,0.0,phi_non_canon_no_rad,phidot_non_canon_no_rad,phidotdot_non_canon_no_rad,H_no_rad)
#                print IP_variables_guess
#                print "we use these IP values to find g from Omega and obtain:"
#                g_val = float(get_g_from_Omega(Omega_test,IP_variables_guess[3],IP_variables_guess[4]))
#                #print g_val
#                if g_val > 100:
#                    print "g was really big"        
#                    break
#                if (g_guess-0.0001) <= g_val and g_val <= (g_guess+0.0001):
#                    print "g iteration converged in ", (p+1), " iterations, g = ", g_val
#                    IP_variables = IP_variables_guess
#                    print "IP_variables are: ", IP_variables                
#                    break
#                else:
#                    print "g values were not close enough, recalculating IP with this new g"        
#                    g_guess = g_val      
#                    print "g_guess is now, " , g_guess
#                    p+=1     
#                        
#            print "from this we can deduce Omega: ", 
#            Omega_new = get_Omega_from_phi_F(phi_F,IP_variables[3])
#            
#            print "Omega = ", Omega_new
#            print "This should match original Omega: ", Omega_test
#            print "from original Omega with these IP values we can find g"
#            g_val = get_g_from_Omega(Omega_test,IP_variables[3],IP_variables[4])
#            print g_val
#            print "To check, this g value gives Omega = ", get_Omega_from_g(g_val,IP_variables[3],IP_variables[4])        
#            
#            T_reh = get_T_reh_main(get_rho_chi(g_val),get_rho_phi_IP(IP_variables,get_rho_chi(g_val)))#Omega_new,((0.5*IP_variables[4]*IP_variables[4])+V_non_canon(IP_variables[3])))
#            N_real = get_accurate_N(get_T_reh_main(get_rho_chi(g_val),get_rho_phi_IP(IP_variables,get_rho_chi(g_val))))# Uses mpl=100 vals as ratios so doesnt matter #Omega_new,((0.5*IP_variables[4]*IP_variables[4])+V_non_canon(IP_variables[3]))))
#            print "This was calculated using N = ", N_guess
#            print "Calculating T_reh, main eq, from the outputs we find : ", T_reh
#            print "In dimensionless units this is T_reh / mpl: ", T_reh/mpl
#            print "Alternatively, in 'real' numbers, our T_reh / mpl(=100) * mpl(=2.14*10**18) = ", T_reh / mpl * (2.43e18)   ,"GeV"     
#            print "This gives N = ", N_real
#            if N_real >= N_guess - 0.01 and N_real <= N_guess + 0.01:
#                print "This N value matches. Inflationary observables are: "
#                print get_infl_obvs(N_guess)
#                print "M is ", M/mpl *2.43e18, " GeV"                
#                break
#            else:
#                N_guess = N_real
#                q += 1        
#####     For finding min n:  
##        min_allowed_g = (np.sqrt(8.0/9.0)*(V_non_canon(phi_end_non_canon)**0.5)*20.0*np.pi)/(mpl*mpl)
##        print "\n min allowed g is: ", min_allowed_g 
##        print " g is ", g_val
##        if g_val > min_allowed_g:
##            print "\n n = ", n, " alpha = ", alpha, " is allowed."
##            n_max = n
##            n = (n_min+n_max)/2
##            print "n_max is ", n_max
##            print "n is ", n
##        else: 
##            print "\n n = ", n, " alpha = ", alpha, " is NOT allowed."
##            n_min = n
##            n = (n_min+n_max)/2
##            print "n_min is ", n_min
##            print "n is ", n
##    print "\n for alpha = ", alpha, " min allowed n is ", n
####        
####   ###  For finding max n using g<1:        
#        if g_val < 1.0:
#            br = backreaction(IP_variables[3],IP_variables[4],IP_variables[6],dV_non_canon,IP_variables[7],g_val)
#
#            print "\n n = ", n, " alpha = ", alpha, " is allowed."
#            n_min = n
#            n = (n_min+n_max)/2
#            print "n_min is ", n_min
#            print "n is ", n
#        else: 
#            br = backreaction(IP_variables[3],IP_variables[4],IP_variables[6],dV_non_canon,IP_variables[7],g_val)
#            print "\n n = ", n, " alpha = ", alpha, " is NOT allowed."
#            n_max = n
#            n = (n_min+n_max)/2
#            print "n_max is ", n_max
#            print "n is ", n
#    print "\n for alpha = ", alpha, " max allowed n is ", n
### For finding n for specific g 
##            min_allowed_g = np.sqrt((32.0*V_non_canon(phi_end_non_canon))/((0.09*np.pi*mpl**4)+(8.0*V_non_canon(phi_end_non_canon))))
##            print " min allowed g is: ", min_allowed_g 
##            new_min_allowed_g = np.sqrt(0.88190346*(V_non_canon(phi_end_non_canon))/(mpl**4))
##            print "new min allowed g is : ", new_min_allowed_g       
##            new_min_allowed_g_no_pi = np.sqrt(0.88190346*np.pi*(V_non_canon(phi_end_non_canon))/(mpl**4))
##            print "new min allowed g no pi is : ", new_min_allowed_g_no_pi       
##    
##            if min_allowed_g < g_ideal:
##                if g_val > g_ideal:
##                    print "\n n = ", n, " alpha = ", alpha, " makes g too big."
##                    n_max = n
##                    n = (n_min+n_max)/2
##                    print "n_max is ", n_max
##                    print "n is ", n
##                else: 
##                    print "\n n = ", n, " alpha = ", alpha, " makes g too small."
##                    n_min = n
##                    n = (n_min+n_max)/2
##                    print "n_min is ", n_min
##                    print "n is ", n 
##            else:
##                print "the minimum allowed value of g is lower than the g you want to find"
##                break
##        print "for alpha = ", alpha, " for you choice of g, ", g_ideal, " n is ", n
##       #

### cycle through n
#        ns = get_infl_obvs(N_guess)[0]
#        r = get_infl_obvs(N_guess)[1]
#        br = backreaction(IP_variables[3],IP_variables[4],IP_variables[6],dV_non_canon,IP_variables[7],g_val)
#        
#        filename = 'C:\Users\Charlotte\Documents\Lancaster 2\Instant Preheating in Quintessential Inflation\IP_at_zero\Python_Data_Files/results_a_003.lst'
#        print "alpha, n, g, T_reh, N*, ns, r, backreation?"
#        print results
#        if os.path.exists(filename):
#            f = open(filename,'rb')
#            print "file opened, results are: "
#            [alpha_old,n_old,g_old,T_old,N_old,ns_old,r_old,br_old] = pickle.load(f)
#            print [alpha_old,n_old,g_old,T_old,N_old,ns_old,r_old,br_old]
#            f.close()    
#            alpha_vals = np.append(alpha_old,alpha)
#            n_vals = np.append(n_old,n)
#            g_vals = np.append(g_old,g_val)
#            T_vals = np.append(T_old,T_reh)
#            N_vals = np.append(N_old,N_real)
#            ns_vals = np.append(ns_old,ns)
#            r_vals = np.append(r_old,r)
#            br_vals = np.append(br_old,br)            
#            f = open(filename,'wb')
#            pickle.dump([alpha_vals,n_vals,g_vals,T_vals,N_vals,ns_vals,r_vals],f)
#            f.close()
#        else:
#            f = open(filename,'wb')
#            pickle.dump([alpha,n,g_val,T_reh,N_real,ns,r,br],f)
#            f.close()
#        
#        n += 1
#
#
###
###    
## For finding max n using RD constraint
##        rho_chi = (g_val**2 * IP_variables_guess[4]**2) / (8.0*np.pi*np.pi*np.pi)
##        rho_phi_before = 0.5*IP_variables_guess[4]**2 + V_non_canon(IP_variables_guess[3])
##        if rho_chi < rho_phi_before - 2*(V_non_canon(IP_variables_guess[3])):
##            print "\n n = ", n, " alpha = ", alpha, " is allowed."
##            n_min = n
##            n = (n_min+n_max)/2
##            print "n_min is ", n_min
##            print "n is ", n
##        else: 
##            print "\n n = ", n, " alpha = ", alpha, " is NOT allowed."
##            n_max = n
##            n = (n_min+n_max)/2
##            print "n_max is ", n_max
##            print "n is ", n
#   
##    print "\n for alpha = ", alpha, " min allowed n is ", n
##
#    
##### For finding max n:        
##        if get_phi_F_from_n_and_alpha(n,alpha) > phi_end[0]:
##            print "\n n = ", n, " alpha = ", alpha, " is allowed."
##            n_min = n
##            n = (n_min+n_max)/2
##            print "n_min is ", n_min
##            print "n is ", n
##        else: 
##            print "\n n = ", n, " alpha = ", alpha, " is NOT allowed."
##            n_max = n
##            n = (n_min+n_max)/2
##            print "n_max is ", n_max
##            print "n is ", n
#
##
##



#### NO LOOP RUN ##########
n = 126.0
alpha = 4.2
N_guess = 60.0
q = 0
g_guess = 1.0
while q < 5:
    print "\n Loop at begining, n is: ", n
    print "alpha is: ", alpha
    print "N_star is: ", N_guess
    [M,Vfun,V,dV,ddV,ep_minus_1,phi_star,integral] = Get_Equations(n,N_guess)
    phi_end = find_phi_e(ep_minus_1)
    print "phi_*/mpl canonical is: ", phi_star/mpl
    print "phi_end/mpl canonical is at : ", phi_end[0]/mpl
    print "M/mpl^4 is : ", M/(mpl**4)
    print "Beginning iteration in canonical variables: "
    [phi_no_rad,phidot_no_rad,phidotdot_no_rad,H_no_rad,Hdot_no_rad] = Integration_no_rad(Get_IC_no_rad(phi_star),N_guess)
    
    V0 = M**4*(sp.exp(-n))
    print "Switching to non-canonical variables"
    print "V0 is ", V0
    [phi_non_canon_no_rad,phi_end_non_canon,phidot_non_canon_no_rad,V_non_canon,dV_non_canon,phidotdot_non_canon_no_rad] = get_non_canon_variables(phi_no_rad,phi_end,phidot_no_rad,phidotdot_no_rad)
    print "We have the non_canonical variables, phi_end/mpl non-canon = ", phi_end_non_canon/mpl, ", we can now look at instant preheating."
    print "We have the non_canonical variables, we can now look at instant preheating."
    print "\n Instant Preheating in non-canonical variables: "
#    kin_for_graph = (phidot_non_canon_no_rad*phidot_non_canon_no_rad*0.5)
#    plt.plot(T1,kin_for_graph,label=r'kin')
#    plt.plot(T1,V_non_canon(phi_non_canon_no_rad),label=r'pot')
#    plt.legend()
    
    IP_variables = get_IP_variables(g_guess,0.0,phi_non_canon_no_rad,phidot_non_canon_no_rad,phidotdot_non_canon_no_rad,H_no_rad)
    print IP_variables    
    print "phi_IP /mpl in canonical : ", np.sqrt(6*alpha)*mpl*np.arctanh(IP_variables[3]/(np.sqrt(6*alpha)*mpl)) /mpl
    print " from n and alpha, we can calculate the distance the field will roll to in the canonical values" 
    phi_F = get_phi_F_from_n_and_alpha(n,alpha)

    print "phi_F/mpl = " , phi_F/mpl
    print "from this we can deduce Omega: ", 
    Omega_test = get_Omega_from_phi_F(phi_F,IP_variables[3])
    print "Omega = ", Omega_test
    print "This omega is from canonical-values, but it should be the same for non-canonical"

    print '\n Entering iteration loop for g'
    p=0
    g_guess = 0.001
    nu=0.0
    while p < 5:
        print "\n p = ", p
        print "For g= ", g_guess, ", instant preheating in non-canonical variables: "
        IP_variables_guess = get_IP_variables(g_guess,nu,phi_non_canon_no_rad,phidot_non_canon_no_rad,phidotdot_non_canon_no_rad,H_no_rad)
        print IP_variables_guess
        print "we use these IP values to find g from Omega and obtain:"
        g_val = float(get_g_from_Omega(Omega_test,IP_variables_guess[3],IP_variables_guess[4]))
        print g_val
        print "Type of g is ", type(g_val)
        if g_val > 100:
            print "g was really big"        
            break
        if (g_guess-0.0001) <= g_val and g_val <= (g_guess+0.0001):
            print "g iteration converged in ", (p+1), " iterations."
            break
        else:
            print "g values were not close enough, recalculating IP with this new g"        
            g_guess = g_val      
            print "g_guess is now, " , g_guess
            p+=1     
    
    print "Obviously this g is less than 1, but what about RD constraints? "
    rho_chi = (g_val**2 * IP_variables_guess[4]**2) / (8.0*np.pi*np.pi*np.pi)
    rho_phi_before = 0.5*IP_variables_guess[4]**2 + V_non_canon(IP_variables_guess[3])
    print "rho_chi = ", rho_chi
    print "This must be less than: ", rho_phi_before - 2*(V_non_canon(IP_variables_guess[3]))
    print "and it is!"
    
    print "\n for the lower g constraint: "
    min_allowed_g = np.sqrt((32.0*V_non_canon(phi_end_non_canon))/((0.09*np.pi*mpl**4)+(8.0*V_non_canon(phi_end_non_canon))))
    print "the min allowed g is ", min_allowed_g
    new_min_allowed_g = (np.sqrt(8.0/9.0)*(V_non_canon(phi_end_non_canon)**0.5)*20.0*np.pi)/(mpl*mpl)
    print "new min allowed g is ", new_min_allowed_g
    #min_allowed_g_full = np.sqrt((32.0*V_non_canon(phi_end_non_canon)*((IP_variables_guess[4]**2)+V_non_canon(IP_variables_guess[3])))/(((0.09*np.pi*mpl**4)+(8.0*V_non_canon(phi_end_non_canon)))*(IP_variables_guess[4]**2)))
    #print "to test assumption in notes, full g_min_calc gives: ", min_allowed_g_full

    T_reh = get_T_reh_main(get_rho_chi(g_val),get_rho_phi_IP(IP_variables,get_rho_chi(g_val)))#Omega_new,((0.5*IP_variables[4]*IP_variables[4])+V_non_canon(IP_variables[3])))
    N_real = get_accurate_N(get_T_reh_main(get_rho_chi(g_val),get_rho_phi_IP(IP_variables,get_rho_chi(g_val))))# Uses mpl=100 vals as ratios so doesnt matter #Omega_new,((0.5*IP_variables[4]*IP_variables[4])+V_non_canon(IP_variables[3]))))
    print "This was calculated using N = ", N_guess
    print "Calculating T_reh, main eq, from the outputs we find : ", T_reh
    print "In dimensionless units this is T_reh / mpl: ", T_reh/mpl
    print "Alternatively, in 'real' numbers, our T_reh / mpl(=100) * mpl(=2.14*10**18) = ", T_reh / mpl * (2.43e18)   ,"GeV"     
    print "This gives N = ", N_real
    print "\n n is ", n
    print "g is ", g_val
    backreaction(IP_variables[3],IP_variables[4],IP_variables[6],dV_non_canon,IP_variables[7],g_val)    
    print phi_non_canon_no_rad
    if N_real >= N_guess - 0.01 and N_real <= N_guess + 0.01:
        print "This N value matches. Inflationary observables are: "
        print get_infl_obvs(N_guess)
        print "M is ", M/mpl * 2.43e18, " GeV"
        break
    else:
        N_guess = N_real
        q += 1    
#    f = open('C:\Users\Charlotte\Documents\Lancaster 2\Instant Preheating in Quintessential Inflation\IP_at_zero\Python_Data_Files/QI_IP_no_rad_n_100_alpha_003_example_for_plotting_phidot_graph.lst','w')
#    pickle.dump([phi_no_rad,phidot_no_rad,phi_non_canon_no_rad,phidot_non_canon_no_rad], f)
#    f.close()    
#    
#    maximum5 = max(phidot_non_canon_no_rad)
#    fig = plt.figure(figsize=(15,15))
#    plt.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.10)
#    ax = fig.add_subplot(1, 1, 1)
#    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
#    ax.spines['right'].set_position('center')
#    #ax.spines['bottom'].set_position('center')
#    # Eliminate upper and right axes
#    ax.spines['left'].set_color('none')
#    ax.spines['top'].set_color('none')
#    # Show ticks in the left and lower axes only
#    ax.xaxis.set_ticks_position('bottom')
#    ax.set_yticks([1.0,1.1])
#    ax.set_yticklabels(['1.0',r'$\frac{\dot{\phi}}{\dot{\phi}_{max}}$'])
#    ax.yaxis.set_ticks_position('right')
#    ax.yaxis.set_label_position('right')
#    plt.plot(phi_no_rad/mpl,phidot_non_canon_no_rad/maximum,linewidth=2.0)
#    plt.xlabel(r'$\frac{\varphi}{m_P}$',labelpad=15)
#    #plt.ylabel(r'$\frac{\dot{\phi}}{\dot{\phi}_{max}}$',labelpad=50)
#    plt.xlim(-4,4)
#    plt.ylim(0,1.1)
#    plt.savefig('phidot_graph_2.pdf')
    #plot_IP_variables_on_graph(IP_variables_guess,g_val,nu)   
          
    #plot_V(phi_non_canon_no_rad/mpl,V_non_canon(phi_non_canon_no_rad),'$V(\phi)$',phidot_non_canon_no_rad,0.1,'Adi constraint')   
    
    #plot_V_t(T1,phi_non_canon_no_rad,V_non_canon(phi_non_canon_no_rad),'$\phi$',phidot_non_canon_no_rad,0.1,'$V(\phi)$','Adi constraint')
    #Plot_V_vs_phidot()
    
    #plot_canon_non_canon_phi()
