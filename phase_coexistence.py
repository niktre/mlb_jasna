#
#  Maxwell construction for phase coexistence
#
#  Copyright (c) 2015 Ulf D. Schiller <mail@ulfschiller.de>
#  All rights reserved.


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plot
from scipy import integrate, optimize, misc


titlefontsize=24
labelfontsize=20
tickfontsize=12
textfontsize=16
linewidth=2

mpl.rc('font', size=titlefontsize)
mpl.rc('axes', labelsize=labelfontsize)
mpl.rc(('xtick','ytick'), labelsize=tickfontsize)


R1 = 0.5
R2 = 1.0
R3 = 2.0*R2-R1
S1 = 0.36
S2 = 0.44
W  = 0.25*(R2-R1)

Tmin = 0.43
Tmax = 0.47

npoints = 100
delta_rho = (R3-R1)/npoints
epsilon   = 1.e-2


def cs2(x, S2, S1=S1, W=W):
    cs2 = (S2-S1)*np.exp(-(x-R2)*(x-R2)/(2.0*W*W))+S1
    return cs2


def pressure(x, S2, S1=S1, W=W):
    p = x*cs2(x, S2, S1, W)
    return p


def psi(x, S2, S1=S1, W=W):    
    func = lambda x: np.log(cs2(x, S2, S1, W))
    psi = misc.derivative(func,x,delta_rho)
    return psi


def free_energy(rhomin, rhomax, S2, S1=S1, W=W):
    func = lambda x: pressure(1./x, S2, S1, W)
    res = integrate.quad(func, 1./rhomin, 1./rhomax)
    energy = res[0]
    return energy


def critical_density(S2, S1=S1, W=W):
    func = lambda x: psi(x, S2, S1, W)
    res = optimize.minimize_scalar(func)
    rhocrit = res.x
    return rhocrit


def stability_limits(S2, S1=S1, W=W):
    func = lambda x: psi(x, S2, S1, W) + 1./x
    rhocrit = critical_density(S2, S1, W)
    rhomin = optimize.fsolve(func, rhocrit-epsilon)[0]
    rhomax = optimize.fsolve(func, rhocrit+epsilon)[0]
    return rhomin, rhomax


def densities(p, S2, S1=S1, W=W):
    func = lambda x: pressure(x, S2, S1, W) - p
    rhomin, rhomax = stability_limits(S2, S1, W)
    rho1 = optimize.fsolve(func, rhomin-delta_rho)[0]
    rho2 = optimize.fsolve(func, 0.5*(rhomin+rhomax))[0]
    rho3 = optimize.fsolve(func, rhomax+delta_rho)[0]
    return rho1, rho2, rho3


def maxwell_construction(p, S2, S1=S1, W=W):
    rho1, rho2, rho3 = np.sort(densities(p, S2, S1, W))
    area = free_energy(rho1,rho3,S2,S1,W) - p * (1./rho3 - 1./rho1)
    return area


if __name__ == '__main__':

    fig = plot.figure(figsize=(8.,8.))
    ax1 = plot.subplot(223)
    ax2 = plot.subplot(224)
    ax3 = plot.subplot(221)
    ax4 = plot.subplot(222)

    rs = np.arange(R1,R3,delta_rho)
    rg = []
    rl = []
    ts = []

    for T in np.arange(Tmin,Tmax,(Tmax-Tmin)/npoints):

        ps = pressure(rs, S2=T) 
        #es = np.array( [ free_energy(R1,rho,S2=T) for rho in rs ] )

        rhocrit = critical_density(S2=T)
        rhomin, rhomax = stability_limits(S2=T)

        pcrit = pressure(rhocrit, S2=T)
        pmin  = pressure(rhomax, S2=T)
        pmax  = pressure(rhomin, S2=T)
    
        a1 = maxwell_construction(pmin, S2=T)
        a2 = maxwell_construction(pmax, S2=T)

        if (a1*a2 < 0.):

            pv = optimize.bisect(maxwell_construction, pmin, pmax, (T,S1,W))

            rhos = np.array(densities(pv,S2=T))
            vps  = [ pressure(x,S2=T) for x in rhos ]
            
            rhogas    = rhos[0]
            rholiquid = rhos[2]

            rg.append(rhogas)
            rl.append(rholiquid)
            ts.append(T)

    ax1.plot(1/rs, ps, lw=linewidth)
    ax1.plot(1/rhos, vps,"o")
    ax1.plot(1/rs, np.zeros_like(rs)+pv)   
    
    ax1.set_xlabel(r'$1/\rho$')
    ax1.set_ylabel(r'$p(\rho)$')

    ax2.plot(rg,ts,".b")
    ax2.plot(rl,ts,".b")

    ax2.set_xlabel(r'$\rho$')
    ax2.set_ylabel(r'$A$')
    
    ax3.plot(rs, cs2(rs, S2=S2), lw=linewidth)
    
    ax3.set_xlabel(r'$\rho$')
    ax3.set_ylabel(r'$c_s^2(\rho)$')

    ax4.plot(rs, psi(rs, S2=S2), lw=linewidth)   
    ax4.plot(rs, 1/rs, "k:", lw=linewidth)
    ax4.plot(rs, -1/rs, "k:", lw=linewidth)

    ax4.set_xlabel(r'$\rho$')
    ax4.set_ylabel(r'$\psi(\rho)$')

    fig.tight_layout()

    fig.savefig('phase_coexistence.png')

    data = zip(ts,rg,rl)

    np.savetxt('phasdiag_data.txt', data)

    plot.show()
