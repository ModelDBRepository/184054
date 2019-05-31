COMMENT

Author: Ching-Lung Hsu (adapted from Hines and Carnevale, 2000)

Calcium ion handling mechanisms with endogenous buffer, radial & longitudinal diffusion and pump
A mobile calcium indicator OGB-1 is included

ENDCOMMENT

NEURON {
  SUFFIX cdp
  USEION ca READ cao, cai, ica WRITE cai, ica
  RANGE ica_pmp
  GLOBAL vrat, TotalBuffer, TotalPump 
  GLOBAL TotalIndicator
    : vrat must be GLOBAL--see INITIAL block
    : however TotalBuffer and TotalPump may be RANGE
}

DEFINE Nannuli 4

UNITS {
  (mol)   = (1)
  (molar) = (1/liter)
  (mM)    = (millimolar)
  (um)    = (micron)
  (mA)    = (milliamp)
  FARADAY = (faraday)  (10000 coulomb)
  PI      = (pi)       (1)
}

PARAMETER {
  DCa   = 0.6 (um2/ms)
  
  k1buf = 100 (/mM-ms) : Yamada et al. 1989 : forward rate        : Kd = 1/KD = k1/k2
  k2buf = 0.1 (/ms) : backward rate
  TotalBuffer = 0.003  (mM)
  
  Dogb = 0.3 (um2/ms) : smaller than DCa (Schmidt et al., 2003; Cornelisse et al., 2007)
                      : An unexpected syntax error -
                      : Don't use "DIndicator" as the parameter name for indicator's diffusion constant.
                      : Because after cdp.mod is compiled into cdp.c, the C code will automatically use DIndicator
                      : as an array whenever the state variable Indicator is involved in any REACTION
                      : (including the radial diffusion and the calcium-indicator binding reaction);
                      : otherwise, the following error message will show up in the shell:
                      : "cdp.c: In fuction _ode_specl:
                      :  error: subscripted value is neither array nor pointer"
  k1ind = 430 (/mM-ms)      : Schmidt et al. (2003)
  k2ind = 0.14 (/ms)        : Schmidt et al. (2003)
  TotalIndicator = 0.05 (mM) : take 50% of the pipette concentration as the effective concentration at equilibrium 

  k1    = 1       (/mM-ms)
  k2    = 0.005   (/ms)
  k3    = 1       (/ms)
  k4    = 0.005   (/mM-ms)
  : to eliminate pump, set TotalPump to 0 in hoc
  TotalPump = 1e-11 (mol/cm2)
}

ASSIGNED {
  diam      (um)
  ica       (mA/cm2)
  ica_pmp   (mA/cm2)
  ica_pmp_last   (mA/cm2)
  parea     (um)     : pump area per unit length
  cai       (mM)
  cao       (mM)
  vrat[Nannuli]  (1) : dimensionless
                     : numeric value of vrat[i] equals the volume 
                     : of annulus i of a 1um diameter cylinder
                     : multiply by diam^2 to get volume per um length
  Kd        (/mM)
  B0        (mM)
  Kd_ind    (/mM)    : Ching-Lung
  B0_ind    (mM)     :
}

CONSTANT { volo = 1e10 (um2) }

STATE {
  : ca[0] is equivalent to cai
  : ca[] are very small, so specify absolute tolerance
  : let it be ~1.5 - 2 orders of magnitude smaller than baseline level
  ca[Nannuli]       (mM) <1e-7>
  CaBuffer[Nannuli] (mM) <1e-5>
  Buffer[Nannuli]   (mM) <1e-5>

  CaIndicator[Nannuli] (mM) <1e-5> :
  Indicator[Nannuli]   (mM) <1e-5> :
  
  pump              (mol/cm2) <1e-15>
  pumpca            (mol/cm2) <1e-15>
}

BREAKPOINT {
  SOLVE state METHOD sparse
  ica_pmp_last = ica_pmp
  ica = ica_pmp
}

LOCAL factors_done

INITIAL {
   if (factors_done == 0) {  : flag becomes 1 in the first segment
      factors_done = 1       :   all subsequent segments will have
      factors()              :   vrat = 0 unless vrat is GLOBAL
   }

  Kd = k1buf/k2buf
  B0 = TotalBuffer/(1 + Kd*cai)
  
  Kd_ind = k1ind/k2ind                      
  B0_ind = TotalIndicator/(1 + Kd_ind*cai)

  FROM i=0 TO Nannuli-1 {
    ca[i] = cai
    Buffer[i] = B0
    CaBuffer[i] = TotalBuffer - B0
    Indicator[i] = B0_ind                       
    CaIndicator[i] = TotalIndicator - B0_ind
  }

  parea = PI*diam

: Manually computed initalization of pump
: assumes that cai has been equal to cai0_ca_ion for a long time
  pump = TotalPump/(1 + (cai*k1/k2))   : a better initialization
  pumpca = TotalPump - pump            
: If possible, instead of using formulas to calculate pump and pumpca,
: let NEURON figure them out--just uncomment the following four statements
  ica=0
  ica_pmp = 0
  ica_pmp_last = 0
  SOLVE state STEADYSTATE sparse
: This requires that pump and pumpca be constrained by the CONSERVE
: statement in the STATE block.

}
COMMENT
As suggested by Ted Carnevale,  
the initialization style below may work best

  pump = TotalPump/(1 + (cai*k1/k2))
  pumpca = TotalPump - pump
  SOLVE state STEADYSTATE sparse
ENDCOMMENT


LOCAL frat[Nannuli]  : scales the rate constants for model geometry
PROCEDURE factors() {
  LOCAL r, dr2
  r = 1/2                : starts at edge (half diam)
  dr2 = r/(Nannuli-1)/2  : full thickness of outermost annulus,
                         : half thickness of all other annuli
  vrat[0] = 0
  frat[0] = 2*r
  FROM i=0 TO Nannuli-2 {
    vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  : interior half
    r = r - dr2
    frat[i+1] = 2*PI*r/(2*dr2)  : outer radius of annulus
                                : div by distance between centers
    r = r - dr2
    vrat[i+1] = PI*(r+dr2/2)*2*dr2  : outer half of annulus
  }
}

LOCAL dsq, dsqvol  : can't define local variable in KINETIC block
                   :   or use in COMPARTMENT statement

KINETIC state {
  : COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer}
  COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer CaIndicator Indicator}
  COMPARTMENT (1e10)*parea {pump pumpca}
  COMPARTMENT volo {cao}
  
  LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}        : Longitudinal diffusion of free Ca
                                                              : Longitudinal diffusion of the Indicator,
                                                              : assuming CaIndicator and Indicator have the same mobility
  LONGITUDINAL_DIFFUSION i, Dogb*diam*diam*vrat[i] {CaIndicator Indicator}

  :pump                                                       : a calcium sink
  ~ ca[0] + pump <-> pumpca  (k1*parea*(1e10), k2*parea*(1e10)) : Nannuli = 0 is the outermost annulus
  ~ pumpca <-> pump + cao    (k3*parea*(1e10), k4*parea*(1e10))
  CONSERVE pump + pumpca = TotalPump * parea * (1e10)
  ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

  : all currents except pump
  : ica is Ca efflux
  ~ ca[0] << (-(ica - ica_pmp_last)*PI*diam/(2*FARADAY))       : calcium source
  
  FROM i=0 TO Nannuli-2 {                                      : Radial diffusion of free Ca
    ~ ca[i] <-> ca[i+1]  (DCa*frat[i+1], DCa*frat[i+1])        
    ~ Indicator[i] <-> Indicator[i+1]  (Dogb*frat[i+1], Dogb*frat[i+1])      : Radial diffusion of the Indicator,
    ~ CaIndicator[i] <-> CaIndicator[i+1]  (Dogb*frat[i+1], Dogb*frat[i+1])  : assuming CaIndicator and Indicator have the same mobility
  }
  
  dsq = diam*diam                                              : calcium buffering
  FROM i=0 TO Nannuli-1 {
    dsqvol = dsq*vrat[i]
    ~ ca[i] + Buffer[i] <-> CaBuffer[i]  (k1buf*dsqvol, k2buf*dsqvol)
    ~ ca[i] + Indicator[i] <-> CaIndicator[i]  (k1ind*dsqvol, k2ind*dsqvol)
  }
  
  cai = ca[0]
  
}
