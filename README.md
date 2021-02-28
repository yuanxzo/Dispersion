# **Dispersion**
**This module is used to calculate Rayleigh surface wave dispersion curve and dispersion function value and other related content.**

## Introduction
1. This module contains two forward algorithms:
    - Fast vector-transfer algorithm, corresponding to the prefix of program name: FVTA,
    - Generalized reflection and transmission coefficient algorithm, corresponding prefix: GRTA.

    But only FVTA can be used to get the complete Rayleigh wave dispersion curves.
    
2. This module contains the following subroutines:
    |Number|Subroutine name                    |input                         |output            |notes|
    |:----:|:----------------------------------|:-----------------------------|:-----------------|:----|
    |1  |FVTA_c(Vr,mods,freq,nf,nv)            |mods, freq, nf and nv         |Vr(nf,nv)         |This program can be called externally.|
    |2  |FVTA_s(Fx,mods,f,v)                   |mods, f and v                 |Fx|               |This program can be called externally.|
    |3  |FVTA_Searchroot(root,mods,f,nroot)    |mods, f and nroot             |root(nroot)       |     |
    |4  |FVTA_bisection(flag,x,a,b,f,mods)     |mods, f and Search range[a,b] |flag and a        |If there is root, flag = 1, root is x, otherwise flag = 0, x = a.|
    |5  |GRTA_s(Fx,mods,f,v)                   |mods, f and v                 |Fx(0:mods%ceng-1) |This program can be called externally.|
    |6  |GRTA_e(eigen,mods,f,nroot)            |mods, f and nroot             |eigen             |This program can be called externally.|
    |7  |Haskell_s(Fx,mods,f,v)                |mods, f and v                 |Fx                |This program can be called externally.|
    |8  |Crfinder(vs1,vp1)                     |vs1 and vp1                   |rayv              |     |
    |9  |Rayleigh(R,DR,c,v1,v2)                |c and vs1 and vp1             |R and DR          |Rayleigh wave equation. It is a subroutine of Crfinder.|
    |10 |rough_distance(dc,mods,f)             |mods and f                    |dc                |     |
    |11 |fine_distance(dcout,ndc,mods,f,c1,dc) |mods, f, dc and c1            |ndc and dcout(ndc)|     |
    |12 |ncf(mods,f,v)                         |mods, f and v                 |ncf               |     |
    |13~|Other mathematical functions          |                              |                  |     |
