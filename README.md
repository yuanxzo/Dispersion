# **Dispersion**
**This module is used to calculate Rayleigh surface wave dispersion curve and dispersion function value and other related content.**

## Introduction
1. This module contains two forward algorithms:
    - Fast vector-transfer algorithm, corresponding to the prefix of program name: FVTA,
    - Generalized reflection and transmission coefficient algorithm, corresponding prefix: GRTA.

    But only FVTA can be used to get the complete Rayleigh wave dispersion curves.
    
2. This module contains the following subroutines:
    |Number|Subroutine name|input|output|
    |:----:|:--------------|:----|:-----|
    |[1]|FVTA_c(Vr,mods,freq,nf,nv)            |mods, freq, nf and nv         |Vr(nf,nv)|
    |[2]|FVTA_s(Fx,mods,f,v)                   |mods, f and v                 |output: Fx|
    |3  |FVTA_Searchroot(root,mods,f,nroot)    |mods, f and nroot             |root(nroot)|
    |4  |FVTA_bisection(flag,x,a,b,f,mods)     |mods, f and Search range[a,b] |dichotomy search root result, if there is root, flag = 1, root is x, otherwise flag = 0, x = a.|
    |5  |GRTA_s(Fx,mods,f,v)                   |mods, f and v                 |Fx(0:mods%ceng-1)|
    |6  |GRTA_e(eigen,mods,f,nroot)            |mods, f and nroot             |eigen|
