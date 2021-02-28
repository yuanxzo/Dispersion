# **Dispersion**
> **This module is used to calculate Rayleigh wave dispersion curve and dispersion function value and other related content.**

## Introduction
1. This module contains three forward algorithms:
>   - Fast vector-transfer algorithm, corresponding to the prefix of program name: FVTA,
>   - Generalized reflection and transmission coefficient algorithm, corresponding to the prefix of program name: GRTA,
>   - Modified Thomson-Haskell algorithm，corresponding to the prefix of program name：Haskell.
>    
>   But only FVTA can be used to get the complete Rayleigh wave dispersion curves.
    
2. This module contains the following subroutines:
    |Number|Subroutine name                      |input                         |output            |notes|
    |:----:|:------------------------------------|:-----------------------------|:-----------------|:----|
    |1     |FVTA_c(Vr,mods,freq,nf,nv)           |mods, freq, nf and nv         |Vr(nf,nv)         |This program can be called externally.|
    |2     |FVTA_s(Fx,mods,f,v)                  |mods, f and v                 |Fx|               |This program can be called externally.|
    |3     |FVTA_Searchroot(root,mods,f,nroot)   |mods, f and nroot             |root(nroot)       |     |
    |4     |FVTA_bisection(flag,x,a,b,f,mods)    |mods, f and Search range[a,b] |flag and a        |If there is root, flag = 1, root is x, otherwise flag = 0, x = a.|
    |5     |GRTA_s(Fx,mods,f,v)                  |mods, f and v                 |Fx(0:mods%ceng-1) |This program can be called externally.|
    |6     |GRTA_e(eigen,mods,f,nroot)           |mods, f and nroot             |eigen             |This program can be called externally.|
    |7     |Haskell_s(Fx,mods,f,v)               |mods, f and v                 |Fx                |This program can be called externally.|
    |8     |Crfinder(vs1,vp1)                    |vs1 and vp1                   |rayv              |     |
    |9     |Rayleigh(R,DR,c,v1,v2)               |c, v1 and v2                  |R and DR          |Rayleigh wave equation. It is a subroutine of Crfinder.|
    |10    |rough_distance(dc,mods,f)            |mods and f                    |dc                |     |
    |11    |fine_distance(dcout,ndc,mods,f,c1,dc)|mods, f, dc and c1            |ndc and dcout(ndc)|     |
    |12    |ncf(mods,f,v)                        |mods, f and v                 |ncf               |     |
    |13~   |Other mathematical functions         |                              |                  |     |

3. Parameter-name description:
    |Number|Name |Description|
    |:----:|:----|:----------|
    |1     |mods |A structure of layered elastic stratum model, including the number of layers, S-wave velocity, P-wave velocity, density and layer thickness parameters.|
    |2     |Vr   |A matrix of modal phase velocities.|
    |3     |freq |A vector of frequencies to be calculated.|
    |4     |nf   |A scalar number of frequency points included in freq.|
    |5     |nv   |A scalar number of dispersion curves to be calculated.|
    |6     |Fx   |Dispersion function value.|
    |7     |f    |A scalar value of frequency.|
    |8     |v    |A scalar value of phase velocity.|
    |9     |root |Roots of dispersion function at f frequency.|
    |10    |nroot|A scalar number of roots to calculate.|
    |11    |eigen|A matrix of displacement-stress vectors.|
    |12    |vs1  |A scalar value of S-wave velocity of the first layer of stratum model.|
    |13    |vp1  |A scalar value of P-wave velocity of the first layer of stratum model.|
    |14    |rayv |Rayleigh wave velocity of the first layer in high frequency approximation, rayv is the global variable.|
    |15    |dc   |Search interval of rough search.|
    |16    |c1   |Starting point of root searching.|
    |17    |ndc  |A scalar number of interval points in dc.|
    |18    |dcout|Distance between fine search interval points.|
    |19    |ncf  |Prediction number of roots at frequency(f) and phase velocity(v).|

## Usage
> 1. Intel(R) Visual Fortran is recommended.
> 2. Other compilers can be used, but the subroutine GRTA_e may not be compiled, you can consider commenting or modifying the subroutine's write statement to file "Eigen.txt" to make the compilation successful.
> 3. "Call FVTA_c(Vr,mods,freq,nf,nv)" is a Fortran program statement that calls the module to obtain the model dispersion curve. Other usages are similar. But before that, you must declare each variable in accordance with the definition of the module.

## Licence
> Dispersion is licensed under the GPL v3.0. See the [GNU General Public License](https://github.com/yuanxzo/Dispersion/blob/main/LICENSE) for more details.
