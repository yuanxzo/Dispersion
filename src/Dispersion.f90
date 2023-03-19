!  Copyright (C) 2021  Bo Yang (seism.yang@foxmail.com)
!  
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <https://www.gnu.org/licenses/>


Module Dispersion

!# **Dispersion**
!    This module is used to calculate Rayleigh wave dispersion curve
!       and dispersion function value and other related content.
!
!## Introduction
!1. This module contains three forward algorithms:
!    - Fast vector-transfer algorithm, corresponding to the prefix of program name: FVTA,
!    - Generalized reflection and transmission coefficient algorithm, corresponding to the prefix of program name: GRTA,
!    - Modified Thomson-Haskell algorithm，corresponding to the prefix of program name：Haskell.
!    
!    But only FVTA can be used to get the complete Rayleigh wave dispersion curves.
!    
!2. This module contains the following subroutines:
!    |Number|Subroutine name                      |input                         |output            |notes|
!    |:----:|:------------------------------------|:-----------------------------|:-----------------|:----|
!    |1     |FVTA_c(Vr,mods,freq,nf,nv)           |mods, freq, nf and nv         |Vr(nf,nv)         |This program can be called externally.|
!    |2     |FVTA_s(Fx,mods,f,v)                  |mods, f and v                 |Fx|               |This program can be called externally.|
!    |3     |FVTA_Searchroot(root,mods,f,nroot)   |mods, f and nroot             |root(nroot)       |     |
!    |4     |FVTA_bisection(flag,x,a,b,f,mods)    |mods, f and Search range[a,b] |flag and a        |If there is root, flag = 1, root is x, otherwise flag = 0, x = a.|
!    |5     |GRTA_s(Fx,mods,f,v)                  |mods, f and v                 |Fx(0:mods%ceng-1) |This program can be called externally.|
!    |6     |GRTA_e(eigen,mods,f,nroot)           |mods, f and nroot             |eigen             |This program can be called externally.|
!    |7     |Haskell_s(Fx,mods,f,v)               |mods, f and v                 |Fx                |This program can be called externally.|
!    |8     |Crfinder(vs1,vp1)                    |vs1 and vp1                   |rayv              |     |
!    |9     |Rayleigh(R,DR,c,v1,v2)               |c, v1 and v2                  |R and DR          |Rayleigh wave equation. It is a subroutine of Crfinder.|
!    |10    |rough_distance(dc,mods,f)            |mods and f                    |dc                |     |
!    |11    |fine_distance(dcout,ndc,mods,f,c1,dc)|mods, f, dc and c1            |ndc and dcout(ndc)|     |
!    |12    |ncf(mods,f,v)                        |mods, f and v                 |ncf               |     |
!    |13~   |Other mathematical functions         |                              |                  |     |
!
!3. Parameter-name description:
!    |Number|Name |Description|
!    |:----:|:----|:----------|
!    |1     |mods |A structure of layered elastic stratum model, including the number of layers, S-wave velocity, P-wave velocity, density and layer thickness parameters.|
!    |2     |Vr   |A matrix of modal phase velocities.|
!    |3     |freq |A vector of frequencies to be calculated.|
!    |4     |nf   |A scalar number of frequency points included in freq.|
!    |5     |nv   |A scalar number of dispersion curves to be calculated.|
!    |6     |Fx   |Dispersion function value.|
!    |7     |f    |A scalar value of frequency.|
!    |8     |v    |A scalar value of phase velocity.|
!    |9     |root |Roots of dispersion function at f frequency.|
!    |10    |nroot|A scalar number of roots to calculate.|
!    |11    |eigen|A matrix of displacement-stress vectors.|
!    |12    |vs1  |A scalar value of S-wave velocity of the first layer of stratum model.|
!    |13    |vp1  |A scalar value of P-wave velocity of the first layer of stratum model.|
!    |14    |rayv |Rayleigh wave velocity of the first layer in high frequency approximation, rayv is the global variable.|
!    |15    |dc   |Search interval of rough search.|
!    |16    |c1   |Starting point of root searching.|
!    |17    |ndc  |A scalar number of interval points in dc.|
!    |18    |dcout|Distance between fine search interval points.|
!    |19    |ncf  |Prediction number of roots at frequency(f) and phase velocity(v).|


    Implicit None
    Integer,parameter::layer=100        !The maximum number of layers of the model that this module allows to calculate.User can adjust.
    Integer,parameter::npoint=1000      !The number of discrete points of eigendisplacement and eigentraction curve.User can adjust.
    Real(kind=8),parameter::ednn=0.15d0 !The minimum difference of the prediction formula ncf between the two phase velocities when searching for roots, the user can adjust it to be smaller when some roots are not searched. The ednn value range is (0, 0.25]
    Real(kind=8),parameter::pi=acos(-1.0d0)
    Real(kind=8)::rayv,minvs,maxvs      !The first layer Rayleigh wave velocity, minimum and maximum shear wave velocity of the model
    
    ! Generalized reflection transmission coefficients
    Complex(kind=8)::Rud(2,2,0:layer),Rdu(2,2,layer),Tdd(2,2,layer),Tuu(2,2,layer),eye(2,2)
    Complex(kind=8)::rs(layer),rp(layer),tt(layer),E11(2,2,layer),E12(2,2,layer),E21(2,2,layer),E22(2,2,layer)

    ! The data structure of the layered model defined in the module.
    Type Model
        Integer::ceng !Number of layers
        Real(kind=8)::vs(layer),vp(layer),dns(layer),thk(layer-1) !S-wave velocity, P-wave velocity, density and layer thickness
    End Type
    
    ! det Function：find the determinant of [real number, complex number] square matrix (n * n)
    Interface det
        Module procedure r_det,z_det
    End interface
    ! sqrt Function: Extended function for finding the square root. It can find the square root of integer and negative number.
    Interface sqrt
        Module procedure int_sqrt,neg_sqrt
    End interface
Contains
    
    ! The complete Rayleigh wave dispersion curves is obtained by FVTA
    Subroutine FVTA_c(Vr,mods,freq,nf,nv)
        Type(Model)             :: mods
        Integer,intent(in)      :: nf,nv
        Real(kind=8),intent(in) :: freq(nf)
        Real(kind=8),intent(out):: Vr(nf,nv)
        Integer::i
        
        Vr(1:nf,1:nv)=0.0d0
        Call Crfinder(mods%vs(1),mods%vp(1))
        minvs=minval(mods%vs(1:mods%ceng))
        maxvs=maxval(mods%vs(1:mods%ceng))
        !$OMP PARALLEL DO
        Do i=1,nf
            Call FVTA_Searchroot(Vr(i,1:nv),mods,freq(i),nv)
        End Do
        !$OMP End PARALLEL DO
    End Subroutine FVTA_c
    
    ! Fast vector-transfer algorithm(FVTA)
    Subroutine FVTA_s(Fx,mods,f,v)
        Type(model):: mods
        Real(kind=8),intent(out)::Fx  ! The value of the dispersion function of the FVTA
        Real(kind=8),intent(in) ::f
        
        Real(kind=8)::w,u(mods%ceng),l,ks,v
        Complex(kind=8)::rp,rs,rr,EE(5),FE(5),FF(5,5)
        Complex(kind=8)::ic,g,p,q,a,b,d,e,r,s
        Integer::i
        
        Fx=0.0
        ic=(0.0d0,1.0d0)
        w=2.0*pi*f
        u=(mods%vs(1:mods%ceng)**2.0d0)*mods%dns(1:mods%ceng)
        If (minval(abs(v-mods%vs(1:mods%ceng)))==0.0 .or. minval(abs(v-mods%vp(1:mods%ceng)))==0.0) Then
            v=v-1.0e-9
        End If
        
        rp=sqrt(cmplx(v**2/mods%vp(mods%ceng)**2-1.0))
        rs=sqrt(cmplx(v**2/mods%vs(mods%ceng)**2-1.0))
        rr=1.0-v**2/(2.0*mods%vs(mods%ceng)**2)   
        EE(1)=1.0+rp*rs
        EE(2)=rr+rp*rs
        EE(3)=ic*rs*(1.0-rr)
        EE(4)=ic*rp*(rr-1.0)
        EE(5)=-rr**2-rp*rs
        
        Do i=mods%ceng-1,1,-1
            rp=sqrt(cmplx(v**2/mods%vp(i)**2-1.0))
            rs=sqrt(cmplx(v**2/mods%vs(i)**2-1.0))
            rr=1.0-v**2/(2.0*mods%vs(i)**2)
            l=u(i+1)/u(i)
            g=1.0-rr
            ks=w/v
            p=rp*ks*mods%thk(i)
            q=rs*ks*mods%thk(i)
            a=cos(p)
            b=cos(q)
            
            e=sin(p)/rp
            d=sin(q)/rs
            
            If (abs(a)>10.0d0**40) Then
                a=1.0
                e=-1.0/(ic*rp)
                
                b=1.0
                d=-1.0/(ic*rs)
            End If
            
            r=rp**2
            s=rs**2
            
            FF(1,1)=(a*b*(1.0+rr**2)-e*d*(r*s+rr**2)-2.0*rr)/l
            FF(1,2)=2.0*(1.0-a*b)*(1.0+rr)+2.0*e*d*(r*s+rr)
            FF(1,3)=-(a*d+b*e*r)*g
            FF(1,4)=(a*d*s+b*e)*g
            FF(1,5)=(2.0*(1.0-a*b)+e*d*(r*s+1.0))*l
            
            FF(2,1)=((a*b-1.0)*(rr+rr**2)-e*d*(r*s+rr**3))/l
            FF(2,2)=-4.0*a*b*rr+2.0*e*d*(r*s+rr**2)+(1.0+rr)**2
            FF(2,3)=-(a*d*rr+b*e*r)*g
            FF(2,4)=(a*d*s+b*e*rr)*g
            FF(2,5)=((1.0-a*b)*(1.0+rr)+e*d*(r*s+rr))*l
            
            FF(3,1)=(a*d*s+b*e*rr**2)*g/l
            FF(3,2)=-2.0*(a*d*s+b*e*rr)*g
            FF(3,3)=a*b*g**2
            FF(3,4)=e*d*s*g**2
            FF(3,5)=-(a*d*s+b*e)*g*l
            
            FF(4,1)=-(a*d*rr**2+b*e*r)*g/l
            FF(4,2)=2.0*(a*d*rr+b*e*r)*g
            FF(4,3)=e*d*r*g**2
            FF(4,4)=a*b*g**2
            FF(4,5)=(a*d+b*e*r)*g*l
            
            FF(5,1)=(2.0*(1.0-a*b)*rr**2+e*d*(r*s+rr**4))/l
            FF(5,2)=2.0*(a*b-1.0)*(rr+rr**2)-2.0*e*d*(r*s+rr**3)
            FF(5,3)=(a*d*rr**2+b*e*r)*g
            FF(5,4)=-(a*d*s+b*e*rr**2)*g
            FF(5,5)=(a*b*(1.0+rr**2)-e*d*(r*s+rr**2)-2.0*rr)*l
            
            FE(1:5)=EE(1:5)
            
            EE(1)=FF(1,1)*FE(1)+ FF(1,2)*FE(2)+ FF(1,3)*FE(3)+ FF(1,4)*FE(4)+ FF(1,5)*FE(5)
            EE(2)=FF(2,1)*FE(1)+ FF(2,2)*FE(2)+ FF(2,3)*FE(3)+ FF(2,4)*FE(4)+ FF(2,5)*FE(5)
            EE(3)=FF(3,1)*FE(1)+ FF(3,2)*FE(2)+ FF(3,3)*FE(3)+ FF(3,4)*FE(4)+ FF(3,5)*FE(5)
            EE(4)=FF(4,1)*FE(1)+ FF(4,2)*FE(2)+ FF(4,3)*FE(3)+ FF(4,4)*FE(4)+ FF(4,5)*FE(5)
            EE(5)=FF(5,1)*FE(1)+ FF(5,2)*FE(2)+ FF(5,3)*FE(3)+ FF(5,4)*FE(4)+ FF(5,5)*FE(5)

        End Do
        Fx=Real(EE(5))
    End Subroutine FVTA_s
    
    ! Calculate the 1 to nroot roots of the model mods at a certain frequency f.
    ! The roots sorted by phase velocity from small to large.
    Subroutine FVTA_Searchroot(root,mods,f,nroot)
        Type(Model)        ::mods
        Integer,intent(in) ::nroot
        Real(kind=8),intent(in) ::f
        Real(kind=8),intent(out)::root(nroot)
        
        Real(kind=8)::c1,dc,dcout(20),c2,x
        Integer::t,ndc,flag,i
        
        root(nroot)=0
        c1=0.88*minvs
        t=1
        Call rough_distance(dc,mods,f)
        Call fine_distance(dcout,ndc,mods,f,c1,dc)

        c2=c1+dcout(1)
        Do while (c2<=maxvs)
            Do i=1,ndc
                c2=c1+dcout(i)
                Call FVTA_bisection(flag,x,c1,c2,f,mods)
                If (flag==1 .and. t<=nroot) Then
                    root(t)=x
                    t=t+1
                End If
                c1=c2
                If (t>nroot) Then
                    Return
                End If
            End Do
            If (t>nroot .or. c1==maxvs .or. abs(c1-maxvs)<1.0e-6) Then
                Return
            End If
            Call fine_distance(dcout,ndc,mods,f,c1,dc)
            c2=c1+dcout(1)
        End Do
    End Subroutine FVTA_Searchroot
    
    ! Bisection search root program
    Subroutine FVTA_bisection(flag,x,a,b,f,mods)
        Type(model):: mods
        Real(kind=8)::a,b,f,x,xeps,yeps,fa,fb,fc,dx,a1,a2,c,erro,fa1,fa2
        Integer::flag
        
        xeps=1.0e-6
        yeps=1.0e-10
        dx=yeps*(b-a)

        Call FVTA_s(fa,mods,f,a)
        do while (isnan(fa))
            a=a+dx;Call FVTA_s(fa,mods,f,a)
        end do
        Call FVTA_s(fb,mods,f,b)
        do while (isnan(fb))
            b=b-dx;Call FVTA_s(fb,mods,f,b)
        end do

        flag=0;x=a
        If (fa==0.0d0) Then
            a1=a-dx;a2=a+dx
            Call FVTA_s(fa1,mods,f,a1)
            Call FVTA_s(fa2,mods,f,a2)
            If (fa1*fa2<0.0d0) Then
                flag=1
                x=a
            Else
                flag=0
                x=a
            End If
            Return
        End If
        If (fb==0.0d0) Then
            flag=0
            x=a
            Return
        End If
        If (fa*fb>0.0d0) Then
            flag=0
            x=a
            Return
        End If
        If (fa*fb<0.0d0) Then
            flag=1
            c=(a+b)/2
            Call FVTA_s(fc,mods,f,c)
            do while (isnan(fc))
                c=c+dx;Call FVTA_s(fc,mods,f,c)
            end do
            erro=b-a
            Do while (erro>xeps)
                If (fc*fa<0.0d0) Then
                    b=c
                    fb=fc
                    c=(a+b)/2
                    erro=b-a
                    Call FVTA_s(fc,mods,f,c)
                    do while (isnan(fc))
                        c=c+dx;Call FVTA_s(fc,mods,f,c)
                    end do
                End If
                If (fc*fb<0.0d0) Then
                    a=c
                    fa=fc
                    c=(a+b)/2
                    erro=b-a
                    Call FVTA_s(fc,mods,f,c)
                    do while (isnan(fc))
                        c=c+dx;Call FVTA_s(fc,mods,f,c)
                    end do
                End If
            End Do
	    x=c
        End If
    End Subroutine FVTA_bisection
    
    ! Generalized reflection and transmission coefficient algorithm
    Subroutine GRTA_s(Fx,mods,f,v)
        Type(model)::mods
        Complex(kind=8),intent(out)::Fx(0:mods%ceng-1) !Secular function values for all layers
        Integer::i
        Real(kind=8) ::f,v,k,miu(mods%ceng),w,w1
        Complex(kind=8)::Dad(2,2,mods%ceng-1)
        
        Complex(kind=8)::A(4,4),B(4,4),C(4,4),D(4,4)
        Complex(kind=8)::td(2,2,mods%ceng-1),ru(2,2,mods%ceng-1),rd(2,2,mods%ceng-1),tu(2,2,mods%ceng-1)
        
        eye=(0.0d0,0.0d0);eye(1,1)=(1.0d0,0.0d0);eye(2,2)=(1.0d0,0.0d0);
        
        w=2.0d0*pi*f;  k=w/v;  w1=1.0d0/w
        
        rs(1:mods%ceng)=sqrt(cmplx(k**2-(w/mods%vs(1:mods%ceng))**2.0d0))
        rp(1:mods%ceng)=sqrt(cmplx(k**2-(w/mods%vp(1:mods%ceng))**2.0d0))
        tt(1:mods%ceng)=k**2+rs**2
        miu(1:mods%ceng)=(mods%vs(1:mods%ceng)**2.0d0)*mods%dns(1:mods%ceng);
	!miu=Real(mods%ceng)*miu/sum(miu)    
     
        E11(1,1,1:mods%ceng)=mods%vp(1:mods%ceng)*k*w1
        E11(1,2,1:mods%ceng)=mods%vs(1:mods%ceng)*rs*w1
        E11(2,1,1:mods%ceng)=mods%vp(1:mods%ceng)*rp*w1
        E11(2,2,1:mods%ceng)=mods%vs(1:mods%ceng)*k*w1
        
        E12(1,1,1:mods%ceng)=mods%vp(1:mods%ceng)*k*w1
        E12(1,2,1:mods%ceng)=mods%vs(1:mods%ceng)*rs*w1
        E12(2,1,1:mods%ceng)=-mods%vp(1:mods%ceng)*rp*w1
        E12(2,2,1:mods%ceng)=-mods%vs(1:mods%ceng)*k*w1
            
        E21(1,1,1:mods%ceng)=-2.0d0*mods%vp(1:mods%ceng)*miu*k*rp*w1
        E21(1,2,1:mods%ceng)=-mods%vs(1:mods%ceng)*miu*tt*w1
        E21(2,1,1:mods%ceng)=-mods%vp(1:mods%ceng)*miu*tt*w1
        E21(2,2,1:mods%ceng)=-2.0d0*mods%vs(1:mods%ceng)*miu*k*rs*w1
        
        E22(1,1,1:mods%ceng)=2.0d0*mods%vp(1:mods%ceng)*miu*k*rp*w1
        E22(1,2,1:mods%ceng)=mods%vs(1:mods%ceng)*miu*tt*w1
        E22(2,1,1:mods%ceng)=-mods%vp(1:mods%ceng)*miu*tt*w1
        E22(2,2,1:mods%ceng)=-2.0d0*mods%vs(1:mods%ceng)*miu*k*rs*w1
        
        Dad=(0.0d0,0.0d0)
        Dad(1,1,1:mods%ceng-1)=exp(-rp(1:mods%ceng-1)*mods%thk(1:mods%ceng-1))
        Dad(2,2,1:mods%ceng-1)=exp(-rs(1:mods%ceng-1)*mods%thk(1:mods%ceng-1))
        
        Rud(:,:,0)=matmul(matmul(-inv_2(E21(1:2,1:2,1)),E22(1:2,1:2,1)),Dad(1:2,1:2,1))
        
        td(:,:,1:mods%ceng-1)=(0.0d0,0.0d0);  ru(:,:,1:mods%ceng-1)=(0.0d0,0.0d0)
        rd(:,:,1:mods%ceng-1)=(0.0d0,0.0d0);  tu(:,:,1:mods%ceng-1)=(0.0d0,0.0d0)
        
        Do i=1,mods%ceng-1
            A(1:2,1:2)=E11(:,:,i+1);  A(1:2,3:4)=-E12(:,:,i)
            A(3:4,1:2)=E21(:,:,i+1);  A(3:4,3:4)=-E22(:,:,i)
            A=inv_4(A)
            If (i<mods%ceng-1) Then
                B(1:2,1:2)=E11(:,:,i);  B(1:2,3:4)=-E12(:,:,i+1)
                B(3:4,1:2)=E21(:,:,i);  B(3:4,3:4)=-E22(:,:,i+1)
                    
                C=(0.0d0,0.0d0);  C(1:2,1:2)=Dad(:,:,i);  C(3:4,3:4)=Dad(:,:,i+1)
                D=matmul(matmul(A,B),C)
                td(:,:,i)=D(1:2,1:2);ru(:,:,i)=D(1:2,3:4);
                rd(:,:,i)=D(3:4,1:2);tu(:,:,i)=D(3:4,3:4);
            End If 
            If (i==mods%ceng-1) Then
                C=(0.0d0,0.0d0); D=(0.0d0,0.0d0)
                C(1:2,1:2)=matmul(E11(:,:,mods%ceng-1),Dad(:,:,mods%ceng-1))
                C(3:4,1:2)=matmul(E21(:,:,mods%ceng-1),Dad(:,:,mods%ceng-1))
                    
                D(1:4,1:2)=matmul(A(1:4,1:4),C(1:4,1:2))
                td(:,:,mods%ceng-1)=D(1:2,1:2)
                rd(:,:,mods%ceng-1)=D(3:4,1:2)
            End If 
        End Do
            
        Rud(:,:,1:mods%ceng-1)=(0.0d0,0.0d0);  Rdu(:,:,1:mods%ceng)  =(0.0d0,0.0d0)
        Tdd(:,:,1:mods%ceng-1)=(0.0d0,0.0d0);  Tuu(:,:,1:mods%ceng-1)=(0.0d0,0.0d0)

        Do i=mods%ceng-1,1,-1
            Tdd(:,:,i)=matmul(inv_2(eye-matmul(ru(:,:,i),Rdu(:,:,i+1))),td(:,:,i))
            Rdu(:,:,i)=rd(:,:,i)+matmul(matmul(tu(:,:,i),Rdu(:,:,i+1)),Tdd(:,:,i))
        End Do
        Do i=1,mods%ceng-1
            Tuu(:,:,i)=matmul(inv_2(eye-matmul(rd(:,:,i),Rud(:,:,i-1))),tu(:,:,i))
            Rud(:,:,i)=ru(:,:,i)+matmul(matmul(td(:,:,i),Rud(:,:,i-1)),Tuu(:,:,i))
        End Do
        
        Fx(0)=(det(E21(:,:,1)+matmul(matmul(E22(:,:,1),Dad(:,:,1)),Rdu(:,:,1)),2)) !at surface
	Do i=1,mods%ceng-1
            Fx(i)=(det(eye-matmul(Rud(:,:,i-1),Rdu(:,:,i)),2))
        End Do
    End Subroutine GRTA_s
    
    Subroutine GRTA_e(eigen,mods,f,nroot)
        Type(Model)       ::mods
        Integer,intent(in)::nroot
        Real(kind=8),intent(in)::f
        Real(kind=8),intent(out)::eigen(npoint,4,nroot)
        
        Integer::i,j,k,ii
        Real(kind=8)::w,dz,z(npoint),H(0:mods%ceng),Vr(nroot),freq(1),etemp(5*nroot)
        Complex(kind=8)::Fx(0:mods%ceng-1),a(2,2),b(2,2),Dud(2,2),Cd(2,mods%ceng),Cu(2,mods%ceng)
        Complex(kind=8)::EE(4,4),lmda(4,4),CC(4,1),temp(4,1)
        
        Open(13,file='Eigen.txt')
        w=2.0d0*pi*f
        H(0)=0.0d0
        Do i=1,mods%ceng-1
            H(i)=H(i-1)+mods%thk(i)
        End Do
        H(mods%ceng)=2*sum(mods%thk(1:mods%ceng-1))
        dz=H(mods%ceng)/Real(npoint-1)
        z(1)=0.0d0
        Do i=2,npoint
            z(i)=z(i-1)+dz
        End Do
        Vr(1:nroot)=0.0d0
        freq(1)=f
        Call FVTA_c(Vr,mods,freq,1,nroot)
        ii=1
        Do i=1,nroot
            If (Vr(i)==0.0d0) Then
                ii=i-1
                exit
            End If
        End Do
        Do i=1,ii
            Call GRTA_s(Fx,mods,f,Vr(i))
            a=E21(:,:,1)/w
            Dud(1:2,1:2)=(0.0d0,0.0d0);Dud(1,1)=exp(-rp(1)*mods%thk(1));Dud(2,2)=exp(-rs(1)*mods%thk(1))
            b=matmul(matmul(E22(:,:,1),Dud),Rdu(:,:,1))/w
            Cd(1,1)=a(1,2)+b(1,2); Cd(2,1)=-a(1,1)-b(1,1)
            Cd(:,1)=Cd(:,1)/(sqrt(abs(a(1,1)+b(1,1))**2+abs(a(1,2)+b(1,2))**2))
            Do j=1,mods%ceng-1
                Cu(:,j)=matmul(Rdu(:,:,j),Cd(:,j))                
                Cd(:,j+1)=matmul(Tdd(:,:,j),Cd(:,j))
            End Do
            Do j=1,npoint
                Do k=1,mods%ceng
                    If (z(j)>=H(k-1) .and. z(j)<H(k)) Then
                        exit
                    End If
                End Do
                EE(1:2,1:2)=E11(:,:,k);EE(1:2,3:4)=E12(:,:,k)
                EE(3:4,1:2)=E21(:,:,k);EE(3:4,3:4)=E22(:,:,k)
                lmda(1:4,1:4)=(0.0d0,0.0d0)
                If (k<mods%ceng) Then
                    lmda(1,1)=exp(-rp(k)*(z(j)-H(k-1)))
                    lmda(2,2)=exp(-rs(k)*(z(j)-H(k-1)))
                    lmda(3,3)=exp(-rp(k)*(H(k)-z(j)))
                    lmda(4,4)=exp(-rs(k)*(H(k)-z(j)))
                    CC(1:2,1)=Cd(:,k);CC(3:4,1)=Cu(:,k)
                    temp=matmul(matmul(EE,lmda),CC)
                    eigen(j,1:4,i)=Real(temp(1:4,1))
                Else
                    lmda(1,1)=exp(-rp(mods%ceng)*(z(j)-H(mods%ceng-1)))
                    lmda(2,2)=exp(-rs(mods%ceng)*(z(j)-H(mods%ceng-1)))
                    temp(1:2,1)=matmul(matmul(E11(:,:,mods%ceng),lmda(1:2,1:2)),Cd(:,mods%ceng))
                    temp(3:4,1)=matmul(matmul(E21(:,:,mods%ceng),lmda(1:2,1:2)),Cd(:,mods%ceng))
                    eigen(j,1:4,i)=Real(temp(1:4,1))
                End If
            End Do
        End Do
        Write(*,*) "<Test of Dispersion>  Copyright (C) <2021>  <Bo Yang. Email:yuanxzo@qq.com>"
        Write(*,*) "This program comes with ABSOLUTELY NO WARRANTY."
        Write(*,*) "This is free software, and you are welcome to redistribute it"
        Write(*,*) "under certain conditions; see <https://www.gnu.org/licenses/> for details."
        Write(*,*)
        Write(*,*) "Press any key to confirm..."
        Read(*,*)
        
        Do i=1,ii
            Write(13,'(<npoint>f15.6)') z
            Do j=1,4
                Write(13,'(<npoint>f15.6)') eigen(:,j,i)
            End Do
        End Do
        Write(*,*) "TThe eigendisplacement and eigentraction have been output in the file: Eigen.txt!"

        !! Matlab code for drawing eigenfunction curve:******************************************
        !load('Eigen.txt')
        ![m,n]=size(Eigen);
        !%% Plot the eigendisplacement (normalized in the drawing)
		!%% the solid line is the vertical component, and the dashed line is the horizontal component.
        !figure
        !subplot(1,2,1)
        !for i=1:m/5
        !    pppv=Eigen(5*i-2,:)/max(2*abs(Eigen(5*i-2,:)))+i;
        !    ppph=Eigen(5*i-3,:)/max(2*abs(Eigen(5*i-3,:)))+i;
        !    plot(Eigen(1,:),pppv,'-');hold on;grid on
        !    plot(Eigen(1,:),ppph,'--');hold on
        !End
        !xlabel('Depth（m）');ylabel('Rayleigh modes');title('Normalized eigendisplacement')
        !%% Plot the eigentraction (normalized in the drawing)
		!%% the solid line is the vertical component, and the dashed line is the horizontal component.
        !subplot(1,2,2)
        !for i=1:m/5 
        !    pppv=Eigen(5*i,:)/max(2*abs(Eigen(5*i,:)))+i;
        !    ppph=Eigen(5*i-1,:)/max(2*abs(Eigen(5*i-1,:)))+i;
        !    plot(Eigen(1,:),pppv,'-');hold on;grid on
        !    plot(Eigen(1,:),ppph,'--');hold on
        !End
        !xlabel('Depth（m）');ylabel('Rayleigh modes');title('Normalized eigentraction')
        !! Matlab code for drawing eigenfunction curve:******************************************
        
    End Subroutine GRTA_e
    
    ! Modified Thomson-Haskell algorithm.
    Subroutine Haskell_s(Fx,mods,f,v)
        type(model)        ::mods
        Real(kind=8),intent(out)::Fx
        Real(kind=8),intent(in)::f,v
        Real(kind=8)::miu(mods%ceng),k,RR
        Complex(kind=8)::iu,s,t,r,EE(4,4),MM(4,4),PP(4,4),TT(4,4),VV(4,2),XX(2,4),temp(2,2)
        integer::i

        Fx=0.0d0 ! To initialize the output variable
        XX(1,1:4)=[0.0d0,0.0d0,1.0d0,0.0d0] ! To initialize the boundary matrix determined by free surface
        XX(2,1:4)=[0.0d0,0.0d0,0.0d0,1.0d0]
        iu=(0.0d0,1.0d0) ! To initialize the imaginary unit
        EE=(0.0d0,0.0d0);MM=(0.0d0,0.0d0);PP=(0.0d0,0.0d0);TT=(0.0d0,0.0d0) ! To initialize the propagation matrices

        miu=mods%dns(1:mods%ceng)*mods%vs(1:mods%ceng)**2
        miu=real(mods%ceng)*miu/sum(miu)
        k=2.0d0*pi*f/v

        do i=1,mods%ceng,1
            if (v<mods%vs(i)) then
                s=sqrt(1.0d0-(v/mods%vs(i))**2)
            else
                s=iu*sqrt((v/mods%vs(i))**2-1.0d0)
            end if
            if (v<mods%vp(i)) then
                r=sqrt(1.0d0-(v/mods%vp(i))**2)
            else
                r=iu*sqrt((v/mods%vp(i))**2-1.0d0)
            end if
            t=2.0d0-(v/mods%vs(i))**2
            if (i<mods%ceng) then
                RR=abs(exp(0.5*k*mods%thk(i)*(r+s)))
                !RR=1.0D0

                EE(1,1)=exp(k*r*mods%thk(i))
                EE(2,2)=exp(-k*r*mods%thk(i))
                EE(3,3)=exp(k*s*mods%thk(i))
                EE(4,4)=exp(-k*s*mods%thk(i))
                EE=EE/RR

                MM(1,1)=1.0d0
                MM(2,2)=1.0d0
                MM(3,3)=miu(i)
                MM(4,4)=miu(i)

                PP(1,1:4)=[(1.0d0,0.0d0),(1.0d0,0.0d0),s,-s]
                PP(2,1:4)=[r,-r,(1.0d0,0.0d0),(1.0d0,0.0d0)]
                PP(3,1:4)=[2.0d0*r,-2.0d0*r,t,t]
                PP(4,1:4)=[t,t,2.0d0*s,-2.0d0*s]

                TT=matmul(matmul(matmul(matmul(MM,PP),EE),inv_4(PP)),inv_4(MM))
                XX=matmul(XX,TT)
            end if

        end do
		
	! To set the boundary matrix determined by buried surface
        VV(1,1)=1.0d0;                  VV(1,2)=s    
        VV(2,1)=r;                      VV(2,2)=1.0d0              
        VV(3,1)=2.0d0*miu(mods%ceng)*r; VV(3,2)=miu(mods%ceng)*t
        VV(4,1)=miu(mods%ceng)*t;       VV(4,2)=2.0d0*miu(mods%ceng)*s

        temp=matmul(XX,VV)

        Fx=abs(real(temp(1,1)*temp(2,2)-temp(1,2)*temp(2,1)))
    End subroutine Haskell_s

    ! Get the Rayleigh wave velocity(rayv) under the high-frequency approximation of the first layer of medium
    Subroutine Crfinder(vs1,vp1)
        Real(kind=8)::vs1,vp1,tol_cr0,c,R,DR
        
        tol_cr0 = 1.0e-7
        c = 0.84*vs1
        R=1;DR=1
        Do while (abs(R/DR)>=tol_cr0)
            Call Rayleigh(R,DR,c,vs1,vp1)
            c = c - R/DR
        End Do
        rayv = c
    End Subroutine crfinder
    Subroutine Rayleigh(R,DR,c,v1,v2)
        Real(kind=8)::R,DR,c,v1,v2,ps,pp,p,p2,ps2,pp2,sps,spp
        ps  = 1.0/v1
        pp  = 1.0/v2
        p   = 1.0/c
        p2  = p*p
        ps2 = ps*ps
        pp2 = pp*pp
        sps = sqrt( p2 - ps2 )
        spp = sqrt( p2 - pp2 )
        R   = ( ps2 - 2.0*p2 )**2 - 4.0*p2*sps*spp
        DR  = p2*( 8.0*p*( ps2 - 2.0*p2 ) + 8.0*p*sps*spp+ 4.0*(p2*p)*( spp/sps + sps/spp ))
    End Subroutine Rayleigh

    ! Get the rough search interval
    Subroutine rough_distance(dc,mods,f)
        Type(Model):: mods
        Real(kind=8)::f,dc,fnv
        fnv=ceiling(ncf(mods,f,mods%vs(mods%ceng)))+1
        dc=(maxvs-minvs)/real(fnv)
    End Subroutine rough_distance

    ! The fine search interval is determined according to the rough search interval
    Subroutine fine_distance(dcout,ndc,mods,f,vin,dcin)
        Type(model)::mods
        Real(kind=8),intent(out)::dcout(20)
        Integer,intent(out)::ndc
        Real(kind=8),intent(in) ::vin,dcin,f
        
        Real(kind=8)::vout,e1,e2
    
        vout=vin+dcin
        If (vout>maxvs) Then
            vout=maxvs
        End If
        e1=ncf(mods,f,vin)
        e2=ncf(mods,f,vout)

        Do while (abs(e1-e2)>ednn)
            vout=vin+(vout-vin)*0.618
            e2=ncf(mods,f,vout)
        End Do
        
        If (rayv>vin .and. rayv<vout) Then
            dcout(1) =(rayv-vin) *(1.0/2.0)
            dcout(2) =(rayv-vin) *(1.0/2.0-1.0/20.0)
            dcout(3) =(rayv-vin) *(1.0/20.0-1.0/200.0)
            dcout(4) =(rayv-vin) *(1.0/200.0-1.0/2000.0)
            dcout(5) =(rayv-vin) *(1.0/2000.0)
            dcout(6) =(vout-rayv)*(1.0/2000.0)
            dcout(7) =(vout-rayv)*(1.0/200.0-1.0/2000.0)
            dcout(8) =(vout-rayv)*(1.0/20.0-1.0/200.0)
            dcout(9) =(vout-rayv)*(1.0/2.0-1.0/20.0)
            dcout(10)=(vout-rayv)*(1.0/2.0)
            ndc=10
        Else
            ndc=10
            dcout(1:ndc)=(vout-vin)/Real(ndc)
        End If
    End Subroutine fine_distance
    
    ! Prediction number of roots(ncf) at frequency(f) and phase velocity(v)
    Function ncf(mods,f,v)
        Type(model)::mods
        Real(kind=8)::ncf,v,f
        Complex(kind=8)::rp(mods%ceng),rs(mods%ceng)
        rp=sqrt(cmplx((v/mods%vp(1:mods%ceng))**2-1.0))
        rs=sqrt(cmplx((v/mods%vs(1:mods%ceng))**2-1.0))
        ncf=Real(f*sum(2.0*(rp(1:mods%ceng-1)+rs(1:mods%ceng-1))*mods%thk(1:mods%ceng-1))/v)
    End Function ncf
    
    
    ! Determinant of x
    Function r_det(x,n)
        Real(kind=8)::x(n,n),r_det
        Integer::n
        If (n==2) Then
            r_det=x(1,1)*x(2,2)-x(1,2)*x(2,1)
        ElseIf (n==3) Then
            r_det=x(1,1)*x(2,2)*x(3,3)+ x(1,2)*x(2,3)*x(3,1)+ x(2,1)*x(3,2)*x(1,3)- &
                  x(1,3)*x(2,2)*x(3,1)- x(1,1)*x(3,2)*x(2,3)- x(2,1)*x(1,2)*x(3,3)
        ElseIf (n==4) Then
            r_det=x(1,1)*x(2,2)*x(3,3)*x(4,4)- x(1,2)*x(2,3)*x(3,4)*x(4,1)+ x(1,3)*x(2,4)*x(3,1)*x(4,2)- &
                  x(1,4)*x(2,1)*x(3,2)*x(4,3)+ x(4,1)*x(3,2)*x(2,3)*x(1,4)- x(4,2)*x(3,3)*x(2,4)*x(1,1)+ &
                  x(4,3)*x(3,4)*x(2,1)*x(1,2)- x(4,4)*x(3,1)*x(2,2)*x(1,3)+ x(1,1)*x(2,3)*x(3,4)*x(4,2)- &
                  x(1,3)*x(2,4)*x(3,2)*x(4,1)+ x(1,4)*x(2,2)*x(3,1)*x(4,3)- x(1,2)*x(2,1)*x(3,3)*x(4,4)+ &
                  x(4,1)*x(3,3)*x(2,4)*x(1,2)- x(4,3)*x(3,4)*x(2,2)*x(1,1)+ x(4,4)*x(3,2)*x(2,1)*x(1,3)- &
                  x(4,2)*x(3,1)*x(2,3)*x(1,4)+ x(1,1)*x(2,4)*x(3,2)*x(4,3)- x(1,4)*x(2,2)*x(3,3)*x(4,1)+ &
                  x(1,2)*x(2,3)*x(3,1)*x(4,4)- x(1,3)*x(2,1)*x(3,4)*x(4,2)+ x(4,1)*x(3,4)*x(2,2)*x(1,3)- &
                  x(4,4)*x(3,2)*x(2,3)*x(1,1)+ x(4,2)*x(3,3)*x(2,1)*x(1,4)- x(4,3)*x(3,1)*x(2,4)*x(1,2)     
        End If
    End Function r_det
    Function z_det(x,n)
        Complex(kind=8)::x(n,n),z_det
        Integer::n
        If (n==2) Then
            z_det=x(1,1)*x(2,2)-x(1,2)*x(2,1)
        ElseIf (n==3) Then
            z_det=x(1,1)*x(2,2)*x(3,3)+ x(1,2)*x(2,3)*x(3,1)+ x(2,1)*x(3,2)*x(1,3)- &
                  x(1,3)*x(2,2)*x(3,1)- x(1,1)*x(3,2)*x(2,3)- x(2,1)*x(1,2)*x(3,3)
        ElseIf (n==4) Then
            z_det=x(1,1)*x(2,2)*x(3,3)*x(4,4)- x(1,2)*x(2,3)*x(3,4)*x(4,1)+ x(1,3)*x(2,4)*x(3,1)*x(4,2)- &
                  x(1,4)*x(2,1)*x(3,2)*x(4,3)+ x(4,1)*x(3,2)*x(2,3)*x(1,4)- x(4,2)*x(3,3)*x(2,4)*x(1,1)+ &
                  x(4,3)*x(3,4)*x(2,1)*x(1,2)- x(4,4)*x(3,1)*x(2,2)*x(1,3)+ x(1,1)*x(2,3)*x(3,4)*x(4,2)- &
                  x(1,3)*x(2,4)*x(3,2)*x(4,1)+ x(1,4)*x(2,2)*x(3,1)*x(4,3)- x(1,2)*x(2,1)*x(3,3)*x(4,4)+ &
                  x(4,1)*x(3,3)*x(2,4)*x(1,2)- x(4,3)*x(3,4)*x(2,2)*x(1,1)+ x(4,4)*x(3,2)*x(2,1)*x(1,3)- &
                  x(4,2)*x(3,1)*x(2,3)*x(1,4)+ x(1,1)*x(2,4)*x(3,2)*x(4,3)- x(1,4)*x(2,2)*x(3,3)*x(4,1)+ &
                  x(1,2)*x(2,3)*x(3,1)*x(4,4)- x(1,3)*x(2,1)*x(3,4)*x(4,2)+ x(4,1)*x(3,4)*x(2,2)*x(1,3)- &
                  x(4,4)*x(3,2)*x(2,3)*x(1,1)+ x(4,2)*x(3,3)*x(2,1)*x(1,4)- x(4,3)*x(3,1)*x(2,4)*x(1,2)   
        End If
    End Function z_det
    
    ! int_sqrt,neg_sqrt
    Function int_sqrt(x)
        Complex(kind=8) ::int_sqrt
        Integer::x
        If (x>0) Then
            int_sqrt=sqrt(Real(x))
        ElseIf (x==0) Then
            int_sqrt=(0.0d0,0.0d0)
        Else
            int_sqrt=sqrt(cmplx(x))
        End If
    End Function int_sqrt
    Function neg_sqrt(x)
        Complex(kind=8) ::neg_sqrt
        Real::x
        neg_sqrt=sqrt(cmplx(x))
    End Function neg_sqrt
    
    Function inv_2(a0)
        Complex(kind=8)::d, a(2, 2), y(2, 2)
        Complex(kind=8)::inv_2(2,2),a0(2, 2)
        Real(kind=8)::epsi02

        epsi02=1.e-8
        a=a0

        d = a(1,1)*a(2,2)-a(1,2)*a(2,1)
        If ( abs(d).lt.epsi02 ) Then
            d = epsi02
        End If

        inv_2(1,1) = a(2,2)/d
        inv_2(2,2) = a(1,1)/d
        inv_2(1,2) =-a(1,2)/d
        inv_2(2,1) =-a(2,1)/d
        return
    End Function inv_2


    Function inv_4(a0)
        Integer::i, j
        Complex(kind=8)::inv_4(4,4),a0(4,4)
	Complex(kind=8)::d, a(4, 4), y(4, 4)
	Complex(kind=8)::a12a23, a12a24, a13a24, a11a23, a11a24, a11a22
	Complex(kind=8)::a33a44, a32a44, a32a43, a31a44, a31a43, a31a42
	Real(kind=8)::epsi04
	
	epsi04=1.e-8
        a=a0

	a12a23 = a(1,2)*a(2,3)-a(1,3)*a(2,2)
	a12a24 = a(1,2)*a(2,4)-a(1,4)*a(2,2)
	a13a24 = a(1,3)*a(2,4)-a(1,4)*a(2,3)
	a11a23 = a(1,1)*a(2,3)-a(1,3)*a(2,1)
	a11a24 = a(1,1)*a(2,4)-a(1,4)*a(2,1)
	a11a22 = a(1,1)*a(2,2)-a(1,2)*a(2,1)
	a33a44 = a(3,3)*a(4,4)-a(3,4)*a(4,3)
	a32a44 = a(3,2)*a(4,4)-a(3,4)*a(4,2)
	a32a43 = a(3,2)*a(4,3)-a(3,3)*a(4,2)
	a31a44 = a(3,1)*a(4,4)-a(3,4)*a(4,1)
	a31a43 = a(3,1)*a(4,3)-a(3,3)*a(4,1)
	a31a42 = a(3,1)*a(4,2)-a(3,2)*a(4,1)

	y(1,1) = a(2,2)*a33a44-a(2,3)*a32a44+a(2,4)*a32a43
	y(2,1) =-a(2,1)*a33a44+a(2,3)*a31a44-a(2,4)*a31a43
	y(3,1) = a(2,1)*a32a44-a(2,2)*a31a44+a(2,4)*a31a42
	y(4,1) =-a(2,1)*a32a43+a(2,2)*a31a43-a(2,3)*a31a42

	y(1,2) =-a(1,2)*a33a44+a(1,3)*a32a44-a(1,4)*a32a43
	y(2,2) = a(1,1)*a33a44-a(1,3)*a31a44+a(1,4)*a31a43
	y(3,2) =-a(1,1)*a32a44+a(1,2)*a31a44-a(1,4)*a31a42
	y(4,2) = a(1,1)*a32a43-a(1,2)*a31a43+a(1,3)*a31a42

	y(1,3) = a(4,4)*a12a23-a(4,3)*a12a24+a(4,2)*a13a24
	y(2,3) =-a(4,4)*a11a23+a(4,3)*a11a24-a(4,1)*a13a24
	y(3,3) = a(4,4)*a11a22-a(4,2)*a11a24+a(4,1)*a12a24
	y(4,3) =-a(4,3)*a11a22+a(4,2)*a11a23-a(4,1)*a12a23

	y(1,4) =-a(3,4)*a12a23+a(3,3)*a12a24-a(3,2)*a13a24
	y(2,4) = a(3,4)*a11a23-a(3,3)*a11a24+a(3,1)*a13a24
	y(3,4) =-a(3,4)*a11a22+a(3,2)*a11a24-a(3,1)*a12a24
	y(4,4) = a(3,3)*a11a22-a(3,2)*a11a23+a(3,1)*a12a23	

	d = a(1,1)*y(1,1)+a(1,2)*y(2,1)+a(1,3)*y(3,1)+a(1,4)*y(4,1)
	If ( abs(d).lt.epsi04 ) Then
	    d = epsi04
	End If

	Do j = 1, 4
	   Do i = 1, 4
	      inv_4(i, j) = y(i,j)/d
	   End Do
	End Do

    End Function inv_4
    
End Module Dispersion
