!  Copyright (C) <2021>  <Bo Yang. Email:yuanxzo@qq.com>
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


Program Main
    Use Dispersion
    
    Type(Model):: mods
    Real(kind=8)::freq(100),Vr(100,10)
    Integer::i
    
    mods.ceng=3
    mods.vs(1:3) =[150.0d0,200.0d0,380.0d0]
    mods.vp(1:3) =[367.0d0,437.0d0,755.0d0]
    mods.dns(1:3)=[1.800d0,2.000d0,2.100d0]
    mods.thk(1:2)=[3.000d0,6.000d0]
    
    ! 1 to 100 Hz
    Do i=1,100
        freq(i)=i*1.000d0
    End Do
    
    ! Obtain the dispersion curves of mods
    ! Vr is a matrix of nf * nv, in this case, nf = 100, nv = 10. So Vr is an array of 100 rows and 10 columns.
    Call FVTA_c(Vr,mods,freq,100,10)
    Write(*,*) "<Test of Dispersion>  Copyright (C) <2021>  <Bo Yang. Email:yuanxzo@qq.com>"
    Write(*,*) "This program comes with ABSOLUTELY NO WARRANTY."
    Write(*,*) "This is free software, and you are welcome to redistribute it"
    Write(*,*) "under certain conditions; see <https://www.gnu.org/licenses/> for details."
    Write(*,*) 
    
    Do i=1,100
        Write(*,"(I3,10(f18.9))") i,Vr(i,:)
    End Do
    
    Write(*,*) 
    Write(*,*) "<Test of Dispersion>  Copyright (C) <2021>  <Bo Yang. Email:yuanxzo@qq.com>"
    Write(*,*) "This program comes with ABSOLUTELY NO WARRANTY."
    Write(*,*) "This is free software, and you are welcome to redistribute it"
    Write(*,*) "under certain conditions; see <https://www.gnu.org/licenses/> for details."
    
    Read(*,*)
End
