      PROGRAM CHARACTERISTIC
	  use netcdf

C     This program calculates the rainfall recycling ratio
C     using the characteristic method.
C     Change the domain dimensions in the main program as well as
C     in the rk4 subroutine! 

C     IY : domain grid dimension in I direction (north-south)
C     JX : domain grid dimension in J direction (east-west)
C     KT : Total time steps 
C     IY0: upper-left corner of the domain
C     JX0: upper-left corner of the domain
C     KT0: Initial time (=0 by default) 
C     DY : grid size (m) of I direction (north-south)
C     DX : grid size (m) of J direction (east-west)
C     DA : grid area = DX*DY
C     DT : time interval (second) of input data 
C     U1 : Integrated x-water vapor flux at current time step.
C     U0 : Integrated x-water vapor flux at previous time step.
C     V1 : Integrated y-water vapor flux at current time step.
C     V0 : Integrated y-water vapor flux at previous time step.
C     W1 : Precipitable water content at current time step.
C     W0 : Precipitable water content at previous time step.
C     E1 : Surface evaporation at current time step.
C     E0 : Surface evaporation at previous time step.
C     RR : (Transformed) rainfall recycling ratio.
C     INDFL : inflow/outflow/interior identifier
C            = +1 inflow boundary (RR=0)
C            = -1 outflow boundary 
C            =  0 interior points

!     INTEGER, PARAMETER :: IY=23,JX=47,KT=1460,IY0=11,JX0=79,KT0=0   !(NorthAmerica)
!     INTEGER, PARAMETER :: IY=24,JX=25,KT=1460,IY0=32,JX0=110,KT0=0 !(Amazon)
!     INTEGER, PARAMETER :: IY=17,JX=33,KT=124,IY0=13,JX0=110,KT0=0
      INTEGER, PARAMETER :: IY=25,JX=25,KT=1460,IY0=33,JX0=109,KT0=0 !(SouthAmerica)
      
      integer :: idu, idv, idq, idlh, idpw, idps, idland
      
      INTEGER LAND(JX,IY), LAND2(IY,JX)!Yao Tang Added
      INTEGER INDFL(IY,JX)
      REAL U1(IY,JX),U0(IY,JX),V1(IY,JX),V0(IY,JX)
      REAL W1(IY,JX),W0(IY,JX),E1(IY,JX),E0(IY,JX)
      REAL U11(IY,JX),U00(IY,JX),V11(IY,JX),V00(IY,JX)
      REAL W11(IY,JX),W00(IY,JX),E11(IY,JX),E00(IY,JX)
      REAL RR(IY,JX), RR0(IY,JX), rk

      data DT, DX, DY, DA/2.16E+4, 2.5E+5, 2.5E+5, 6.25E+10/
      
            
      COMMON /C1/ U1, V1, W1, E1
      COMMON /C0/ U0, V0, W0, E0
      COMMON /C11/ U11, V11, W11, E11
      COMMON /C00/ U00, V00, W00, E00
      COMMON /D1/ DT, DX, DY, DA
      COMMON /D2/ RR0 
      COMMON /IDS/idu, idv, idq, idlh, idpw, idps, idland
      
      
      PARAMETER(NC=1)
      DT = DT/float(NC)

      status = nf90_open("uwnd.2003.nc", nf90_nowrite, idu)
      status = nf90_open("vwnd.2003.nc", nf90_nowrite, idv)
      status = nf90_open("shum.2003.nc", nf90_nowrite, idq) 
      status = nf90_open("lhtfl.sfc.gauss.2003.nc", nf90_nowrite, idlh)
      status = nf90_open("pr_wtr.eatm.2003.nc", nf90_nowrite, idpw)
      status = nf90_open("pres.sfc.2003.nc", nf90_nowrite, idps) 
      status = nf90_open("land.nc", nf90_nowrite, idland) 

      OPEN(110,file='RR_ncep03.out',status='unknown',form='formatted')

      OPEN(31,file='U1_ncep03.out',status='unknown',form='formatted')
      OPEN(32,file='V1_ncep03.out',status='unknown',form='formatted')
      OPEN(33,file='W1_ncep03.out',status='unknown',form='formatted')
      OPEN(34,file='E1_ncep03.out',status='unknown',form='formatted')
      
      !======================Yao Tang Add============================!
      !====================consider Land ============================!
      status = nf90_get_var(idland,4,LAND,(/JX0,IY0,1/),(/JX,IY,1/))
      do i=1,IY
          do j=1,JX
          LAND2(i,j) = LAND(j,IY+1-i)
          enddo
      enddo
      !=====================Yao Tang Add=============================!
      
      call init_cond(RR0,U0,V0,W0,E0,IY,JX,IY0,JX0,KT0)
      
      do k=2,KT
         
         call input(k,U1,V1,E1,W1,IY,JX,IY0,JX0,KT0)

            do i=1,IY
               write(31,100) (U1(i,j),j=1,JX)
               write(32,100) (V1(i,j),j=1,JX)
               write(33,100) (W1(i,j),j=1,JX)
               write(34,100) (E1(i,j),j=1,JX)
            enddo
         do i=1,IY
            do j=1,JX
		       U00(i,j)=U0(i,j)
               V00(i,j)=V0(i,j)
               W00(i,j)=W0(i,j)
               E00(i,j)=E0(i,j)
            enddo
         enddo

          do n=1,NC
              do i=1,IY
                  do j=1,JX
                      U11(i,j)=(1-n/NC)*U0(i,j)+n/NC*U1(i,j)
                      V11(i,j)=(1-n/NC)*V0(i,j)+n/NC*V1(i,j)
                      W11(i,j)=(1-n/NC)*W0(i,j)+n/NC*W1(i,j)
                      E11(i,j)=(1-n/NC)*E0(i,j)+n/NC*E1(i,j)
                  enddo
              enddo
                       
              call flowbdy(INDFL,U11,V11,IY,JX)             

              do i=1,IY
                  do j=1,JX
                   if(INDFL(i,j).gt.0) then
                          RR(i,j) = 0.0
                   else 
                          call recycl_cal(rk,i,j)
                          RR(i,j)=rk
                   endif           
                   enddo
              enddo
               do i=1,IY
                 do j=1,JX
                    RR0(i,j)=RR(i,j)
     		        U00(i,j)=U11(i,j)
                    V00(i,j)=V11(i,j)
                    W00(i,j)=W11(i,j)
                    E00(i,j)=E11(i,j)
                 enddo
              enddo
           enddo
 
         do i=1,IY
            do j=1,JX
               RR0(i,j)=RR(i,j)
		       U0(i,j)=U1(i,j)
               V0(i,j)=V1(i,j)
               W0(i,j)=W1(i,j)
               E0(i,j)=E1(i,j)
            enddo
         enddo
      !======Yao Tang Add=========!
          do i=1,IY
            do j=1,JX
               RR(i,j)=RR(i,j)*LAND2(i,j)
            enddo
         enddo
      !======Yao Tang Add=========!   
         do i=1,IY
            write(110,100) (RR(i,j),j=1,JX)
         enddo

      enddo

 100  format(1x,25(e10.4,1x))

      stop
      end



C---- calculating recycling ratio with two steps:
C
C     Step 1: find the characteristic, C, passing through (i,j) at time k,
C             over the period of [k-1,k],
C     Step 2: computing the C-integral of E/W.
C
C     NC : no of time steps of characteristic calculation over [k-1,k]
C     XC : x-coordinate of the characteristic
C     YC : y-coordinate of the characteristic
C     SC : source function, E/W, along the characteristic, C
C     rl : local recycling ratio at (i,j) at current time, k. 
C     r0 : local recycling ratio at (i,j) at previous time, k-1. 


      subroutine recycl_cal(rl,i0,j0)

      PARAMETER(NC=2)

      INTEGER i0, j0
      REAL XC(NC), YC(NC), SC(NC), rl, r0

      call charact(i0,j0,r0,SC,NC)

      call cintegral(rl,r0,SC,NC)

      return
      end



C---- Step 1 subroutine: get the characteristic, C.

      subroutine charact(i,j,r0,SC,NC)

      INTEGER NC, i, j, k, m, m0      
      REAL    DT, DX, DY, DA
      REAL    xc(NC), yc(NC), sc(NC), xm, ym, sm, h, r0

      COMMON /D1/ DT, DX, DY, DA

      h=DT/float(NC-1)

      do m=1,NC
         
         if (m.eq.1) then
           call rk4(xm,ym,sm,r0,0.0,0.0,m,i,j,h,NC)
         else
           call rk4(xm,ym,sm,r0,xc(m-1),yc(m-1),m,i,j,h,NC)
         endif

         if(abs(xm).gt.1.0E+08.or.abs(ym).gt.1.0E+08) then
            do m0=m,NC
               sc(m0)= 0.0
            enddo
            return
         endif

         xc(m)=xm
         yc(m)=ym
         sc(m)=sm

      enddo

      return
      end



C---- sub 1 of Step 1 :
C     Runge-Kutta method for calculating the characteristic

      subroutine rk4(xm,ym,sm,r0,xm0,ym0,m,i,j,h,NC)

      PARAMETER(IY=25,JX=25)

      REAL U11(IY,JX),U00(IY,JX),V11(IY,JX),V00(IY,JX)
      REAL W11(IY,JX),W00(IY,JX),E11(IY,JX),E00(IY,JX)

      REAL RR0(IY,JX)
      REAL DT, DX, DY, DA

      real kx(4), ky(4), al(4), cf(4), w(4)
      real xm0,ym0,h,x1,y1,r0
      integer il(4), jl(4)
      integer ii, i, j, m, id, jd, NC
      integer imax, imin, jmax, jmin

      COMMON /C11/ U11, V11, W11, E11
      COMMON /C00/ U00, V00, W00, E00
      COMMON /D1/ DT, DX, DY, DA
      COMMON /D2/ RR0

      data cf /0., 0.5, 0.5, 1.0/
      data w  /6., 3.0, 3.0, 6.0/

      if (m.eq.1) then
         xm=xm0
         ym=ym0
         sm=E11(i,j)
	     r0=0.0
         return
      elseif (m.eq.2) then
         call localno(il,jl,id,jd,i,j,U11(i,j),V11(i,j))
      else
         call distantno(il,jl,id,jd,i,j,-xm0,-ym0)
      endif   

      imax = max(il(1),il(2),il(3),il(4))
      imin = min(il(1),il(2),il(3),il(4))
      jmax = max(jl(1),jl(2),jl(3),jl(4))
      jmin = min(jl(1),jl(2),jl(3),jl(4))

      if (imin.lt.1) then
         sm = 0.0
         xm = 0.0
         ym = -5.0E+08
         r0 = 0.0
         return
       elseif (imax.gt.IY) then
         sm = 0.0
         xm = 0.0
         ym = 5.0E+08
         r0 = 0.0
         return
       elseif (jmin.lt.1) then
         sm = 0.0
         xm = -5.0E+08
         ym = 0.0
         r0 = 0.0
         return
       elseif (jmax.gt.JX) then
         sm = 0.0
         xm = 5.0E+08
         ym = 0.0
         r0 = 0.0
         return
       endif

      sjd= float(-jd)
      sid= float(-id)

      x1 = xm0
      y1 = ym0

      call weight_interp(al,x1,y1,sid,sjd)

      call space_inter(u1m,il,jl,al,IY,JX,U11)
      call space_inter(v1m,il,jl,al,IY,JX,V11)

      kx(1) = -u1m*h
      ky(1) = -v1m*h

      do ii=2,4
         x1=xm0+cf(ii)*kx(ii-1)
         y1=ym0+cf(ii)*ky(ii-1)

         call weight_interp(al,x1,y1,sid,sjd)

         call space_inter(u1m,il,jl,al,IY,JX,U11)
         call space_inter(v1m,il,jl,al,IY,JX,V11)
         call space_inter(u0m,il,jl,al,IY,JX,U00)
         call space_inter(v0m,il,jl,al,IY,JX,V00)
      
         ta = (float(m-2)+cf(ii))*h/DT
      
         kx(ii) = -(u1m*(1.0-ta)+u0m*ta)*h
         ky(ii) = -(v1m*(1.0-ta)+v0m*ta)*h

      enddo

      xm=0.0
      ym=0.0
      do ii=1,4
       xm=xm+kx(ii)/w(ii)
       ym=ym+ky(ii)/w(ii)
      enddo
      xm=xm+xm0
      ym=ym+ym0

      x1=xm
      y1=ym

      call weight_interp(al,x1,y1,sid,sjd)

      ta = float(m-1)*h/DT

      call space_inter(e1m,il,jl,al,IY,JX,E11)
      call space_inter(e0m,il,jl,al,IY,JX,E00)

      sm = e1m*(1.0-ta)+e0m*ta

      if (m.eq.NC) then
         call space_inter(r0,il,jl,al,IY,JX,RR0)
      else
         r0 = 0.0
      endif

      return
      end


C---- sub 2 of step 1: spatial interpolation

      subroutine space_inter(z,il,jl,al,IY,JX,V)

      integer IY, JX, il(4), jl(4)
      real    al(4), V(IY,JX), z

      z = V(il(1),jl(1))*al(3)+V(il(2),jl(2))*al(4)+
     &    V(il(3),jl(3))*al(1)+V(il(4),jl(4))*al(2)

      return
      end



C---- sub 3 of step 1: spatial interpolation weights

      subroutine weight_interp(al,x1,y1,sid,sjd)

      real al(4)
      real x1, y1, sid, sjd, DX, DY, DA, xl, yl

      COMMON /D1/ DT, DX, DY, DA

      xl=x1*sjd-floor(x1*sjd/DX)*DX
      yl=y1*sid-floor(y1*sid/DY)*DY

      al(1)=xl*yl/DA
      al(2)=(DX-xl)*yl/DA
      al(3)=(DX-xl)*(DY-yl)/DA
      al(4)=xl*(DY-yl)/DA

      return
      end


C---- sub 4 of step 1: get local node numbers at (i, j)

      subroutine localno(il,jl,id,jd,i,j,u,v)

      real u, v
      integer il(4), jl(4), id, jd, i, j

      if (u.eq.0.0) then
         jd = 0
      else
         jd = int(u/abs(u))
      endif

      if (v.eq.0.0) then
         id = 0
      else
         id = int(v/abs(v))
      endif

      il(1) = i
      il(2) = i
      il(3) = i-id
      il(4) = i-id

      jl(1) = j     
      jl(2) = j-jd
      jl(3) = j-jd
      jl(4) = j

      return
      end



C---- sub 5 of step 1: get distant node numbers at (i, j)
C                      for the case the trajectory goes beyond
C                      the neighboring grids of (i,j)

      subroutine distantno(il,jl,id,jd,i,j,u,v)

      real u, v
      integer il(4), jl(4), id, jd, i, j, i0, j0, is, js
      real DT, DX, DY, DA

      COMMON /D1/DT, DX, DY, DA

      if (u.eq.0.0) then
         jd = 0
      else
         jd = int(u/abs(u))
      endif

      if (v.eq.0.0) then
         id = 0
      else
         id = int(v/abs(v))
      endif

      is = -id*floor(abs(v)/DY)
      js = -jd*floor(abs(u)/DX)

      i0 = i+is
      j0 = j+js

      il(1) = i0
      il(2) = i0
      il(3) = i0-id
      il(4) = i0-id

      jl(1) = j0     
      jl(2) = j0-jd
      jl(3) = j0-jd
      jl(4) = j0

      return
      end


C---- Step 2 subroutine: characteristic integral


      subroutine cintegral(r,r0,SC,NC)

      INTEGER NC
      REAL DT, DX, DY, DA
      REAL SC(NC),SCC(NC), h, r, r0

      COMMON /D1/ DT, DX, DY, DA

      h=DT/float(NC-1)
      
      do i=1,NC
          SCC(i)=SC(NC+1-i)
      enddo
      do i=1,NC
          SC(i)=SCC(i)
      enddo
      
      r = r0
      do i=1,NC-1
         r = (r+0.5*(SC(i+1)+(1-r)*SC(i))*h)/(1+0.5*SC(i+1)*h)
      enddo
      
      return
      end



C---- assign initial value of RR ---------------
C     dim_1 = lon  (x-axis, JX)
C     dim_2 = lat  (y-axis, IY)
C     dim_3 = lev  (p-axis or z-axis, NZ)
C     dim_4 = time (t-axis, KT)
C     dim_5 = var  


      subroutine init_cond(R,U0,V0,W0,E0,IY,JX,IY0,JX0,KT0)
      use netcdf

      REAL R(IY,JX),U0(IY,JX),V0(IY,JX),W0(IY,JX),E0(IY,JX)

      integer :: IY, JX, IY0, JX0, KT0
      integer :: idu, idv, idq, idlh, idpw, idps, idland
      integer, parameter :: NZ=8, lon=192, lat=94
      integer, allocatable, dimension(:,:) :: LAND, LAND2 
      real, allocatable, dimension(:,:,:)  :: U, V, Q
      real, allocatable, dimension(:,:)    :: PW, PS, LH
      real, allocatable, dimension(:,:)    :: U2, V2, E2, W2, Q0
      real, allocatable, dimension(:)      :: xlon, ylat, ylat1
      real, dimension(73,144) :: EG

      COMMON /IDS/idu, idv, idq, idlh, idpw, idps, idland

      allocate(U(JX,IY,NZ), V(JX,IY,NZ), Q(JX,IY,NZ))
      allocate(PW(JX,IY), PS(JX,IY), LH(lon,lat), LAND(JX,IY))
      allocate(xlon(JX), ylat(IY), ylat1(IY) )
      allocate(Q0(IY,JX))

      status = nf90_get_var(idu,5,U,(/JX0,IY0,1,KT0+1/),(/JX,IY,NZ,1/))
      status = nf90_get_var(idv,5,V,(/JX0,IY0,1,KT0+1/),(/JX,IY,NZ,1/))
      status = nf90_get_var(idq,5,Q,(/JX0,IY0,1,KT0+1/),(/JX,IY,NZ,1/))

      status = nf90_get_var(idq,2,ylat,(/IY0/),(/IY/))
      status = nf90_get_var(idq,3,xlon,(/JX0/),(/JX/))

      status = nf90_get_var(idlh,4,LH,(/1,1,KT0+1/),(/lon,lat,1/))
      status = nf90_get_var(idpw,4,PW,(/JX0,IY0,KT0+1/),(/JX,IY,1/))
      status = nf90_get_var(idps,4,PS,(/JX0,IY0,KT0+1/),(/JX,IY,1/))
      status = nf90_get_var(idland,4,LAND,(/JX0,IY0,1/),(/JX,IY,1/))

      call rescale(U,V,Q,PW,PS,LH,JX,IY,NZ,lon,lat)

      call int_vert(U0,V0,Q0,U,V,Q,PS,IY,JX,NZ)

      call Emap(EG,LH,73,144,1,1,lon,lat)

      do i=1,IY
         do j=1,JX
            E0(i,j)=EG(IY0+i-1,JX0+j-1)
         enddo
      enddo

      do i=1,IY
         do j=1,JX
            R(i,j)=1;
            W0(i,j) = Q0(i,j)
         enddo
      enddo

      allocate(U2(IY,JX),V2(IY,JX),E2(IY,JX),W2(IY,JX),LAND2(IY,JX))

      do i=1,IY
         ylat1(i) = ylat(IY+1-i)
         do j=1,JX
            U2(i,j) = U0(IY+1-i,j)
            V2(i,j) = V0(IY+1-i,j)
            E2(i,j) = E0(IY+1-i,j)
            W2(i,j) = W0(IY+1-i,j)
            LAND2(i,j) = LAND(j,IY+1-i)
         enddo
      enddo

      do i=1,IY
         ylat(i) = ylat1(i)
         do j=1,JX
            U0(i,j) = U2(i,j)/W2(i,j)
            V0(i,j) = V2(i,j)/W2(i,j)
            E0(i,j) = E2(i,j)/W2(i,j)
            W0(i,j) = W2(i,j)
            !========Yao Tang added========!
            E0(i,j) = E0(i,j) * LAND2(i,j)
            R(i,j) = R(i,j) * LAND2(i,j)
            !========Yao Tang added========!
         enddo
      enddo
    
      OPEN(11,file='ylat_ncep03.out',status='unknown',form='formatted')
      OPEN(12,file='xlon_ncep03.out',status='unknown',form='formatted')
      OPEN(13,file='land_ncep03.out',status='unknown',form='formatted')
         write(11,100)(ylat(i),i=1,IY)
         write(12,100)(xlon(j),j=1,JX)
         do i=1,IY
            write(13,200)(LAND2(i,j),j=1,JX)
         enddo
      close(11)
      close(12)
      close(13)

      deallocate(ylat,ylat1,xlon)
      deallocate(U2,V2,E2,W2,LAND2,Q0)
      deallocate(U, V, Q, PW, PS, LH, LAND)

 100  format(50(1x,E14.7))
 200  format(50(1x,I2))

      return
      end



C---- transform the raw netCDF data back to their true values 
C     using scale and offset factors

      subroutine rescale(U,V,Q,PW,PS,LH,JX,IY,NZ,lon,lat)

      integer :: JX, IY, NZ, lon, lat, i, j, k
      real, dimension(JX,IY,NZ) :: U, V, Q
      real, dimension(JX,IY)    :: PW, PS
      real, dimension(lon,lat)  :: LH
      real :: add_q,  scl_q, add_uv, scl_uv, add_lh, scl_lh
      real :: add_pw, scl_pw, add_ps, scl_ps,cc

      data add_q,  scl_q /3.2666E-02, 1.0E-06/ 
      data add_uv, scl_uv/202.66, 1.0E-02/
      data add_lh, scl_lh/1276.5, 0.1/
      data add_pw, scl_pw/277.65, 1.0E-02/
      data add_ps, scl_ps/367650., 10./
      data cc/4.0E-07/

      do i=1,IY
         do j=1,JX
            do k=1,NZ
              U(j,i,k)=U(j,i,k)*scl_uv+add_uv
              V(j,i,k)=V(j,i,k)*scl_uv+add_uv
              Q(j,i,k)=Q(j,i,k)*scl_q+add_q
            end do
              PW(j,i)=PW(j,i)*scl_pw+add_pw
              PS(j,i)=PS(j,i)*scl_ps+add_ps
         enddo
      enddo   
        
      do i=1,lat
         do j=1,lon
            LH(j,i)=LH(j,i)*scl_lh+add_lh
            LH(j,i)=LH(j,i)*cc! Yao Tang Add
         end do
      end do

      return
      end




C---- read input data of H2O vapor fluxes, prec. water and evap
C     from input.dat at time step k

      subroutine input(k,U1,V1,E1,W1,IY,JX,IY0,JX0,KT0)
      use netcdf

      integer, parameter :: NZ=8, lon=192, lat=94
      integer :: idu, idv, idq, idlh, idpw, idps
      integer :: k, IY, JX, IY0, JX0, KT0
      REAL U1(IY,JX),U0(IY,JX),V1(IY,JX),V0(IY,JX)
      REAL W1(IY,JX),W0(IY,JX),E1(IY,JX),E0(IY,JX)
      real, dimension(73,144) :: EG

      real, allocatable, dimension(:,:,:) :: U, V, Q
      real, allocatable, dimension(:,:)   :: PW, PS, LH
      real, allocatable, dimension(:,:)   :: U2, V2, E2, Q1
      integer, allocatable, dimension(:,:) :: LAND, LAND2 !add

      COMMON /IDS/idu, idv, idq, idlh, idpw, idps, idland

      allocate(U(JX,IY,NZ), V(JX,IY,NZ), Q(JX,IY,NZ))
      allocate(PW(JX,IY), PS(JX,IY), LH(lon,lat), LAND(JX,IY))
      allocate(Q1(IY,JX))

      status = nf90_get_var(idu,5,U,(/JX0,IY0,1,KT0+k/),(/JX,IY,NZ,1/))
      status = nf90_get_var(idv,5,V,(/JX0,IY0,1,KT0+k/),(/JX,IY,NZ,1/))
      status = nf90_get_var(idq,5,Q,(/JX0,IY0,1,KT0+k/),(/JX,IY,NZ,1/))

      status = nf90_get_var(idlh,4,LH,(/1,1,KT0+k/),  (/lon,lat,1/))
      status = nf90_get_var(idpw,4,PW,(/JX0,IY0,KT0+k/),(/JX,IY,1/))
      status = nf90_get_var(idps,4,PS,(/JX0,IY0,KT0+k/),(/JX,IY,1/))
      
      status = nf90_get_var(idland,4,LAND,(/JX0,IY0,1/),(/JX,IY,1/))

      call rescale(U,V,Q,PW,PS,LH,JX,IY,NZ,lon,lat)
      call int_vert(U1,V1,Q1,U,V,Q,PS,IY,JX,NZ)
      call Emap(EG,LH,73,144,1,1,lon,lat)

      do i=1,IY
         do j=1,JX
            E1(i,j)=EG(IY0+i-1,JX0+j-1)
         enddo
      enddo

      do i=1,IY
         do j=1,JX

            if (abs(PW(j,i)).lt.1.0E-05) then
               write(*,*)'i=',i,' j=',j,' : PW(j,i)=',PW(j,i)
            endif

            W1(i,j) = Q1(IY+1-i,j)
         enddo
      enddo

      allocate(U2(IY,JX),V2(IY,JX),E2(IY,JX), LAND2(IY,JX))

      do i=1,IY
         do j=1,JX
            U2(i,j) = U1(IY+1-i,j)
            V2(i,j) = V1(IY+1-i,j)
            E2(i,j) = E1(IY+1-i,j)
            LAND2(i,j) = LAND(j,IY+1-i)
         enddo
      enddo

      do i=1,IY
         do j=1,JX

            if (abs(W1(i,j)).lt.1.0E-05) then
               write(*,*)'i=',i,' j=',j,' : W1(i,j)=',W1(i,j)
               W1(i,j) = 20.0
            endif

            U1(i,j) = U2(i,j)/W1(i,j)
            V1(i,j) = V2(i,j)/W1(i,j)
            E1(i,j) = E2(i,j)/W1(i,j)
            !======Yao Tang=============!
            E1(i,j) = E1(i,j)*LAND2(i,j)
            !======Yao Tang=============!
         enddo
      enddo

      deallocate(U2,V2,E2,Q1, LAND2)
      deallocate(U, V, Q, PW, PS, LH, LAND)
     
      return
      end



C---- water vapor fluxes integrated vertically

      subroutine int_vert(UQ,VQ,QQ,U,V,Q,PS,IY,JX,NZ)

      integer :: IY, JX, NZ, i, j, k

      real P(8)
      real, dimension(IY,JX)          :: UQ, VQ, QQ
      real, dimension(JX,IY)          :: PS
      real, dimension(JX,IY,NZ)       :: U, V, Q
      real :: P0, g, dP, dP0, c1, c2, u0, v0, q0

      data P/1.E+5,.925E+5,.85E+5,.7E+5,.6E+5,.5E+5,.4E+5,.3E+5/
      data P0, g/1.0E+05, 9.8/
      
      do i=1,IY
         do j=1,JX
            UQ(i,j) = 0.0
            VQ(i,j) = 0.0
            QQ(i,j) = 0.0
         end do
      end do
      
      do i=1,IY
         do j=1,JX

            if (PS(j,i).ge.P(1)) then
               k0 = 1
               P(1) = PS(j,i)
               goto 101
            endif

            do k=2,NZ
               if (PS(j,i).le.P(k-1).and.PS(j,i).gt.P(k)) then
                  k0 = k-1
                  goto 101
               end if
            end do
            return

 101        do k=k0+1,NZ-1
               dP = (P(k)-P(k+1))/g
               UQ(i,j)=UQ(i,j)+0.5*dP*(Q(j,i,k)*U(j,i,k) +
     &                                 Q(j,i,k+1)*U(j,i,k+1))
               VQ(i,j)=VQ(i,j)+0.5*dP*(Q(j,i,k)*V(j,i,k) + 
     &                                 Q(j,i,k+1)*V(j,i,k+1))
               QQ(i,j)=QQ(i,j)+0.5*dP*(Q(j,i,k) + Q(j,i,k+1))
            end do

            dP0= P(k0)-P(k0+1)
            c1 = (PS(j,i)-P(k0+1))/dP0
            c2 = (P(k0)-PS(j,i))/dP0
            q0 = c1*Q(j,i,k0)+c2*Q(j,i,k0+1)
            u0 = c1*U(j,i,k0)+c2*U(j,i,k0+1)
            v0 = c1*V(j,i,k0)+c2*V(j,i,k0+1)

            dP0= (PS(j,i)-P(k0+1))/g

            UQ(i,j)=UQ(i,j)+0.5*dP0*(q0*u0 + Q(j,i,k0+1)*U(j,i,k0+1))
            VQ(i,j)=VQ(i,j)+0.5*dP0*(q0*v0 + Q(j,i,k0+1)*V(j,i,k0+1))
            QQ(i,j)=QQ(i,j)+0.5*dP0*(q0 + Q(j,i,k0+1))

         end do
      end do

      return
      end


C---- Mapping E on T62 grid onto 144x73 grid

      subroutine Emap(E,LH,IY,JX,IY0,JX0,lon,lat)

      integer :: IY,JX,IY0,JX0, lon,lat, i, j, i1, j1
      integer,    dimension(4) :: il, jl
      real,       dimension(4) :: a
      real, dimension(IY,JX)   :: E
      real, dimension(lon,lat) :: LH
      real :: dlon, dlat, d0, d, ci, cj, da 

      data dlon, dlat, d0, d/1.87487, 1.904086, 1.46, 2.5/

      da = dlon*dlat

      do i=1,IY

        if (i.eq.1.or.i.eq.IY) then
           do j=1,JX
               E(i,j)=0.0
           end do
        else

           do j=1,JX

              i1 = ceiling((float(IY0+i-2)*d-d0)/dlat)
              if (JX0+j-1.eq.1) then
                 j1=1
              else   
                 j1 = ceiling(float(JX0+j-2)*d/dlon)
              endif

              ci = float(IY0+i-2)*d/dlat-float(i1-1)-d0/dlat
              cj = float(JX0+j-2)*d/dlon-float(j1-1)

              a(1) = ci*cj
              a(2) = ci*(1.-cj)
              a(3) = (1.-ci)*(1.-cj)
              a(4) = (1.-ci)*cj

              il(1) = i1
              jl(1) = j1
              il(2) = i1
              jl(2) = j1+1
              il(3) = i1+1
              jl(3) = j1+1
              il(4) = i1+1
              jl(4) = j1

              E(i,j)=LH(jl(1),il(1))*a(3)+LH(jl(2),il(2))*a(4)
     &             + LH(jl(3),il(3))*a(1)+LH(jl(4),il(4))*a(2)
           end do
        endif
      end do
      return
      end





C---- find inflow/outflow boundary points at current time

      subroutine flowbdy(INDFL,U1,V1,IY,JX)

      INTEGER INDFL(IY,JX)
      REAL U1(IY,JX),V1(IY,JX),u,v

      do 200 i=1,IY
         do 100 j=1,JX
            u=U1(i,j)
            v=V1(i,j)

            if(i.eq.1.and.j.eq.1) then
               if(u.lt.0.0.and.v.lt.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(i.eq.1.and.j.gt.1.and.j.lt.JX) then
               if(v.lt.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(i.eq.1.and.j.eq.JX) then
               if(u.ge.0.0.and.v.le.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(j.eq.1.and.i.gt.1.and.i.lt.IY) then
               if(u.lt.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(j.eq.JX.and.i.gt.1.and.i.lt.IY) then
               if(u.gt.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(j.eq.1.and.i.eq.IY) then
               if(u.le.0.0.and.v.ge.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(i.eq.IY.and.j.gt.1.and.j.lt.JX) then

               if(v.gt.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif
             
             if(i.eq.IY.and.j.eq.JX) then
               if(u.ge.0.0.and.v.ge.0.0) then
                  INDFL(i,j)=-1
                  goto 100
               else
                  INDFL(i,j)=+1
                  goto 100
               endif
             endif             

             INDFL(i,j)=0 

 100      enddo
 200   enddo
           
      return
      end


