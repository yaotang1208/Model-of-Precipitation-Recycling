      PROGRAM CHARACTERISTIC
C	  use netcdf

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

!     INTEGER, PARAMETER :: IY=17,JX=41,KT=1460,IY0=13,JX0=81,KT0=0   (US)
      INTEGER, PARAMETER :: IY=41,JX=61,KT=240 !(Central US)
!     INTEGER, PARAMETER :: IY=17,JX=33,KT=124,IY0=13,JX0=110,KT0=0

      integer :: idu, idv, idq, idlh, idpw, idps, idland

      INTEGER INDFL(IY,JX)
      REAL U1(IY,JX),U0(IY,JX),V1(IY,JX),V0(IY,JX)
      REAL W1(IY,JX),W0(IY,JX),E1(IY,JX),E0(IY,JX)
      REAL P0(IY,JX),P1(IY,JX)
      REAL RR(IY,JX), RR0(IY,JX), rk

      data DT, DX, DY, DA/1.08E+4, 2.21E+4, 2.77E+4, 6.12E+8/
                  
      COMMON /C1/ U1, V1, W1, E1
      COMMON /C0/ U0, V0, W0, E0
      COMMON /D1/ DT, DX, DY, DA
      COMMON /D2/ RR0 
      COMMON /IDS/ idu, idv, idq, idlh, idpw, idps, idland
            
      OPEN(21,file='E.txt',status='old')
      OPEN(22,file='PCP.txt',status='old')
      OPEN(23,file='U199806.txt',status='old')
      OPEN(24,file='V199806.txt',status='old') 
      OPEN(25,file='Q199806.txt',status='old')       

      OPEN(110,file='RR_ncep03.out',status='unknown',form='formatted')
      
      OPEN(31,file='U1_ncep03.out',status='unknown',form='formatted')
      OPEN(32,file='V1_ncep03.out',status='unknown',form='formatted')
      OPEN(33,file='W1_ncep03.out',status='unknown',form='formatted')
      OPEN(34,file='E1_ncep03.out',status='unknown',form='formatted')
      
      OPEN(41,file='U0.dat',status='unknown',form='formatted')
      OPEN(42,file='V0.dat',status='unknown',form='formatted')
      OPEN(43,file='W0.dat',status='unknown',form='formatted')
      OPEN(44,file='E0.dat',status='unknown',form='formatted')
      
      do i=1,IY
            read(21,*) (E0(i,j),j=1,JX)
            read(22,*) (P0(i,j),j=1,JX)
            read(23,*) (U0(i,j),j=1,JX)
            read(24,*) (V0(i,j),j=1,JX)
            read(25,*) (W0(i,j),j=1,JX)
      enddo
      
      write(41,100) ((U0(i,j),j=1,JX),i=1,IY)
      write(42,100) ((V0(i,j),j=1,JX),i=1,IY)
      write(43,100) ((W0(i,j),j=1,JX),i=1,IY)
      write(44,100) ((E0(i,j),j=1,JX),i=1,IY)

      do i=1,IY
        do j=1,JX
		    RR0(i,j)=0
        enddo
      enddo
      !===Yao Tang added=====!
      do i=1,IY
        do j=1,JX
		    U0(i,j)=U0(i,j)/W0(i,j)
		    V0(i,j)=V0(i,j)/W0(i,j)
		    E0(i,j)=E0(i,j)/W0(i,j)
		    !E0(i,j)=E0(i,j)/W0(i,j)/DT!!
        enddo
      enddo
      !===Yao Tang added=====!
      do k=2,KT
         
           do i=1,IY
                 read(21,*) (E1(i,j),j=1,JX)
                 read(22,*) (P1(i,j),j=1,JX)
                 read(23,*) (U1(i,j),j=1,JX)
                 read(24,*) (V1(i,j),j=1,JX)
                 read(25,*) (W1(i,j),j=1,JX)
           enddo
           
            do i=1,IY
               write(31,100) (U1(i,j),j=1,JX)
               write(32,100) (V1(i,j),j=1,JX)
               write(33,100) (W1(i,j),j=1,JX)
               write(34,100) (E1(i,j),j=1,JX)
            enddo
            
       !===Yao Tang added=====!     
           do i=1,IY
             do j=1,JX
     		    U1(i,j)=U1(i,j)/W1(i,j)
     		    V1(i,j)=V1(i,j)/W1(i,j)
     		    E1(i,j)=E1(i,j)/W1(i,j)
     		    !E1(i,j)=E1(i,j)/W0(i,j)/DT!!
             enddo
           enddo
       !===Yao Tang added=====!
       
         call flowbdy(INDFL,U1,V1,IY,JX)             

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
		       U0(i,j)=U1(i,j)
               V0(i,j)=V1(i,j)
               W0(i,j)=W1(i,j)
               E0(i,j)=E1(i,j)
            enddo
         enddo
  
         do i=1,IY
            write(110,100) (RR(i,j),j=1,JX)
         enddo

      enddo

 100  format(1x,61(e10.4,1x))

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

      PARAMETER(IY=41,JX=61)

      REAL U1(IY,JX),U0(IY,JX),V1(IY,JX),V0(IY,JX)
      REAL W1(IY,JX),W0(IY,JX),E1(IY,JX),E0(IY,JX)

      REAL RR0(IY,JX)
      REAL DT, DX, DY, DA

      real kx(4), ky(4), al(4), cf(4), w(4)
      real xm0,ym0,h,x1,y1,r0
      integer il(4), jl(4)
      integer ii, i, j, m, id, jd, NC
      integer imax, imin, jmax, jmin

      COMMON /C1/ U1, V1, W1, E1
      COMMON /C0/ U0, V0, W0, E0
      COMMON /D1/ DT, DX, DY, DA
      COMMON /D2/ RR0

      data cf /0., 0.5, 0.5, 1.0/
      data w  /6., 3.0, 3.0, 6.0/

      if (m.eq.1) then
         xm=xm0
         ym=ym0
         sm=E1(i,j)
	     r0=0.0
         return
      elseif (m.eq.2) then
         call localno(il,jl,id,jd,i,j,U1(i,j),V1(i,j))
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

      call space_inter(u1m,il,jl,al,IY,JX,U1)
      call space_inter(v1m,il,jl,al,IY,JX,V1)

      kx(1) = -u1m*h
      ky(1) = -v1m*h

      do ii=2,4
         x1=xm0+cf(ii)*kx(ii-1)
         y1=ym0+cf(ii)*ky(ii-1)

         call weight_interp(al,x1,y1,sid,sjd)

         call space_inter(u1m,il,jl,al,IY,JX,U1)
         call space_inter(v1m,il,jl,al,IY,JX,V1)
         call space_inter(u0m,il,jl,al,IY,JX,U0)
         call space_inter(v0m,il,jl,al,IY,JX,V0)
      
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

      call space_inter(e1m,il,jl,al,IY,JX,E1)
      call space_inter(e0m,il,jl,al,IY,JX,E0)

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