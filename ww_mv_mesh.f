c----------------------------------------------------------------------
      subroutine make_wire_wrap(nrings,Pptch,Wptch,Dwire,Rcan,Rf0,sunk0,
     &           PBID,HCID)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT' 

c     Inputs
      integer nrings !number of pin layers/rings
      integer PBID   !boundary ID of the pins
      integer HCID   !boundary ID of the hex-can
      real Pptch !center to center distance of pins
      real Wptch !height of a single wire helix (negative for LH)
      real Dwire !diameter of the wire
      real Rf0   !radius of the wire/pin filet
      real Rcan  !radius for rounding the corners of the hexcan
      real sunk0 !distance to sink the wire into the pin

c     variables for adding the wires
      integer ipin,jpin,lpins,niters
      integer i,i0,i1,j,j0,j1,k,k0,k1,iel,ifc,n,ilayer
      parameter(lpins=271) 

      real pio3,pio6 !pi/3,pi/6
      real Rp,Rw,Rf,sunk,d1,d2,d3,theta1,theta2,theta3 !geometric constants
      real xmn,xmx,ymn,ymx,zmn,zmx,apoth,nturns !additional geometric constants
      real delT,T1,T2,T1_x,T1_y,T2_x,T2_y,T3,T4,Ttot,psi !arc length constants for step 1
      real S3,S4,S5,S6,Stot !arc length constants for step 2
      real xx,yy,zz,rr,thw,xpt,ypt,zpt,theta,thloc,xnew,ynew,xr,yr !local variables

      integer pindx(lx1,ly1,lz1,lelt) !pin ID number
      real xxc(lpins),yyc(lpins) !pin center coordinates
      real spt(lx1,ly1,lz1,lelt),sptl !arc length
      real delx(lx1,ly1,lz1,lelt) !x-displacement 
      real dely(lx1,ly1,lz1,lelt) !y-displacement
      real delz(lx1,ly1,lz1,lelt) !z-displacement (not needed)
      real thetawire(lx1,ly1,lz1,lelt) !local wire angle

      logical dopin(0:lpins)
      logical ifoutiters

c     variables for hexcan modification
      logical dohexcan
      real alpha,xi1,xi2,fact

c     variables for wire trimming (not yet supported)
      logical dotrim,perhex
      real cpinangles(5)
      real spinangles(5)
      real xtrim,ytrim,dtrim,rtrim,rnew
      integer iang

c     data initialization

      data dopin /.false.,lpins*.true./

      if(nrings.gt.10) then 
        if(nio.eq.0) write(*,*) 
     &    "Error: too many pin layers. Increase lpins in make_wire_wrap"
        call exitt
      endif

      n = lx1*ly1*lz1*nelt
      niters = 20 !number of iterations for deforming the mesh
      ifoutiters = .false. !output mesh at every iteration
      pio3 = pi/3.
      pio6 = pi/6.

      call izero(pindx,n)
      call rzero(delx,n)
      call rzero(dely,n)
      call rzero(delz,n)

c     trim angles for corner pins
      cpinangles(1) = pio6
      cpinangles(2) = 2*pio3
      cpinangles(3) = pi
      cpinangles(4) = 4*pio3
      cpinangles(5) = 11.*pio6

c     trim angles for side pins
      spinangles(1) = pio6
      spinangles(2) = 2.*pio3
      spinangles(3) = pi
      spinangles(4) = 4.*pio3
      spinangles(5) = 5.*pio3

      dotrim = .false. !trim the wires (not supported)
      perhex = .false. !single pin in fully periodic domain
      dohexcan = .false. !modify the hexcan
      if(Rcan.gt.0.) dohexcan = .true.

c     skip adding wire to individual pins
c     dopin(2)= .false.
c     dopin(4)= .false.
c     dopin(6)= .false.

      call domain_size(xmn,xmx,ymn,ymx,zmn,zmx)
      call pincenters(xxc,yyc,nrings,Pptch)

      nturns =  (zmx-zmn)/Wptch   !negative for left-handed wires
      call copy(thetawire,zm1,n)
      call rescale_x(thetawire,0.0,nturns*2.*pi)
      apoth = (ymx-ymn)/2.    !hexcan apothem

      Rp = 0.5 !pin radius
      Rw = Dwire*0.5
      Rf = Rf0
      if(Rf0.lt.0.) Rf = -Rf0*Rw  !filet radius
      sunk = sunk0
      if(sunk0.lt.0.) sunk = -sunk0*Rw !sink the wire into the pin
c     dtrim = Rp+(Pptch-1.0)-0.05*2.*Rw !not supported

C     what comes next is particularly esoteric...
c     distances between the centerpoints of pin, wire, and filet circles
      d1 = Rp+Rw-sunk      !pin to wire
      d2 = Rp+Rf           !pin to filet
      d3 = Rf+Rw           !filet to wire

c     angles between the centerpoints of pin, wire, and filet
      theta1 = acos((d2**2+d3**2-d1**2)/(2.0*d2*d3)) ! filet angle, opposite d1
      theta2 = acos((d1**2+d3**2-d2**2)/(2.0*d1*d3)) ! wire angle, opposite d2
      theta3 = pi - theta1 - theta2                  ! pin angle, opposite d3
c     if(nio.eq.0) write(6,*)
c    &                      theta1*180./pi,theta2*180./pi,theta3*180./pi

c     arc lengths of the transitions between circles, S4 < S3 < S6 < S5
c     middle of the wire is 0
      S4=Rw*(pi-theta2)       !transition from wire to filet on top
      S3=S4+Rf*theta1         !filet to pin
      S6=S3+2.*(pi-theta3)*Rp !pin to filet on bottom
      S5=S6+Rf*theta1         !filet to wire on bottom
      Stot=S5+S4              !total length of the perimeter

c     arc lengths for the first step
      psi = acos((Rp - Rw)/d1)
      delT = d1*sin(psi)
      T1  = Rw*psi
      T1_x = d1 + Rw*cos(psi)
      T1_y = Rw * sin(psi)
      T2  = T1 + delT
      T2_x = Rp*cos(psi)
      T2_y = Rp*sin(psi)
      Ttot= 2.*(T2 + (pi - psi)*Rp)
      T3  = Ttot - T2
      T4  = Ttot - T1

c     set the cbc array for mesh helmholtz solver
      do iel = 1,nelt
      do ifc = 1,2*ndim
        cbc(ifc,iel,0)=cbc(ifc,iel,1)
      enddo
      enddo
      call setbc(PBID,0,'mv ')
      call setbc(HCID,0,'mv ')

c     First, assign pin index and get the fractional arc length for all points on the pin surface
      do 010 iel = 1,nelt
      do 010 ifc = 1,2*ldim
        if(BoundaryID(ifc,iel).eq.PBID) then
          ipin = 0
          do 011 ilayer=1,nrings
          do 011 jpin = 1,max(1,6*(ilayer-1))
            ipin = ipin+1
            call get_face_m1centroid(xx,yy,zz,rr,iel,ifc)
            rr=sqrt((xx-xxc(ipin))**2+(yy-yyc(ipin))**2) !need radius to pin center, not origin
            if (abs(rr-Rp).lt.5.e-2) then
              call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
              do 020 k=k0,k1
              do 020 j=j0,j1
              do 020 i=i0,i1
                pindx(i,j,k,iel) = ipin
                if(dopin(ipin)) then
                  xpt=xm1(i,j,k,iel)-xxc(ipin)
                  ypt=ym1(i,j,k,iel)-yyc(ipin)
                  thw = mod(thetawire(i,j,k,iel),2.*pi)
                  theta = atan2(ypt,xpt) - thw + 4.*pi
                  theta = mod(theta,2.*pi)
                  sptl = theta/(2.*pi)
                  spt(i,j,k,iel) = sptl
cc                 experimental code to remap points 
c                 sptl = theta/pi !0 < S < 2
c                 if(sptl.le.1.0) then !0 < S < 1
c                   sptl=(Y0-1.)*sptl*sptl+(2.-Y0)*sptl
c                 elseif(sptl.gt.1.0) then  ! 1 < S < 2
c                   sptl = 2.- sptl
c                   sptl = (Y0-1.)*sptl*sptl+(2.-Y0)*sptl
c                   sptl = 2.- sptl
c                 endif
c                 spt(i,j,k,iel) = sptl/2.0
                endif
 020          continue
            endif
 011      continue
        endif
 010  continue

c     step 1: add a wire with infinite filet radius
      do i=1,n
        ipin=pindx(i,1,1,1)
        if(dopin(ipin)) then !always skip ipin = 0
          xpt=xm1(i,1,1,1)-xxc(ipin)
          ypt=ym1(i,1,1,1)-yyc(ipin)
          thw = thetawire(i,1,1,1)
          sptl = spt(i,1,1,1) * Ttot !load from array
C         Determine boundary displacement
          if(sptl.le.T1) then
            thloc = sptl/T1*psi
            xnew = d1 + Rw*cos(thloc)
            ynew = Rw*sin(thloc)
          elseif(sptl.le.T2) then
            xnew = (sptl-T1)/delT*(T2_x-T1_x)+T1_x
            ynew = (sptl-T1)/delT*(T2_y-T1_y)+T1_y
          elseif(sptl.le.T3) then
            thloc = psi + (sptl-T2)/Rp
            xnew = Rp*cos(thloc)
            ynew = Rp*sin(thloc)
          elseif(sptl.le.T4) then
            xnew = (sptl-T3)/delT*(T1_x-T2_x)+T2_x
            ynew =((sptl-T3)/delT*(T1_y-T2_y)+T2_y)*(-1.)
          else
            thloc = (sptl-Ttot)/T1*psi
            xnew = d1 + Rw*cos(thloc)
            ynew = Rw*sin(thloc)
          endif
          call rotate_point_2d(xnew,ynew,0.0,0.0,thw,xr,yr)
          delx(i,1,1,1) = xr - xpt
          dely(i,1,1,1) = yr - ypt
        endif
      enddo

c     If you're going to modify the hexcan, I recommend you do it in step 1
      if(dohexcan) call rounded_displacement(delx,dely,delz,Rcan,HCID)

      call ww_mv_mesh(delx,dely,delz,niters,ifoutiters,'st1') !execute step 1

c     reset displacement arrays
      call rzero(delx,n)
      call rzero(dely,n)
      call rzero(delz,n)

c     step 2: add the filet
      do i=1,n
        ipin=pindx(i,1,1,1)
        if(dopin(ipin)) then
          xpt=xm1(i,1,1,1)-xxc(ipin)
          ypt=ym1(i,1,1,1)-yyc(ipin)
          thw = thetawire(i,1,1,1)
          sptl = spt(i,1,1,1) * Stot !load from array
          if(sptl.gt.S5) then !on the bottom of the wire
            thloc = (sptl-S5)/Rw+pi+theta2
            xnew=D1+Rw*cos(thloc)
            ynew=Rw*sin(thloc)
          elseif(sptl.gt.S6) then !on the bottom fillet
            thloc = (S5-sptl)/Rf+theta2
            xnew=  D2*cos(theta3)+Rf*cos(thloc)
            ynew= -D2*sin(theta3)+Rf*sin(thloc)
          elseif(sptl.gt.S3) then !on the pin
            thloc = (sptl-S3)/Rp+theta3
            xnew= Rp*cos(thloc)
            ynew= Rp*sin(thloc)
          elseif(sptl.gt.S4) then !on the top fillet
            thloc = (S3-sptl)/Rf+pi+theta3
            xnew= D2*cos(theta3)+Rf*cos(thloc)
            ynew= D2*sin(theta3)+Rf*sin(thloc)
          else  !on the top of the wire
            thloc =sptl/Rw
            xnew= D1+Rw*cos(thloc)
            ynew= Rw*sin(thloc)
          endif
          call rotate_point_2d(xnew,ynew,0.0,0.0,thw,xr,yr)
          delx(i,1,1,1) = xr - xpt
          dely(i,1,1,1) = yr - ypt
        endif
      enddo

      call ww_mv_mesh(delx,dely,delz,niters,ifoutiters,'st2') !execute step 2

      return
      end     
c-----------------------------------------------------------------------
      subroutine rounded_displacement(delx,dely,delz,RR,bid)

c     This subroutine calculates a displacement vector for a fractional
c     arc-length preserving transformation from a sharp-cornered hexagon
c     to a rounded corner hexagon. It assumes a sharp corner is bisected
c     by the x-axis.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      real delx(lx1,ly1,lz1,lelt)  !displacement vector
      real dely(lx1,ly1,lz1,lelt)
      real delz(lx1,ly1,lz1,lelt)
      real RR                      !radius of the rounded corner
      integer bid !boundary ID to modify

      integer ifc,iel,i0,i1,j0,j1,k0,k1,i,j,k,n,icorn
      real x0,xxc,x1,xx0,xx1,y0,yyc,y1,yy0,yy1,s0,s1,theta0,theta1
      real glmax,psi,ddx,ddy,delmax

      if(RR.le.0) return

      n = lx1*ly1*lz1*nelt

      x0 = glmax(xm1,n)      !coordinate of the sharp edge
      y0 = 0.0
      xxc = x0-RR/cos(pi/6.) !center of the rounded edge
      yyc = 0.0
      x1 = RR*cos(pi/6.)+xxc !tangent point between rounded edge and hexcan
      y1 = RR*sin(pi/6.)+yyc
      delmax = RR*(1/cos(pi/6.)-1)
c     if(nio.eq.0) write(*,*) "delmax = ",delmax

      do 10 iel=1,nelv
      do 10 ifc=1,2*ldim
        if(cbc(ifc,iel,0).eq.'mv '.and.BoundaryID(ifc,iel).eq.bid) then
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
          do 20 icorn = 0,5
            psi = real(icorn)*pi/3.
            xx0 = xm1(i,j,k,iel)
            yy0 = ym1(i,j,k,iel)
c           Rotate the reference frame so each corner is aligned with the x-axis
            call rotate_point_2d(xx0,yy0,0.0,0.0,-psi,xx0,yy0)
            theta0 = atan2(yy0,xx0)
            if(abs(theta0).lt.pi/6.) then
              s0 = (xx0-x0)/(x1-x0)
              if(s0.lt.1.0) then
                theta1 = pi/6.*s0
                xx1 = RR*cos(theta1)+xxc
                yy1 = RR*sin(theta1)+yyc
                if(theta0.lt.0.0) yy1=-yy1
                ddx = xx1-xx0
                ddy = yy1-yy1
c               Rotate displacement from aligned reference frame back to original
                call rotate_point_2d(ddx,ddy,0.0,0.0,psi,ddx,ddy)
                delx(i,j,k,iel)=ddx
                dely(i,j,k,iel)=ddy
              endif
            endif
  20      continue
        endif
  10  continue

      return
      end

c----------------------------------------------------------------------
      subroutine ww_mv_mesh(dxo,dyo,dzo,nstps,ifout,na3)

c     This subroutine solves for and applies the overall mesh 
c     displacement when provided with the displacment vector on boundaries
c     with the 'mv ' BC in field 0.

      include 'SIZE'
      include 'TOTAL'
      character*3 na3
      real umeshx(lx1,ly1,lz1,lelt),dxo(lx1,ly1,lz1,lelt)
      real umeshy(lx1,ly1,lz1,lelt),dyo(lx1,ly1,lz1,lelt)
      real umeshz(lx1,ly1,lz1,lelt),dzo(lx1,ly1,lz1,lelt)
      parameter (lt = lx1*ly1*lz1*lelt)
      common /mrthoi/ napprx(2),nappry(2),napprz(2)
      common /mrthov/ apprx(lt,0:mxprev)
     $              , appry(lt,0:mxprev)
     $              , apprz(lt,0:mxprev)
      common /mstuff/ d(lt),h1(lt),h2(lt),mask(lt)
      real mask,pmax,pmin
      real srfbl,volbl,delta,deltap1,deltap2,arg1,arg2
      real zero,one
      integer e,f,nstps,ifield_sv,nbl
      integer icalld
      logical ifxyos,ifout
      save    icalld
      data    icalld /0/

c     use local arrays to avoid relying on lx1m
      real delx(lx1,ly1,lz1,lelt)
      real dely(lx1,ly1,lz1,lelt)
      real delz(lx1,ly1,lz1,lelt)

      real finx(lx1,ly1,lz1,lelt)
      real finy(lx1,ly1,lz1,lelt)
      real finz(lx1,ly1,lz1,lelt)
      real relax

      ifield_sv=ifield
      ifield = 0
      restol(ifield) = 1.0e-3 !solver tolerance

      n = nx1*ny1*nz1*nelv
      nface = 2*ndim
      zero = 0.
      one  = 1.

c     For non-linear relaxation
      call copy(finx,xm1,n)
      call copy(finy,ym1,n)
      call copy(finz,zm1,n)

      call add2(finx,dxo,n)
      call add2(finy,dyo,n)
      call add2(finz,dzo,n)

c     For linear relaxation
      relax = 1.0/real(nstps)
      call copy(umeshx,dxo,n)
      call copy(umeshy,dyo,n)
      call copy(umeshz,dzo,n)
      call cmult(umeshx,relax,n)
      call cmult(umeshy,relax,n)
      call cmult(umeshz,relax,n)

      utx_usr=glamax(dxo,n)
      uty_usr=glamax(dyo,n)  
      utz_usr=glamax(dzo,n)

      if (nid.eq.0) then
        write(6,*) "utx_usr: ",utx_usr
        write(6,*) "uty_usr: ",uty_usr
        write(6,*) "utz_usr: ",utz_usr
      endif

      time = 0.0    
      call save_ioflags
      call clear_ioflags
      ifxyo = .true.
      ifvo = .true.
      if(nstps.gt.1.and.ifout) then !never have 2 files because of Paraview...
        call copy(vx,dxo,n)
        call copy(vy,dyo,n)
        call copy(vz,dzo,n)
        call prepost(.true.,na3)
      endif
      nbl = 0

      do istep = 1,nstps
c       factor for decaying boundary layer preservation
c       fact = real(istep-1)/real(nstps-1)
c       fact = 1.0-fact 
        fact = 1.0

c       for non-linear relaxation
c       call sub3(umeshx,finx,xm1,n)
c       call sub3(umeshy,finy,ym1,n)
c       call sub3(umeshz,finz,zm1,n)
c       relax = 1.0 - (real(istep)/real(nstps)-1.0)**2
cc      relax = real(istep)/real(nstps)
c       if(nio.eq.0) write(*,*) "relaxation factor: ",relax
c       call cmult(umeshx,relax,n)
c       call cmult(umeshy,relax,n)
c       call cmult(umeshz,relax,n)

        napprx(1)=0
        nappry(1)=0
        napprz(1)=0
        nxz   = nx1*nz1
        nxyz  = nx1*ny1*nz1

        if (icalld.eq.0) then
          icalld=1
          call rone(mask,n)
          do e=1,nelv
          do f=1,nface
c           if(cbc(f,e,0).eq.'W  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
c           if(cbc(f,e,0).eq.'W1 ')call facev(mask,e,f,one ,nx1,ny1,nz1)
c           if(cbc(f,e,0).eq.'v  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
c           if(cbc(f,e,0).eq.'O  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'mv ')call facev(mask,e,f,zero,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'mvb')then
              call facev(mask,e,f,zero,nx1,ny1,nz1)
              nbl = 1
            endif
          enddo
          enddo
          call dsop(mask,'*  ',nx1,ny1,nz1)    ! dsop mask
          call opzero(delx,dely,delz)
          nbl = iglsum(nbl,1)
        endif ! icalld

        call rone (h1,n)
        call rzero(h2,n)

        if(nbl.gt.0) then
          srfbl = 0.   ! Surface area of elements in b.l.
          volbl = 0.   ! Volume of elements in boundary layer
          do e=1,nelv
          do f=1,nface
            if (cbc(f,e,0).eq.'mvb') then
              srfbl = srfbl + vlsum(area(1,1,f,e),nxz )
              volbl = volbl + vlsum(bm1 (1,1,1,e),nxyz)
            endif
          enddo
          enddo
          srfbl = glsum(srfbl,1)  ! Sum over all processors
          volbl = glsum(volbl,1)
       
          call cheap_dist_0(d,0,'mvb')

          delta=volbl/srfbl
          if (nid.eq.0) write(6,*) "delta: ",delta
          deltap1 = 2.0*delta! /real(lx1)  ! Protected b.l. thickness
          deltap2 = 2.0*delta

c         magic distribution - it really does a better job of preseving BLs 
          do i=1,n
            arg1   = -(d(i)/deltap1)**2
            arg2   = -(d(i)/deltap2)
            h1(i)  = h1(i) + 
     &                  fact*(5.0*exp(arg1) +1.0*exp(arg2))
          enddo
        endif

        do e=1,nelv
        do f=1,nface
          if (cbc(f,e,0).eq.'mv '.or.cbc(f,e,0).eq.'mvb') then
           call facec(delx,umeshx,e,f,nx1,ny1,nz1,nelv)
           call facec(dely,umeshy,e,f,nx1,ny1,nz1,nelv)
           call facec(delz,umeshz,e,f,nx1,ny1,nz1,nelv)
          endif
        enddo
        enddo
        tol = 1.e-6 

        if (utx_usr.gt.1e-10)
     &    call laplaceww('mshx',delx,h1,h2,mask,vmult,1,tol,
     &    1000,apprx,napprx)
        if (uty_usr.gt.1e-10) 
     &    call laplaceww('mshy',dely,h1,h2,mask,vmult,1,tol,
     &    1000,appry,nappry)
        if (utz_usr.gt.1e-10)
     &    call laplaceww('mshz',delz,h1,h2,mask,vmult,1,tol,
     &    1000,apprz,napprz)

        call dsavg(delx)
        call dsavg(dely)
        call dsavg(delz)

        call add2(xm1,delx,n)
        call add2(ym1,dely,n)
        call add2(zm1,delz,n)

        call copy(vx,delx,n)
        call copy(vy,dely,n)
        call copy(vz,delz,n)

        time = istep
        if(ifout) call prepost(.true.,na3)

        call fix_geom
        call mesh_metrics(.true.)

        djmin = vlmin(JACM1,n)
        djmax = vlmax(JACM1,n)

c       djave = glsum(JACM1,n)/real(nelgv*lx1*ly1*lz1)
c        do i=1,n
c         djave = djave + JACM1(i,1,1,1)
c       enddo
c       djave = glsum(djave,1)/real(nelgv*lx1*ly1*lz1)

        if(nio.eq.0) then 
          write(6,'(A,1p2E9.2)') ' Absolute Jacobian  min/max:',
     &      djmin,djmax
          write(6,*)
        endif

      enddo

      call restore_ioflags
      ifield = ifield_sv

      return
      end
c-----------------------------------------------------------------------
      subroutine laplaceww
     $     (name,u,h1,h2,mask,mult,ifld,tli,maxi,approx,napprox)
c
c     Solve Laplace's equation, with projection onto previous solutions.
c
c     Boundary condition strategy:
c
c     u = u0 + ub
c
c        u0 = 0 on Dirichlet boundaries
c        ub = u on Dirichlet boundaries
c
c        _
c        A ( u0 + ub ) = 0
c
c        _            _
c        A  u0  =   - A ub
c
c        _             _
c       MAM u0  =   -M A ub,    M is the mask
c
c                      _
c        A  u0  =   -M A ub ,  Helmholtz solve with SPD matrix A
c
c        u = u0+ub
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
c
      character*4 name
      real u(1),h1(1),h2(1),mask(1),mult(1),approx (1)
      integer   napprox(1)

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scruz/ r (lt),ub(lt)

      logical ifstdh
      character*4  cname
      character*6  name6

      logical ifwt,ifvec

      call chcopy(cname,name,4)
      call capit (cname,4)

      call blank (name6,6)
      call chcopy(name6,name,4)
      ifwt  = .true.
      ifvec = .false.
      isd   = 1
      imsh  = 1
      nel   = nelfld(ifld)

      n = nx1*ny1*nz1*nel

      call copy (ub,u,n)             ! ub = u on boundary
      call dsavg(ub)                 ! Make certain ub is in H1
                                     !     _
      call axhelm (r,ub,h1,h2,1,1)   ! r = A*ub

      do i=1,n                       !        _
         r(i)=-r(i)*mask(i)          ! r = -M*A*ub
      enddo

      call dssum  (r,nx1,ny1,nz1)    ! dssum rhs

c     call project1
c    $    (r,n,approx,napprox,h1,h2,mask,mult,ifwt,ifvec,name6)

      tol = abs(tli)
      p22=param(22)
      param(22)=abs(tol)
      restol(ifield)=tol
      if (nel.eq.nelv) then
        call hmhzpf (name,u,r,h1,h2,mask,mult,imsh,tol,maxi,isd,binvm1)
      else
        call hmhzpf (name,u,r,h1,h2,mask,mult,imsh,tol,maxi,isd,bintm1)
      endif
      param(22)=p22

c     call project2
c    $     (u,n,approx,napprox,h1,h2,mask,mult,ifwt,ifvec,name6)

      call add2(u,ub,n)

      return
      end
C-----------------------------------------------------------------------
      subroutine cheap_dist_0(d,ifld,b)

c     Finds a pseudo-distance function.

c     INPUT:  ifld - field type for which distance function is to be found.
c             ifld = 1 for velocity
c             ifld = 2 for temperature, etc.

c     OUTPUT: d = "path" distance to nearest wall

c     This approach has a significant advantage that it works for
c     periodict boundary conditions, whereas most other approaches
c     will not.

      include 'SIZE'
      include 'GEOM'       ! Coordinates
      include 'INPUT'      ! cbc()
      include 'TSTEP'      ! nelfld
      include 'PARALLEL'   ! gather-scatter handle for field "ifld"

      real d(lx1,ly1,lz1,lelt)

      character*3 b  ! Boundary condition of interest

      integer e,eg,f

      nel = nelt
      n = lx1*ly1*lz1*nel

      call domain_size(xmin,xmax,ymin,ymax,zmin,zmax)

      xmn = min(xmin,ymin)
      xmx = max(xmax,ymax)
      if (if3d) xmn = min(xmn ,zmin)
      if (if3d) xmx = max(xmx ,zmax)

      big = 10*(xmx-xmn)
      call cfill(d,big,n)

      nface = 2*ldim
      do e=1,nel     ! Set d=0 on walls
      do f=1,nface
        if (cbc(f,e,ifld).eq.b) call facev(d,e,f,0.,lx1,ly1,lz1)
      enddo
      enddo

      do ipass=1,10000
         dmax    = 0
         nchange = 0
         do e=1,nel
           do k=1,lz1
           do j=1,ly1
           do i=1,lx1
             i0=max(  1,i-1)
             j0=max(  1,j-1)
             k0=max(  1,k-1)
             i1=min(lx1,i+1)
             j1=min(ly1,j+1)
             k1=min(lz1,k+1)
             do kk=k0,k1
             do jj=j0,j1
             do ii=i0,i1

              if (if3d) then
               dtmp = d(ii,jj,kk,e) + dist3d(
     $           xm1(ii,jj,kk,e),ym1(ii,jj,kk,e),zm1(ii,jj,kk,e)
     $          ,xm1(i ,j ,k ,e),ym1(i ,j ,k ,e),zm1(i ,j ,k ,e))
              else
               dtmp = d(ii,jj,kk,e) + dist2d(
     $           xm1(ii,jj,kk,e),ym1(ii,jj,kk,e)
     $          ,xm1(i ,j ,k ,e),ym1(i ,j ,k ,e))
              endif

              if (dtmp.lt.d(i,j,k,e)) then
                d(i,j,k,e) = dtmp
                nchange = nchange+1
                dmax = max(dmax,d(i,j,k,e))
              endif
             enddo
             enddo
             enddo

           enddo
           enddo
           enddo
         enddo
         call fgslib_gs_op(gsh_fld(ifld),d,1,3,0) ! min over all elements
         nchange = iglsum(nchange,1)
         dmax = glmax(dmax,1)
         if (nio.eq.0.and.loglevel.gt.2) write(6,1) ipass,nchange,dmax,b
    1    format(i9,i12,1pe12.4,' max distance b: ',a3)
         if (nchange.eq.0) goto 1000
      enddo
 1000 return
      end

c-----------------------------------------------------------------------
C Old code that might come in handy one day
c-----------------------------------------------------------------------

c     old indexing algorithm with wire trimming
c      do 110 iel = 1,nelt
c      do 110 ifc = 1,2*ldim
c        if(BoundaryID(ifc,iel).eq.1) then  !make the wire-wraps
c          ipin = 0
c          do 111 ilayer=1,nrings
c          do 111 jpin = 1,max(1,6*(ilayer-1))
c            ipin = ipin+1
c            if(dopin(ipin)) then
c              call get_face_m1centroid(xx,yy,zz,rr,iel,ifc)
c              rr=sqrt((xx-xxc(ipin))**2+(yy-yyc(ipin))**2)
c              if (abs(rr-0.5).lt.5.e-2) then
c                call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
c                do 120 k=k0,k1
c                do 120 j=j0,j1
c                do 120 i=i0,i1
c                  xpt=xm1(i,j,k,iel)-xxc(ipin)
c                  ypt=ym1(i,j,k,iel)-yyc(ipin)
c                  thw = tmax*zm1(i,j,k,iel)
c                  sptl = spt(i,j,k,iel) * Ttot !load from array
cC                 Determine boundary displacement
c                  if(sptl.le.T1) then
c                    thloc = sptl/T1*psi1
c                    xnew = d1 + Rw*cos(thloc)
c                    ynew = Rw*sin(thloc)
c                  elseif(sptl.le.T2) then
c                    xnew = (sptl-T1)/delT*(T2_x-T1_x)+T1_x
c                    ynew = (sptl-T1)/delT*(T2_y-T1_y)+T1_y
c                  elseif(sptl.le.T3) then
c                    thloc = psi1 + (sptl-T2)/Rp
c                    xnew = Rp*cos(thloc)
c                    ynew = Rp*sin(thloc)
c                  elseif(sptl.le.T4) then
c                    xnew = (sptl-T3)/delT*(T1_x-T2_x)+T2_x
c                    ynew =((sptl-T3)/delT*(T1_y-T2_y)+T2_y)*(-1.)
c                  else
c                    thloc = (sptl-Ttot)/T1*psi1
c                    xnew = d1 + Rw*cos(thloc)
c                    ynew = Rw*sin(thloc)
c                  endif
c                  call rotate_point_2d(xnew,ynew,0.0,0.0,thw,xr,yr)
c
cc   Trim the tips of the wires, if necessary, do this last??
c                  if(dotrim) then
c                    if((nrings.eq.1).and.(.not.perhex)) then !only 1 pin, special case, UNTESTED!!
c                      theta=atan2(yr,xr)
c                      rnew=sqrt(xr*xr+yr*yr)
c                      if(theta.lt.0.0) theta = theta + 2.*pi
c                      alpha=mod(theta,psi)-psi/2.
c                      rtrim = dtrim/cos(alpha)
c                      rnew=min(rnew,rtrim)
c                      xr=rnew*cos(theta)
c                      yr=rnew*sin(theta)
c                    elseif((ilayer.lt.nrings).or.perhex) then !internal layers do all 6 angles
c                      theta=atan2(yr,xr)
c                      rnew=sqrt(xr*xr+yr*yr)
c                      if(theta.lt.0.0) theta = theta + 2.*pi
c                      alpha=mod(theta+psi/2.,psi)-psi/2.
c                      rtrim = dtrim/cos(alpha)
c                      rnew=min(rnew,rtrim)
c                      xr=rnew*cos(theta)
c                      yr=rnew*sin(theta)
c                    elseif(nrings.eq.2.or.
c     &                            mod(jpin,max(1,nrings-1)).eq.1) then !corner pin
c                      theta=atan2(yr,xr)
c                      rnew=sqrt(xr*xr+yr*yr)
c                      if(theta.lt.0.0) theta = theta + 2.*pi
c                      alpha=2.*pi
c                      do iang=1,5
c                        thw=theta-floor(real(jpin-1)/real(nrings-1))
c     &                                                           *psi
c                        if(thw.lt.0.0) thw=thw+2.*pi
c                        alpha=min(alpha,abs(thw-cpinangles(iang)))
c                      enddo
c                      rtrim = dtrim/cos(alpha)
c                      rnew=min(rnew,rtrim)
c                      xr=rnew*cos(theta)
c                      yr=rnew*sin(theta)
c                    else !side pins
c                      theta=atan2(yr,xr)
c                      rnew=sqrt(xr*xr+yr*yr)
c                      if(theta.lt.0.0) theta = theta + 2.*pi
c                      alpha=2.*pi
c                      do iang=1,5
c                        thw=theta-floor(real(jpin-1)/real(nrings-1))
c     &                                                         *psi
c                        if(thw.lt.0.0) thw=thw+2.*pi
c                        alpha=min(alpha,abs(thw-spinangles(iang)))
c                      enddo
c                      rtrim = dtrim/cos(alpha)
c                      rnew=min(rnew,rtrim)
c                      xr=rnew*cos(theta)
c                      yr=rnew*sin(theta)
c                    endif
c                  endif
c                  delx(i,j,k,iel) = xr - xpt
c                  dely(i,j,k,iel) = yr - ypt
c 120            continue
c              endif
c            endif
c 111      continue
c        endif
c 110  continue

c     morph the hexcan if necessary
c     if(dohexcan) then
c       do 310 iel = 1,nelt
c       do 310 ifc = 1,2*ldim
c         if(BoundaryID(ifc,iel).eq.2) then  !modify the hex-can
c           call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
c           do 321 k=k0,k1
c           do 321 j=j0,j1
c           do 321 i=i0,i1
c             xpt=xm1(i,j,k,iel)
c             ypt=ym1(i,j,k,iel)
c             theta=atan2(ypt,xpt)
c             if(theta.lt.0.0) theta = theta + 2.*pi
c             alpha = theta - mod(theta,pio3)
c             call rotate_point_2d(xpt,ypt,0.0,0.0,-alpha,xnew,ynew)
c             sptl = ((xnew-xi1)/(xi2-xi1)+alpha/pio3)/6.
c             theta = 2.*pi*sptl
c             zpt = zm1(i,j,k,iel)/(2.*pi)
c             fact = 1.-(2.*zpt-1.)**2
c             delx(i,j,k,iel) = (xi1*cos(theta) - xpt)*fact*0.65
c             dely(i,j,k,iel) = (xi1*sin(theta) - ypt)*fact*0.65
c321        continue
c         endif
c310    continue
c     endif

c     calculate effective conductivity for mesh solve to preserve BL in the filet, but NOT on the wire
c     thcr=1.1*atan(Rw/(Rp))
c     wdth=4.0
c     do i=1,n
c       xpt = xm1(i,1,1,1)
c       ypt = ym1(i,1,1,1)
c       thw = thetawire(i,1,1,1)
c       thp = mod(thw+pi/6.,pi/3.)-pi/6.
c       thb = thw-thp
c       if(thb.gt.pi) then
c         tho = thb-pi
c       else
c         tho = thb+pi
c       endif
c       tho = tho - thp !theta of the opposing wire
c       theta = atan2(ypt,xpt)
c       thw = thw-theta
c       tho = tho-theta
c       if(thw.gt.pi) thw = thw-2.*pi
c       if(thw.lt.-pi) thw = thw+2.*pi
c       if(tho.gt.pi) tho = tho-2.*pi
c       if(tho.lt.-pi) tho = tho+2.*pi
c       bl(i)=0.25
c    &   *(tanh(wdth*(abs(thw)/thcr-1.))+1.0)
c    &   *(tanh(wdth*(abs(tho)/thcr-1.))+1.0)
c     enddo
