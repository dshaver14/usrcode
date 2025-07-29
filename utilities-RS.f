c-----------------------------------------------------------------------
      real function planar_ave_m1(phi,norm,pt,eps)
      implicit none
C
C     Compute area average of phi() on the
C     plane defined by normal 'norm' and point 'pt'
C
      include 'SIZE'
      include 'TOTAL'

      real phi(lx1*ly1*lz1,1),norm(3),pt(3),eps,phi_ext
      real dpdx(lx1*ly1*lz1)
      real dpdy(lx1*ly1*lz1)
      real dpdz(lx1*ly1*lz1)
      real aa,bb,cc,dd,w1,w2,x0,y0,z0,r0,rr,del,xsa,glsum

      integer e,i,j,k,n

      n=lx1*ly1*lz1

!     make sure norm is a unit vector
      dd=(norm(1)**2+norm(2)**2+norm(3)**2)
      if(dd.gt.0) dd=sqrt(dd)
      aa=norm(1)/dd
      bb=norm(2)/dd
      cc=0.0
      if(if3d) cc=norm(3)/dd
      dd=-1.0*(aa*pt(1)+bb*pt(2)+cc*pt(3))

      w1=0.0
      w2=0.0
      do e=1,nelv
        call gradm11(dpdx,dpdy,dpdz,phi,e)
        do i=1,n
          x0=xm1(i,1,1,e)
          y0=ym1(i,1,1,e)
          z0=zm1(i,1,1,e)
          r0=aa*(x0-pt(1))+bb*(y0-pt(2))+cc*(z0-pt(3)) !signed distance to plane
          rr=min(2.0,abs(r0)*2.0/eps)
          phi_ext=phi(i,e)+r0*(aa*dpdx(i)+bb*dpdy(i)+cc*dpdz(i)) !1st order extrapolation to the plane
          if(rr.gt.1.0) then
            del = 1.0/8.0*(5.0-2.0*rr-sqrt(-7.0+12.0*rr-4.0*rr**2))
          else 
            del = 1.0/8.0*(3.0-2.0*rr+sqrt( 1.0+ 4.0*rr-4.0*rr**2))
          endif
          w1=w1+phi_ext*bm1(i,1,1,e)*del
          w2=w2+bm1(i,1,1,e)*del
        enddo
      enddo
      xsa=glsum(w2,1)
      planar_ave_m1 = glsum(w1,1)/max(xsa,1.0e-8)
      xsa=2.0*xsa/eps  !cross sectional area

      return
      end
C-----------------------------------------------------------------------
      real function planar_ave_m2(phi,norm,pt,eps)
      implicit none
C
C     Compute area average of phi() on the
C     plane defined by normal 'norm' and point 'pt'
C
      include 'SIZE'
      include 'TOTAL'

      real phi(1),norm(3),pt(3),eps
      real aa,bb,cc,dd,w1,w2,x0,y0,z0,r0,rr,del,xsa,glsum

      integer i,j,k,n

      n=lx2*ly2*lz2*nelv

      aa=norm(1)
      bb=norm(2) 
      cc=0.0
      if(if3d) cc=norm(3)
      dd=-1.0*(aa*pt(1)+bb*pt(2)+cc*pt(3))
      w1=0.0
      w2=0.0
      do i=1,n
        x0=xm2(i,1,1,1)
        y0=ym2(i,1,1,1)
        z0=zm2(i,1,1,1)
        r0=(aa*x0+bb*y0+cc*z0+dd)/sqrt(aa**2+bb**2+cc**2)
        rr=min(2.0,abs(r0)*2.0/eps)
        if(rr.gt.1.0) then
          del = 1.0/8.0*(5.0-2.0*rr-sqrt(-7.0+12.0*rr-4.0*rr**2))
        else 
          del = 1.0/8.0*(3.0-2.0*rr+sqrt( 1.0+ 4.0*rr-4.0*rr**2))
        endif
        w1=w1+phi(i)*bm2(i,1,1,1)*del
        w2=w2+bm2(i,1,1,1)*del
      enddo
      xsa=glsum(w2,1)
      planar_ave_m2 = glsum(w1,1)/max(xsa,1.0e-8)
      xsa=2.0*xsa/eps  !cross sectional area

      return
      end
C-----------------------------------------------------------------------
      subroutine div_check(phi)

      real phi(1)

      if(phi(1).ne.phi(1)) call exitt
      return
      end
C-----------------------------------------------------------------------
      subroutine get_wall_distance(wd,itype)
      implicit none

      include 'SIZE'

      integer icalled
      data icalled /0/
      save icalled

      real wd(*)
      real w1(lx1*ly1*lz1*lelv)
      real w2(lx1*ly1*lz1*lelv)
      real w3(lx1*ly1*lz1*lelv)
      real w4(lx1*ly1*lz1*lelv)
      real w5(lx1*ly1*lz1*lelv)
      real w6(lx1*ly1*lz1*lelv)
      common /SCRNS/ w1,w2,w3,w4,w5,w6

      integer n,itype

      if(icalled.eq.0) then
        if(itype.eq.1) then
          call cheap_dist(w1,1,'W  ')
        elseif(itype.eq.2) then
          call distf(w1,1,'W  ',w2,w3,w4,w5,w6)
        else
          if(nio.eq.0) write(*,*) 
     &           "Error in get_wall_distance, unsupported distance type"
        endif
        icalled = 1
      endif

      n=lx1*ly1*lz1*nelv
      call copy(wd,w1,n)

      return
      end
C-----------------------------------------------------------------------
      subroutine set_wwpin_BCs(Nlay,pitch)
      include 'SIZE'
      include 'TOTAL'

      real xxc(271),yyc(271)

      logical iflag

      radius = 0.5
      call pincenters(xxc,yyc,Nlay,pitch)

      do 10 ie=1,nelv
      do 10 ic=1,2*ldim
        iflag = .false.
        cbc(ic,ie,2)=cbc(ic,ie,1)
        if(cbc(ic,ie,1).eq.'W  ') then
          cbc(ic,ie,2)='I  '
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ic)
          do 20 ipin=1,Npin
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            xx=xm1(i,j,k,ie)-xxc(ipin)
            yy=ym1(i,j,k,ie)-yyc(ipin)
            dist=sqrt((xx)**2+(yy)**2)
            if(abs(dist-radius).le.1.0e-3) iflag=.true.
  20      continue
        endif
        if(iflag) cbc(ic,ie,2)='f  '
  10  continue

      return
      end
C-----------------------------------------------------------------------
      subroutine set_wwpin_BCids(Nlay,pitch,iid,idn)
      include 'SIZE'
      include 'TOTAL'
      
      real xxc(271),yyc(271)

      logical iflag

      radius = 0.5
      call pincenters(xxc,yyc,Nlay,pitch)
 
      do 10 ie=1,nelv
      do 10 ic=1,2*ldim
        iflag = .false.
        if(BoundaryID(ic,ie).eq.iid) then
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ic)
          do 20 ipin=1,Npin
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            xx=xm1(i,j,k,ie)-xxc(ipin)
            yy=ym1(i,j,k,ie)-yyc(ipin)
            dist=sqrt((xx)**2+(yy)**2)
            if(abs(dist-radius).le.1.0e-3) iflag=.true.
  20      continue
        endif
        if(iflag) BoundaryID(ic,ie)=idn
  10  continue

      return
      end
C-----------------------------------------------------------------------
      subroutine pincenters(xxc,yyc,Nlay,pitch)

      real xxc(1),yyc(1),pitch,tht,xx,yy,pi
      integer nlayers,ipin,Npin,ilay,j,k

      pi = 4.0*atan(1.0)

      ipin=1
      Npin=1
      xxc(1)=0.0
      yyc(1)=0.0
      do ilay=1,Nlay
        Npin=Npin+6*(ilay-1)
        if(ilay.gt.1) then
          do j= 1,6
            tht = (j-1)*pi/3.
            do k= 1,(ilay-1)
              ipin=ipin+1
              xx=(ilay-1)*pitch-(k-1)*pitch*cos(pi/3.)
              yy=(k-1)*pitch*sin(pi/3.)
              xxc(ipin)= xx*cos(tht)-yy*sin(tht)
              yyc(ipin)= xx*sin(tht)+yy*cos(tht)
            enddo
          enddo
        endif
      enddo

  255 format(a4,2a11)
  256 format(i4,2f11.6)

      return
      end
C-----------------------------------------------------------------------
      real function get_nearest(loc,coord)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'

      integer ipoint,ierr
      real loc,coord(1),ds(2),dsg(2)

      ds(1)=1.0d30
      do ipoint=1,lx1*ly1*lz1*nelv
        if(abs(coord(ipoint)-loc).lt.ds(1)) then
          ds(1)=abs(coord(ipoint)-loc)
          ds(2)=coord(ipoint)
        endif
      enddo
      call MPI_ALLREDUCE(ds,dsg,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC
     &                                             ,MPI_COMM_WORLD,ierr)
      get_nearest=dsg(2)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_point3d(loc1,loc2,loc3,ix,iy,iz,eg)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'

      integer ipoint,ierr,dsi,ix,iy,iz,ie,eg,jx,jy,jz,je
      real loc1,loc2,loc3,ds(2),dsg(2),dist(2)

      dist(2)=1.0d30
      do je=1,nelv
      do jz=1,lz1
      do jy=1,ly1
      do jx=1,lx1
        dist(1)=sqrt((loc1-xm1(jx,jy,jz,je))**2
     &           +(loc2-ym1(jx,jy,jz,je))**2+(loc3-zm1(jx,jy,jz,je))**2)
        if(dist(1).lt.dist(2)) then
          dist(2)=dist(1)
          ix=jx
          iy=jy
          iz=jz
          ie=je
        endif
      enddo
      enddo
      enddo
      enddo

      eg=lglel(ie)

      ds(1)=dist(2)

      ds(2)=dble(ix)
      call MPI_ALLREDUCE(ds,dsg,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC
     &                                             ,MPI_COMM_WORLD,ierr)
      ix=int(dsg(2))

      ds(2)=dble(iy)
      call MPI_ALLREDUCE(ds,dsg,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC
     &                                             ,MPI_COMM_WORLD,ierr)
      iy=int(dsg(2))

      ds(2)=dble(iz)
      call MPI_ALLREDUCE(ds,dsg,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC
     &                                             ,MPI_COMM_WORLD,ierr)
      iz=int(dsg(2))

      ds(2)=dble(eg)
      call MPI_ALLREDUCE(ds,dsg,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC
     &                                             ,MPI_COMM_WORLD,ierr)
      eg=int(dsg(2))

      return
      end
c-----------------------------------------------------------------------
      real function get_nearest_face(loc,coord,norm)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'

      integer ielem,iside,i0,i1,j0,j1,k0,k1,i,j,k,ierr
      real loc,coord(1),ds(2),dsg(2),norm(3),fnorm(3),dp

      ds(1)=1.0d30
      do ielem=1,nelv
        do iside=1,ldim*2
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,iside)
          i=(i0+i1)/2
          j=(j0+j1)/2
          k=(k0+k1)/2
          ipoint=i+(j-1)*lx1+(k-1)*lx1*ly1+lx1*ly1*lz1*(ielem-1)
          call getSnormal(fnorm,i,j,k,iside,ielem)
          dp=fnorm(1)*norm(1)+fnorm(2)*norm(2)
          if(if3d) dp=dp+fnorm(3)*norm(3)
          if((1.0d0-abs(dp)).lt.1.0d-8)then
            if(abs(coord(ipoint)-loc).lt.ds(1)) then
              ds(1)=abs(coord(ipoint)-loc)
              ds(2)=coord(ipoint)
            endif
          endif
        enddo
      enddo
      call MPI_ALLREDUCE(ds,dsg,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC
     &                                             ,MPI_COMM_WORLD,ierr)
      get_nearest_face=dsg(2)
      return
      end
C-----------------------------------------------------------------------
      subroutine weighted_average(phi,wrt,loc,coord,norm,phia)
C
C     Compute planar averages of phi() weighted by wrt() on the
C     plane normal to norm() with intercept coord = loc
C
      include 'SIZE'
      include 'TOTAL'

      integer ielem,iside,i,i0,i1,j,j0,j1,k,k0,k1
      real phi(1),wrt(1),loc,coord(1),norm(3),phia
      real fnorm(3),dp,a1,phia1

      loc=get_nearest_face(loc,coord,norm)

      phia=0.0
      phia1=0.0
      do ielem=1,nelv
        do iside=1,ndim*2
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,iside)
          i=(i0+i1)/2  !just use the point in the middle of the face
          j=(j0+j1)/2
          k=(k0+k1)/2
          ipoint=i+(j-1)*lx1+(k-1)*lx1*ly1+(ielem-1)*lx1*ly1*lz1
          if(abs(coord(ipoint)-loc).lt.1.0d-8)then 
            call getSnormal(fnorm,i,j,k,iside,ielem)
            dp=fnorm(1)*norm(1)+fnorm(2)*norm(2)
            if(if3d) dp=dp+fnorm(3)*norm(3)
            if(abs(1.0d0-dp).lt.1.0d-8) then
              do i=i0,i1
              do j=j0,j1
              do k=k0,k1
                ipoint=i+(j-1)*lx1+(k-1)*lx1*ly1+(ielem-1)*lx1*ly1*lz1
                if    ((iside.eq.1).or.(iside.eq.3)) then
                  a1=area(i,k,iside,ielem)
                elseif((iside.eq.2).or.(iside.eq.4)) then
                  a1=area(j,k,iside,ielem)
                else
                  a1=area(i,j,iside,ielem)
                endif
                phia=phia+phi(ipoint)*wrt(ipoint)*a1
                phia1=phia1+wrt(ipoint)*a1
              enddo
              enddo
              enddo
            endif
          endif
        enddo
      enddo

      phia=glsum(phia,1)
      phia1=glsum(phia1,1)
      phia=phia/phia1

      return
      end
C-----------------------------------------------------------------------
      subroutine planar_average_weighted(phia,phi,wrt,w1,w2)
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'ZPER'

      real phi(nx1*ny1,nz1,nelv),wrt(nx1*ny1,nz1,nelv),phia(nz1,nelz)
      real w1(nz1,nelz),w2(nz1,nelz) !work arrays

      integer e,eg,ez,melxy,nz,i,k
      real zz,aa

      melxy=nelx*nely !number of elements in the plane
      nz=nz1*nelz !number of z-slices

      if(melxy.lt.1)then
        if(nio.eq.0)write(*,256)'nelx*nely'
        return
      elseif(nelz.lt.1) then
        if(nio.eq.0)write(*,256)'nelz'
        return
      elseif(melxy.gt.lelx*lely) then
        if(nio.eq.0)write(*,257)'nelx*nely','lelx*lely'
        return
      elseif(nelz.gt.lelz) then
        if(nio.eq.0)write(*,257)'nelz','lelz'
        return
      endif

 256  format(5x,'ERROR: ',a,' must be at least 1!')
 257  format(5x,'ERROR: ',a,' must be less than ',a,'!')

      call rzero(phia,nz)
      call rzero(w1,nz)

      do e=1,nelt
        eg=lglel(e)
        ez=1+(eg-1)/melxy !z-slice id
        do k=1,nz1
          do i=1,nx1*ny1
            zz=(1.0-zgm1(k,3))/2.0
            aa=zz*area(i,1,5,e)+(1.0-zz)*area(i,1,6,e)
            w1(k,ez)=w1(k,ez)+aa*wrt(i,k,e)
            phia(k,ez)=phia(k,ez)+aa*wrt(i,k,e)*phi(i,k,e)
          enddo
          if(abs(w1(k,ez)).lt.1.0d-10) w1(k,ez)=1.0d-10
        enddo
      enddo

      call gop(phia,w2,'+  ',nz)
      call gop(w1,w2,'+  ',nz)
      call invcol2(phia,w1,nz)

      return
      end
C-----------------------------------------------------------------------
      subroutine x_planar_average(phia,phi,w1,w2)
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'ZPER'

      real phi(nx1,ny1,nz1,nelv),wrt(nx1,ny1,nz1,nelv),phia(nx1,nelx)
      real w1(nx1,nelx),w2(nx1,nelx) !work arrays

      integer e,eg,ex,nx,i,j,estride
      real xx,aa

      if(ldim.gt.2) then
        write(*,'(5x,a)')
     &            "x-average routine only written for 2D genbox meshes!"
        return
      endif

      nx=nx1*nelx !number of z-slices

      call rzero(phia,nx)
      call rzero(w1,nx)

      do e=1,nelt
        eg=lglel(e)
        ex=mod(eg,nelx) !x-slice id
        if(ex.eq.0)ex=nelx
        do i=1,nx1
          do j=1,ny1
            xx=(1.0-zgm1(i,1))/2.0
            aa=zz*area(j,1,4,e)+(1.0-zz)*area(j,1,2,e)
            w1(i,ex)=w1(i,ex)+aa
            phia(i,ex)=phia(i,ex)+aa*phi(i,j,1,e)
          enddo
        enddo
      enddo

      call gop(phia,w2,'+  ',nx)
      call gop(w1,w2,'+  ',nx)
      call invcol2(phia,w1,nx)

      return
      end
C-----------------------------------------------------------------------
      subroutine x_average_weighted(phia,phi,wrt,w1,w2)
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'ZPER'

      real phi(nx1,ny1,nz1,nelv),wrt(nx1,ny1,nz1,nelv),phia(nx1,nelx)
      real w1(nx1,nelx),w2(nx1,nelx) !work arrays

      integer e,eg,ex,nx,i,j,estride
      real xx,aa

      if(lz1.gt.1) then
        write(*,'(5x,a)')
     &                   "x-average routine only written for 2D meshes!"
        return
      endif

      nx=nx1*nelx !number of z-slices

      call rzero(phia,nx)
      call rzero(w1,nx)

      do e=1,nelt
        eg=lglel(e)
        ex=mod(eg,nelx) !x-slice id
        if(ex.eq.0)ex=nelx
        do i=1,nx1
          do j=1,ny1
            xx=(1.0-zgm1(i,1))/2.0
            aa=zz*area(j,1,4,e)+(1.0-zz)*area(j,1,2,e)
            w1(i,ex)=w1(i,ex)+aa*wrt(i,j,1,e)
            phia(i,ex)=phia(i,ex)+aa*wrt(i,j,1,e)*phi(i,j,1,e)
          enddo
          if(abs(w1(i,ex)).lt.1.0d-10) w1(i,ex)=1.0d-10
        enddo
      enddo

      call gop(phia,w2,'+  ',nx)
      call gop(w1,w2,'+  ',nx)
      call invcol2(phia,w1,nx)

      return
      end
C-----------------------------------------------------------------------
      real function q_vol_periodic(ix,iy,iz,ie,ifld)
      implicit none
      include 'SIZE'
      include 'TOTAL'
 
      integer ix,iy,iz,ie,ifld,n,e,f,dir

      logical ifdid(ldimt),ifprintNu(ldimt)
      common /printNu/ ifprintNu

      real dummy,sarea,tarea,time0(ldimt),tcorr
      real f_gm(ldimt),vel_avg,glsum,glsc2,glsc3

      data ifdid /ldimt*.false./
      data time0 /ldimt*-1.0/

      save ifdid,time0,vel_avg,f_gm,dir

      n=nx1*ny1*nz1*nelv

      if(.not.ifdid(ifld-1)) then
        dir=nint(abs(param(54))) !make sure this is an int, for my own sanity
        ifdid(ifld-1) = .true.
        tarea = 0.0
        do e=1,nelv
          do f=1,2*ndim
            if(cbc(f,e,ifld).eq.'f  ') then
              call surface_int(dummy,sarea,xm1,e,f)
              tarea=tarea+sarea
            endif
          enddo
        enddo
        tarea=glsum(tarea,1)
        f_gm(ifld-1)=abs(tarea/volvm1) !probably wrong for CHT
      endif

      ifprintNu(ifld-1)=.true.

      e=lglel(ie) !do nothing in solid region
      if(e.gt.nelgv) then
        q_vol_periodic = 0.0
        return
      endif

      if(time.ne.time0(ifld-1)) then
        time0(ifld-1)=time
        if(dir.eq.1) then
          vel_avg=glsc2(vx,bm1,n)/volvm1
          tcorr = -1.0*glsc3(t(1,1,1,1,ifld-1),vx,bm1,n)
        elseif(dir.eq.2) then
          vel_avg=glsc2(vy,bm1,n)/volvm1
          tcorr = -1.0*glsc3(t(1,1,1,1,ifld-1),vy,bm1,n)
        elseif(dir.eq.3) then
          vel_avg=glsc2(vz,bm1,n)/volvm1
          tcorr = -1.0*glsc3(t(1,1,1,1,ifld-1),vz,bm1,n)
        endif
        tcorr=tcorr/(vel_avg*volvm1)
        call cadd (t(1,1,1,1,ifld-1),tcorr,n)
      endif

      if(dir.eq.1) q_vol_periodic=-f_gm(ifld-1)*vx(ix,iy,iz,ie)/vel_avg
      if(dir.eq.2) q_vol_periodic=-f_gm(ifld-1)*vy(ix,iy,iz,ie)/vel_avg
      if(dir.eq.3) q_vol_periodic=-f_gm(ifld-1)*vz(ix,iy,iz,ie)/vel_avg

      return
      end
c-----------------------------------------------------------------------
      subroutine print_Nusselt
      include 'SIZE'
      include 'TOTAL'

      logical ifdo,ifprintNu(ldimt)
      real rNus(ldimt)

      common /printNu/ ifprintNu

      save ifdo
      data ifdo /.true./

      if(.not.ifdo) return
      if(istep.eq.0) then
        do i=1,ldimt
          ifprintNu(i)=.false.
        enddo 
        ifldsv=ifield
        do ifield=2,nfield
          call makeuq
        enddo
        ifield=ifldsv
      endif
c     if(nsteps.eq.0) then
c       do i=1,ldimt
c         ifprintNu(i)=.true.
c       enddo
c     endif

      if(abs(param(54)).lt.0.1) then
        ifdo=.false.
        if(nio.eq.0) then
          write(*,'(a)') "******************************"
          write(*,'(a)')
     &       "print_Nusselt routine only compatible with forced flow"
          write(*,'(a)') "******************************"
        endif
        return
      endif

      n=nx1*ny1*nz1*nelv

      jfld=0
      do ifld=2,ldimt1
        if(ifprintNu(ifld-1)) then
          jfld=jfld+1
          tarea=0.0
          twall=0.0
          do ie=1,nelt
            do ifa=1,2*ndim
              if(cbc(ifa,ie,ifld).eq.'f  ') then
                call surface_int(swall,sarea,t(1,1,1,1,ifld-1),ie,ifa)
                twall=twall+swall
                tarea=tarea+sarea
              endif
            enddo
          enddo
          tarea=glsum(tarea,1)
          twall=glsum(twall,1)/tarea

          if(nint(abs(param(54))).eq.1) 
     &           tbulk=glsc3(t(1,1,1,1,ifld-1),vx,bm1,n)/glsc2(vx,bm1,n)
          if(nint(abs(param(54))).eq.2) 
     &           tbulk=glsc3(t(1,1,1,1,ifld-1),vy,bm1,n)/glsc2(vy,bm1,n)
          if(nint(abs(param(54))).eq.3) 
     &           tbulk=glsc3(t(1,1,1,1,ifld-1),vz,bm1,n)/glsc2(vz,bm1,n)
          rNus(jfld)=1.0/((twall-tbulk)*cpfld(ifld,1))
        endif
      enddo

      if(nio.eq.0) then
        write(*,'(a24,es15.7,10es15.7)') 
     &                           "time, Nusselt",time,(rNus(i),i=1,jfld)
      endif

      return
      end
c-----------------------------------------------------------------------
      real function bc_area_average(phi,bca,ifld)
c     return the area-weighted average of scalar phi on BC bca
      implicit none
      include 'SIZE'
      include 'INPUT'

      character*3 bca
      integer ifld
      real phi(lx1*ly1*lz1*lelv)

      integer f,e
      real phibc,Abc,dphi,dA
      real glsum

      phibc=0.0
      Abc=0.0

      do 10 e=1,nelt
      do 10 f=1,ndim*2
        if(cbc(f,e,ifld).eq.bca) then
          call surface_int(dphi,dA,phi,e,f)
          phibc=phibc+dphi
          Abc=Abc+dA
        endif
 10   continue
      Abc=glsum(Abc,1)
      phibc=glsum(phibc,1)/Abc

      bc_area_average = phibc

      return
      end
c-----------------------------------------------------------------------
      real function bcID_area_average(phi,iID)
c     return the area-weighted average of scalar phi on boundary iID
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      integer iID
      real phi(lx1*ly1*lz1*lelv)

      integer f,e
      real phibc,Abc,dphi,dA
      real glsum

      phibc=0.0
      Abc=0.0

      do 10 e=1,nelt
      do 10 f=1,ndim*2
        if(boundaryID(f,e).eq.iID.or.boundaryIDt(f,e).eq.iID) then
          call surface_int(dphi,dA,phi,e,f)
          phibc=phibc+dphi
          Abc=Abc+dA
        endif
 10   continue
      Abc=glsum(Abc,1)
      phibc=glsum(phibc,1)/Abc

      bcID_area_average = phibc

      return
      end
c-----------------------------------------------------------------------
      real function bcID_area_integral(phi,iID)
c     return the area-integral of scalar phi on boundary iID
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      integer iID
      real phi(lx1*ly1*lz1*lelv)

      integer f,e
      real phibc,dphi,dA
      real glsum

      phibc=0.0

      do 10 e=1,nelt
      do 10 f=1,ndim*2
        if(boundaryID(f,e).eq.iID.or.boundaryIDt(f,e).eq.iID) then
          call surface_int(dphi,dA,phi,e,f)
          phibc=phibc+dphi
        endif
 10   continue
      phibc=glsum(phibc,1)

      bcID_area_integral = phibc

      return
      end
c-----------------------------------------------------------------------
      real function bc_area(bca,ifld)
c     return the total area of faces with BC bca
      implicit none
      include 'SIZE'
      include 'INPUT'

      character*3 bca
      integer ifld

      integer f,e
      real Abc,dA
      real glsum

      Abc=0.0

      do 10 e=1,nelt
      do 10 f=1,ndim*2
        if(cbc(f,e,ifld).eq.bca) then
          call surface_area(dA,e,f)
          Abc=Abc+dA
        endif
 10   continue
      Abc=glsum(Abc,1)

      bc_area = Abc

      return
      end
c-----------------------------------------------------------------------
      real function bcID_area(iID)
c     return the total area of faces with Boundary ID iID
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      integer iID

      integer f,e
      real Abc,dA
      real glsum

      Abc=0.0

      do 10 e=1,nelt
      do 10 f=1,ndim*2
        if(boundaryID(f,e).eq.iID.or.boundaryIDt(f,e).eq.iID) then
          call surface_area(dA,e,f)
          Abc=Abc+dA
        endif
 10   continue
      Abc=glsum(Abc,1)

      bcID_area = Abc

      return
      end
c-----------------------------------------------------------------------
      real function bc_flux_average(phi,bca,ifld)
c     return the flux-weighted average of phi on BC bca
      implicit none
      include 'SIZE'
      include 'TOTAL'

      character*3 bca
      integer ifld
      real phi(lx1*ly1*lz1,lelv)

      integer f,e,lxyz
      parameter (lxyz=lx1*ly1*lz1)
      real phibc,AA,dphi,dAA,w(lxyz)
      real glsum

      phibc=0.0
      AA=0.0

      do e=1,nelv
        do f=1,ndim*2
          if(cbc(f,e,ifld).eq.bca) then
            call surface_flux2(dphi,phi,e,f)
            call surface_flux(dAA,vx,vy,vz,e,f,w)
            phibc=phibc+dphi
            AA=AA+dAA
          endif
        enddo
      enddo
      AA=glsum(AA,1)
      phibc=glsum(phibc,1)/AA

      bc_flux_average = phibc

      return
      end
c-----------------------------------------------------------------------
      real function bcID_flux_average(phi,bid)
c     return the flux-weighted average of phi on faces with boundary bid
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer bid
      real phi(lx1*ly1*lz1,lelv)

      integer f,e,lxyz
      parameter (lxyz=lx1*ly1*lz1)
      real phibc,AA,dphi,dAA,w(lxyz)
      real glsum

      phibc=0.0
      AA=0.0

      do e=1,nelv
      do f=1,ndim*2
        if(BoundaryID(f,e).eq.bid) then
          call surface_flux2(dphi,phi,e,f)
          call surface_flux(dAA,vx,vy,vz,e,f,w)
          phibc=phibc+dphi
          AA=AA+dAA
        endif
      enddo
      enddo
      AA=glsum(AA,1)
      phibc=glsum(phibc,1)/AA

      bcID_flux_average = phibc

      return
      end
c-----------------------------------------------------------------------
      real function bc_totalflux(phi,ptrans,pdiff,bca,ifld)
c     return the total flux through the faces with BC bca
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real phi(lx1*ly1*lz1,lelt)
      real ptrans(lx1*ly1*lz1,lelt)
      real pdiff(lx1*ly1*lz1,lelt)
      integer ifld
      character*3 bca

      real face_advective_flux,face_diffusive_flux,glsum

      real phibc
      integer f,e

      phibc=0.0
      do e=1,nelv
      do f=1,ldim*2
        if(cbc(f,e,ifld).eq.bca) then
          phibc=phibc+face_advective_flux(phi,ptrans,e,f)
          phibc=phibc+face_diffusive_flux(phi,pdiff,e,f)
        endif
      enddo
      enddo
      phibc=glsum(phibc,1)

      bc_totalflux = phibc

      return
      end
c-----------------------------------------------------------------------
      real function bcID_totalflux(phi,ptrans,pdiff,bcid)
c     return the total flux through the faces with boundaryID bcid
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real phi(lx1*ly1*lz1,lelt)
      real ptrans(lx1*ly1*lz1,lelt)
      real pdiff(lx1*ly1*lz1,lelt)
      integer bcid

      real face_advective_flux,face_diffusive_flux,glsum

      real phibc
      integer f,e

      phibc=0.0
      do e=1,nelv
      do f=1,ndim*2
        if(BoundaryID(f,e).eq.bcid) then
          phibc=phibc+face_advective_flux(phi,ptrans,e,f)
          phibc=phibc+face_diffusive_flux(phi,pdiff,e,f)
        endif
      enddo
      enddo
      phibc=glsum(phibc,1)

      bcId_totalflux = phibc

      return
      end
c-----------------------------------------------------------------------
      real function face_advective_flux(phi,ptrans,e,f)
c     return the advective flux of scalar phi through face f of element e
      implicit none
      include 'SIZE'
      include 'TOTAL' !probably just need SOLN and GEOM

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real phi(lxyz,lelt)
      real ptrans(lxyz,lelt)
      integer e,f

      real wkx(lxyz),wky(lxyz),wkz(lxyz),dphi,wkk(lxyz)
      integer iface,js1,jf1,jskip1,js2,jf2,jskip2

      call col3(wkx,vx(1,1,1,e),phi(1,e),lxyz)
      call col3(wky,vy(1,1,1,e),phi(1,e),lxyz)
      if(if3d) call col3(wkx,vz(1,1,1,e),phi(1,e),lxyz)

      call col2(wkx,ptrans(1,e),lxyz)
      call col2(wky,ptrans(1,e),lxyz)
      if(if3d) call col2(wkz,ptrans(1,e),lxyz)

      call surface_flux_local(dphi,wkx,wky,wkz,e,f,wkk)
      face_advective_flux=dphi

      return
      end
c-----------------------------------------------------------------------
      real function face_diffusive_flux(phi,pdiff,e,f)
c     return the diffusive flux of scalar phi through face f of element e
      implicit none
      include 'SIZE'
      include 'TOTAL' !probably just need SOLN and GEOM

      integer lxyz
      parameter (lxyz=lx1*ly1*lz1)

      real phi(lxyz,lelt)
      real pdiff(lxyz,lelt)
      integer e,f

      real wkx(lxyz),wky(lxyz),wkz(lxyz),dphi,wkk(lxyz)

      call gradm11(wkx,wky,wkz,phi,e)

      call col2(wkx,pdiff(1,e),lxyz)
      call col2(wky,pdiff(1,e),lxyz)
      if(if3d) call col2(wkz,pdiff(1,e),lxyz)

      call chsign(wkx,lxyz)
      call chsign(wky,lxyz)
      if(if3d) call chsign(wkz,lxyz)

      call surface_flux_local(dphi,wkx,wky,wkz,e,f,wkk)
      face_diffusive_flux=dphi

      return
      end
c-----------------------------------------------------------------------
      subroutine surface_flux_local(dq,qx,qy,qz,e,f,w)
c     identical to surface_flux in navier5, but qx,qy,qz are only a single element
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      parameter (l=lx1*ly1*lz1)

      real qx(l),qy(l),qz(l),w(lx1,ly1,lz1)
      integer e,f

      call           faccl3  (w,qx,unx(1,1,f,e),f)
      call           faddcl3 (w,qy,uny(1,1,f,e),f)
      if (if3d) call faddcl3 (w,qz,unz(1,1,f,e),f)

      call dsset(lx1,ly1,lz1)
      iface  = eface1(f)
      js1    = skpdat(1,iface)
      jf1    = skpdat(2,iface)
      jskip1 = skpdat(3,iface)
      js2    = skpdat(4,iface)
      jf2    = skpdat(5,iface)
      jskip2 = skpdat(6,iface)

      dq = 0
      i  = 0
      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         i = i+1
         dq    = dq + area(i,1,f,e)*w(j1,j2,1)
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine average_files(inbase,navg)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'AVG'

      character*(*) inbase
      character*8   ftail
      integer lbase,navg,n,n2,i,j
      logical ifxyo_s

      n=nx1*ny1*nz1*nelv
      n2=nx2*ny2*nz2*nelv

      if(navg.eq.0.or.nsteps.gt.0) return

      atime=0.0
      call rzero(uavg,n)
      call rzero(vavg,n)
      call rzero(wavg,n)
      call rzero(pavg,n2)
      do j=1,ldimt
        call rzero(tavg(1,1,1,1,j),n)
      enddo
      do i=1,navg
        if(i.lt.10) then
          write(ftail,'(a7,i1)')'0.f0000',i
        elseif(i.lt.100) then
          write(ftail,'(a6,i2)')'0.f000',i
        elseif(i.lt.1000) then
          write(ftail,'(a5,i3)')'0.f00',i
        endif
        call blank(initc(1),132)
        initc(1)=trim(inbase)//ftail

        call restart(1)

        atime=atime+time
        call add2s2(uavg,vx,time,n)
        call add2s2(vavg,vy,time,n)
        call add2s2(wavg,vz,time,n)
        call add2s2(pavg,pr,time,n2)
        do j=1,ldimt
          call add2s2(tavg(1,1,1,1,j),t(1,1,1,1,j),time,n)
        enddo
      enddo
      time=atime
      call cmult(uavg,1.0/atime,n)
      call cmult(vavg,1.0/atime,n)
      call cmult(wavg,1.0/atime,n)
      call cmult(pavg,1.0/atime,n2)
      do j=1,ldimt
        call cmult(tavg(1,1,1,1,j),1.0/atime,n)
      enddo

      call copy (vx,uavg,n)
      call copy (vy,vavg,n)
      call copy (vz,wavg,n)
      call copy (pr,pavg,n2)
      do j=1,ldimt
        call copy(t(1,1,1,1,j),tavg(1,1,1,1,j),n)
      enddo

      if(nio.eq.0) write(*,*) "  average data:"
      call print_limits !print out the average data

      ifxyo_s = ifxyo
      ifxyo=.true.

      call prepost (.true.,'AVG')

      ifxyo = ifxyo_s

      return
      end
c-----------------------------------------------------------------------
      subroutine get_face_m1centroid(xx,yy,zz,rr,ie,iface)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ie,iface
      real xx,yy,zz,r0,rr

      integer i,i0,i1,j,j0,j1,k,k0,k1,ifc
      real bm0

      bm0=0.0
      xx=0.0
      yy=0.0
      zz=0.0
      rr=0.0
      ifc=0
      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
      do 10 k=k0,k1
      do 10 j=j0,j1
      do 10 i=i0,i1
        ifc=ifc+1
        bm0=bm0+area(ifc,1,iface,ie)
        xx=xx+xm1(i,j,k,ie)*area(ifc,1,iface,ie)
        yy=yy+ym1(i,j,k,ie)*area(ifc,1,iface,ie)
        zz=zz+zm1(i,j,k,ie)*area(ifc,1,iface,ie)
        r0=sqrt(xm1(i,j,k,ie)**2+ym1(i,j,k,ie)**2)
        rr=rr+r0*area(ifc,1,iface,ie)
 10   continue

      xx=xx/bm0
      yy=yy/bm0
      zz=zz/bm0
      rr=rr/bm0

      return
      end
c-----------------------------------------------------------------------
      subroutine get_elem_m1centroid(xx,yy,zz,ie)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'MASS'

      integer ie
      real xx,yy,zz

      integer ipt
      real bb

      xx=0.0
      yy=0.0
      zz=0.0
      bb=0.0
      do ipt = 1,lx1*ly1*lz1
        xx=xx+xm1(ipt,1,1,ie)*bm1(ipt,1,1,ie)
        yy=yy+ym1(ipt,1,1,ie)*bm1(ipt,1,1,ie)
        zz=zz+zm1(ipt,1,1,ie)*bm1(ipt,1,1,ie)
        bb=bb+bm1(ipt,1,1,ie)
      enddo
      xx=xx/bb
      yy=yy/bb
      zz=zz/bb

      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_point_2d(x1,y1,x0,y0,theta,xo,yo)
      implicit none

C     rotate point x1,y1 about point x0,y0 along the z-axis

      real x1,y1,x0,y0,theta,xo,yo

      xo=(x1-x0)*cos(theta)-(y1-y0)*sin(theta)+x0
      yo=(x1-x0)*sin(theta)+(y1-y0)*cos(theta)+y0

      return
      end
c-----------------------------------------------------------------------
      subroutine flag_bndry(bcc,ifld,phi)
      implicit none
      include 'SIZE'
      include 'INPUT'

      character*3 bcc 
      integer ifld
      real phi(lx1,ly1,lz1,1)

      integer iel,ifc,i,i0,i1,j,j0,j1,k,k0,k1,n
 
      n=lx1*ly1*lz1*nelv
      call rzero(phi,n)

      do 10 iel=1,nelt
      do 10 ifc=1,ndim*2
        if(cbc(ifc,iel,ifld).eq.bcc) then
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            phi(i,j,k,iel)=1.0
 20       continue
        endif
 10   continue
 
      return
      end
c-----------------------------------------------------------------------
      subroutine flag_bid(bid,phi,const)
      implicit none
      include 'SIZE'
      include 'GEOM'

      integer bid,iglsum
      real phi(lx1,ly1,lz1,1),const

      integer iel,ifc,i,i0,i1,j,j0,j1,k,k0,k1,n,ict
 
      n=lx1*ly1*lz1*nelv
      call rzero(phi,n)
      ict=0

      do 10 iel=1,nelt
      do 10 ifc=1,ndim*2
        if(BoundaryID(ifc,iel).eq.bid) then
          ict=ict+1
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            phi(i,j,k,iel)=const
 20       continue
        endif
 10   continue
      n=iglsum(ict,1)
      if(nio.eq.0) write(*,*)bid,n
 
      return
      end
c-----------------------------------------------------------------------
      subroutine walltime(tfin)
      implicit none 
C
C     ends the run and dumps a restart file if the wall time is greater 
C     than tfin (in seconds)
C     CALL BEFORE AVG_ALL
C
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      real time0,timenow,tfin,glmax
      save time0
      data time0 /0.0/

      if(istep.eq.0) then
        time0=dnekclock()
        return
      elseif(istep.gt.0) then
        if(nid.eq.0) then
          timenow=dnekclock()-time0
          if(tfin.lt.timenow) then
            write(6,*)"at wall time limit, last time step..."
            lastep=1
          endif
        endif
        call bcast(lastep,isize)
        if(lastep.eq.1) ifoutfld=.true.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine surface_flux2(dq,q,e,f)

      include 'SIZE'
      include 'TOTAL'
      parameter (l=lx1*ly1*lz1)

      real q(l,1),w(lx1,ly1,lz1)
      integer e,f

      call           faccl4  (w,vx(1,1,1,e),q(1,e),unx(1,1,f,e),f)
      call           faddcl4 (w,vy(1,1,1,e),q(1,e),uny(1,1,f,e),f)
      if (if3d) call faddcl4 (w,vz(1,1,1,e),q(1,e),unz(1,1,f,e),f)

      call dsset(lx1,ly1,lz1)
      iface  = eface1(f)
      js1    = skpdat(1,iface)
      jf1    = skpdat(2,iface)
      jskip1 = skpdat(3,iface)
      js2    = skpdat(4,iface)
      jf2    = skpdat(5,iface)
      jskip2 = skpdat(6,iface)

      dq = 0
      i  = 0
      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         i = i+1
         dq    = dq + area(i,1,f,e)*w(j1,j2,1)
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine faccl4(a,b,c,d,iface1)
C
C     Collocate B with A on the surface IFACE1 of element IE.
C
C         A, B, and C are (NX,NY,NZ) data structures
C         D is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation
C         IFACE  is the dssum notation.
C
      include 'SIZE'
      include 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1,LZ1),D(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(lx1,ly1,lz1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = B(J1,J2,1)*C(J1,J2,1)*D(I,1)
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine faddcl4(a,b,c,d,iface1)
C
C     Collocate B with C and add to A on the surface IFACE1 of element
C     IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,NZ) data structure
C         C is a (NX,NY,NZ) data structure
C         D is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation
C         IFACE  is the dssum notation.
C         29 Jan 1990 18:00 PST   PFF
C
      include 'SIZE'
      include 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1,LZ1),D(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(lx1,ly1,lz1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = A(J1,J2,1) + B(J1,J2,1)*C(J1,J2,1)*D(I,1)
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine loadmesh(string)
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'RESTART'
    
      integer l,ltrunc
      real t_s
      character string*(*)
      character*132 initc_sav

      l=ltrunc(string,len(string))
      if(l.gt.130) call exitti('invalid string length$',l)

      initc_sav=initc(1)
      t_s=time
      call blank(initc(1),132)
      initc(1)=trim(string)//' X'
      call restart(1)
      initc(1)=initc_sav
      time=t_s
 
      return
      end
c-----------------------------------------------------------------------
      subroutine loadfld2(string)
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'RESTART'
    
      integer l,ltrunc
      real t_s
      character string*(*)
      character*132 initc_sav

      l=ltrunc(string,len(string))
      if(l.gt.132) call exitti('invalid string length$',l)

      initc_sav=initc(1)
      t_s=time
      call blank(initc(1),132)
      initc(1)=trim(string)
      call restart(1)
      initc(1)=initc_sav
      time=t_s
 
      return
      end
c-----------------------------------------------------------------------
      subroutine save_ioflags

      implicit none

      include 'SIZE'
      include 'INPUT'

      logical ifbackup(0:ldimt+4)
      integer i

      common /ioflagss/ ifbackup 

      ifbackup(0)=ifxyo
      ifbackup(1)=ifvo
      ifbackup(2)=ifpo
      ifbackup(3)=ifto
      do i = 1,ldimt1
        ifbackup(3+i)=ifpsco(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine restore_ioflags

      implicit none

      include 'SIZE'
      include 'INPUT'

      logical ifbackup(0:ldimt+4)
      integer i

      common /ioflagss/ ifbackup 

      ifxyo=ifbackup(0)
      ifvo =ifbackup(1)
      ifpo =ifbackup(2)
      ifto =ifbackup(3)
      do i = 1,ldimt1
        ifpsco(i)=ifbackup(3+i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine clear_ioflags

      implicit none

      include 'SIZE'
      include 'INPUT'

      integer i

      ifxyo=.false.
      ifvo =.false.
      ifpo =.false.
      ifto =.false.
      do i = 1,ldimt1
        ifpsco(i)=.false.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpmesh(na3in)

      implicit none

      include 'SIZE'
      include 'TOTAL'

      character*3 na3in,na3
      logical ifxyo_s,ifpo_s,ifvo_s,ifto_s,ifpsco_s(ldimt1)

      integer i

      na3='msh'
      if(na3in.ne.'   ') na3=na3in

      call save_ioflags
      call clear_ioflags

      ifxyo=.true.

      call prepost (.true.,na3)

      call restore_ioflags

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpBCIDs(na3in)

      implicit none

      include 'SIZE'
      include 'TOTAL'

      character*3 na3in,na3

      integer iel,ifc,i,n,i0,i1,k0,k1,j0,j1,j,k
      integer nBCIDs, iglsum

      n=lx1*ly1*lz1
      nBCIDs = 0

      na3='bid'
      if(na3in.ne.'   ') na3=na3in

      call save_ioflags
      call clear_ioflags

      ifxyo=.true.
      ifto=.true.

      call rzero(t,n*nelt)
c     call izero(out_mask,nelt)

      do iel=1,nelt
      do ifc=1,2*ldim
        if(BoundaryID(ifc,iel).ne.0) then
c         out_mask(iel) = 1
          nBCIDs=nBCIDs+1
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            t(i,j,k,iel,1)=BoundaryID(ifc,iel)
 20       continue
        endif
      enddo
      enddo

      nBCIDs = iglsum(nBCIDs,1)
      if(nBCIDs.eq.0) then
        if(nio.eq.0) write(*,*) 
     &            "Warning in dumpBCIDs: no non-zero boundary IDs found"
      else
        call prepost (.true.,na3)
      endif

      call restore_ioflags

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpscalars(nsca)

      implicit none

      include 'SIZE'
      include 'TOTAL'

      integer i,nsca

      call save_ioflags
      call clear_ioflags

      ifxyo=.true.

      if(nsca.ge.1) ifto=.true.
      do i=1,nsca-1
        ifpsco(i)=.true.
      enddo
      do i=nsca,ldimt1
        ifpsco(i)=.false.
      enddo

      call prepost (.true.,'sca')

      call restore_ioflags

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpcfl

      implicit none

      include 'SIZE'
      include 'TOTAL'

      real DVC,DV1,DV2,DFC
      COMMON /SCRNS/ DVC  (LX1,LY1,LZ1,LELV),
     $               DV1  (LX1,LY1,LZ1,LELV),
     $               DV2  (LX1,LY1,LZ1,LELV),
     $               DFC  (LX1,LY1,LZ1,LELV)

      integer i,n
      real tsv(lx1,ly1,lz1,lelv),cdum,psv(lx2,ly2,lz2,lelv)

      n = lx1*ly1*lz1*nelv

      call save_ioflags
      call clear_ioflags

      ifxyo=.true.
      ifto=.true.

      if(ifsplit) then
        ifpo=.true.
        call copy(psv,pr,n)

c       call qthermal  !something wrong here
c       call add2 (qtl,usrdiv,n)
        call rzero(qtl,n)

        CALL OPDIV   (DVC,VX,VY,VZ)
        CALL DSSUM   (DVC,lx1,ly1,lz1)
        CALL COL2    (DVC,BINVM1,n)

        CALL SUB3    (DFC,DVC,QTL,n)
        CALL COL3    (pr,DFC,DFC,n)
      endif

      call copy(tsv,t,n)
      call compute_cfl(cdum,vx,vy,vz,dt)
      call copy(t,cflf,n)
      call prepost (.true.,'cfl')
      call copy(t,tsv,n)
      if(ifsplit) call copy(pr,psv,n)

      call restore_ioflags

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpmetrics

      implicit none

      include 'SIZE'
      include 'TOTAL'

      real DVC,DV1,DV2,DFC
      COMMON /SCRNS/ DVC  (LX1,LY1,LZ1,LELV),
     $               DV1  (LX1,LY1,LZ1,LELV),
     $               DV2  (LX1,LY1,LZ1,LELV),
     $               DFC  (LX1,LY1,LZ1,LELV)

      integer i,n
      real tsv(lx1*ly1*lz1*lelv,ldimt),cdum,psv(lx2,ly2,lz2,lelv)
      real wd(lx1*ly1*lz1*lelv)

      n = lx1*ly1*lz1*nelv

      call save_ioflags
      call clear_ioflags

      ifxyo=.true.
      ifto=.true.
      ifpsco(1)=.true.

      if(ifsplit) then
        ifpo=.true.
        call copy(psv,pr,n)

c       call qthermal  !something wrong here
c       call add2 (qtl,usrdiv,n)
        call rzero(qtl,n)

        CALL OPDIV   (DVC,VX,VY,VZ)
        CALL DSSUM   (DVC,lx1,ly1,lz1)
        CALL COL2    (DVC,BINVM1,n)

        CALL SUB3    (DFC,DVC,QTL,n)
        CALL COL3    (pr,DFC,DFC,n)
      endif

      call copy(tsv(1,1),t,n)
      call compute_cfl(cdum,vx,vy,vz,dt)
      call copy(t,cflf,n)

      call copy(tsv(1,2),t(1,1,1,1,2),n)
      call get_wall_distance(wd,2)
      call get_y_p(wd,t(1,1,1,1,2),.true.)
      
      call prepost (.true.,'mtr')
      call copy(t,tsv,n)
      call copy(t(1,1,1,1,2),tsv(1,2),n)
      if(ifsplit) call copy(pr,psv,n)

      call restore_ioflags

      return
      end
c-----------------------------------------------------------------------
      subroutine surface_area(sarea,e,f)

      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'TOPOL'

      integer e,f

      call dsset(lx1,ly1,lz1)

      iface  = eface1(f)
      js1    = skpdat(1,iface)
      jf1    = skpdat(2,iface)
      jskip1 = skpdat(3,iface)
      js2    = skpdat(4,iface)
      jf2    = skpdat(5,iface)
      jskip2 = skpdat(6,iface)

      sarea = 0.
      i     = 0

      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         i = i+1
         sarea = sarea+area(i,1,f,e)
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine store_solution
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      real w1,w2,w3,w4,w5,tsv
      common /STORE/ w1(lx1*ly1*lz1*lelv)
     &              ,w2(lx1*ly1*lz1*lelv)
     &              ,w3(lx1*ly1*lz1*lelv)
     &              ,w4(lx2*ly2*lz2*lelv)
     &              ,w5(lx1*ly1*lz1*lelv,ldimt)
     &              ,tsv

      n1=lx1*ly1*lz1*nelv
      n2=lx2*ly2*lz2*nelv
      nt=lx1*ly1*lz1*nelt

      call copy(w1,vx,n1)
      call copy(w2,vy,n1)
      if(if3d) call copy(w3,vz,n1)
      call copy(w4,pr,n2)
      do i=1,ldimt
        call copy(w5(1,i),t(1,1,1,1,i),nt)
      enddo
      tsv=time

      return
      end
c-----------------------------------------------------------------------
      subroutine reload_solution
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      real w1,w2,w3,w4,w5,tsv
      common /STORE/ w1(lx1*ly1*lz1*lelv)
     &              ,w2(lx1*ly1*lz1*lelv)
     &              ,w3(lx1*ly1*lz1*lelv)
     &              ,w4(lx2*ly2*lz2*lelv)
     &              ,w5(lx1*ly1*lz1*lelv,ldimt)
     &              ,tsv

      n1=lx1*ly1*lz1*nelv
      n2=lx2*ly2*lz2*nelv
      nt=lx1*ly1*lz1*nelt

      call copy(vx,w1,n1)
      call copy(vy,w2,n1)
      if(if3d) call copy(vz,w3,n1)
      call copy(pr,w4,n2)
      do i=1,ldimt
        call copy(t(1,1,1,1,i),w5(1,i),nt)
      enddo
      time=tsv

      return
      end
c-----------------------------------------------------------------------
      subroutine convertvel_rtz
      include 'SIZE'
      include 'TOTAL'

      n=lx1*ly1*lz1*nelv

      if(nio.eq.0) write(6,'(a)') 
     &                    "   converting velocity to radial coordinates"

      do i=1,n
        xx=xm1(i,1,1,1)
        yy=ym1(i,1,1,1)
        tht=atan2(yy,xx)
        ux=vx(i,1,1,1)
        uy=vy(i,1,1,1)
        ur=ux*cos(tht)+uy*sin(tht)
        ut=-ux*sin(tht)+uy*cos(tht)
        vx(i,1,1,1)=ur
        vy(i,1,1,1)=ut
      enddo

      return
      end
C-----------------------------------------------------------------------
      subroutine cfill_face(ifc,iel,phi,cc)
      implicit none
      include 'SIZE'

      integer ifc,iel
      real cc, phi(lx1,ly1,lz1,*)

      integer i0,i1,j0,j1,k0,k1,i,j,k

      call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
      do 100 k=k0,k1
      do 100 j=j0,j1
      do 100 i=i0,i1
        phi(i,j,k,iel) = cc
 100  continue

      return
      end
C-----------------------------------------------------------------------
      subroutine cfill_rface(ifc,iel,phi,cc)
      implicit none
      include 'SIZE'

      integer ifc,iel
      real cc, phi(lx1,ly1,lz1,*)

      integer i0,i1,j0,j1,k0,k1,i,j,k

      call facindr(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
      do 100 k=k0,k1
      do 100 j=j0,j1
      do 100 i=i0,i1
        phi(i,j,k,iel) = cc
 100  continue

      return
      end
C-----------------------------------------------------------------------
      subroutine cadd_face(ifc,iel,phi,cc)
      implicit none
      include 'SIZE'

      integer ifc,iel
      real cc, phi(lx1,ly1,lz1,*),ph0

      integer i0,i1,j0,j1,k0,k1,i,j,k

      call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
      do 100 k=k0,k1
      do 100 j=j0,j1
      do 100 i=i0,i1
        ph0 = phi(i,j,k,iel)
        phi(i,j,k,iel) = ph0+cc
 100  continue

      return
      end
C-----------------------------------------------------------------------
