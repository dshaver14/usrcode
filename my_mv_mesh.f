c----------------------------------------------------------------------
      subroutine my_mv_mesh(dxo,dyo,dzo,nstps)

c     This subroutine solves for and applies the overall mesh 
c     displacement when provided with the displacment vector on boundaries
c     with the 'mv ' BC in field 0.

      include 'SIZE'
      include 'TOTAL'
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
      logical ifxyos
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
c     relax = 1.0/real(nstps)
c     call copy(umeshx,dxo,n)
c     call copy(umeshy,dyo,n)
c     call copy(umeshz,dzo,n)

c     call cmult(umeshx,relax,n)
c     call cmult(umeshy,relax,n)
c     call cmult(umeshz,relax,n)

      utx_usr=glamax(dxo,n)
      uty_usr=glamax(dyo,n)  
      utz_usr=glamax(dzo,n)

      if (nid.eq.0) then
        write(6,*) "utx_usr: ",utx_usr
        write(6,*) "uty_usr: ",uty_usr
        write(6,*) "utz_usr: ",utz_usr
      endif

      time = 0.0    
      ifxyos = ifxyo
      ifxyo = .true.
      call prepost(.true.,'mvm')
      ifxyo = ifxyos
      nbl = 0
c     call dumpmesh('mvm')

      do istep = 1,nstps
c       factor for decaying boundary layer preservation
        fact = real(istep-1)/real(nstps-1)
        fact = 1.0-fact

c       for non-linear relaxation
        call sub3(umeshx,finx,xm1,n)
        call sub3(umeshy,finy,ym1,n)
        call sub3(umeshz,finz,zm1,n)

        relax = 1.0 - (real(istep)/real(nstps)-1.0)**2
c       relax = real(istep)/real(nstps)
        if(nio.eq.0) write(*,*) "relaxation factor: ",relax
        call cmult(umeshx,relax,n)
        call cmult(umeshy,relax,n)
        call cmult(umeshz,relax,n)

        if (icalld.eq.0) then
          icalld=1
          call rone(mask,n)
          do e=1,nelv
          do f=1,nface
            if(cbc(f,e,0).eq.'W  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'W1 ')call facev(mask,e,f,one ,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'v  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'mv ')call facev(mask,e,f,zero,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'O  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
            if(cbc(f,e,0).eq.'mvb')then
              call facev(mask,e,f,zero,nx1,ny1,nz1)
              nbl = 1
            endif
          enddo
          enddo
          call dsop(mask,'*  ',nx1,ny1,nz1)    ! dsop mask
          call opzero(delx,dely,delz)
          nbl = iglsum(nbl,1)
        endif
    
        napprx(1)=0
        nappry(1)=0
        napprz(1)=0
        nxz   = nx1*nz1
        nxyz  = nx1*ny1*nz1

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
          deltap1 = 2.0*delta/real(lx1)  ! Protected b.l. thickness
          deltap2 = 1.0*delta

c       magic distribution - it really does a better job of preseving BLs 
          do i=1,n
            arg1   = -(d(i)/deltap1)**2
            arg2   = -(d(i)/deltap2)**2
            h1(i)  = h1(i) + fact*(100.0*exp(arg1) + 5.0*exp(arg2))
            t(i,1,1,1,1) = h1(i)
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
     &    call laplaceh('mshx',delx,h1,h2,mask,vmult,1,tol,
     &    1000,apprx,napprx)
        if (uty_usr.gt.1e-10) 
     &    call laplaceh('mshy',dely,h1,h2,mask,vmult,1,tol,
     &    1000,appry,nappry)
        if (utz_usr.gt.1e-10)
     &    call laplaceh('mshz',delz,h1,h2,mask,vmult,1,tol,
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

        call print_limits

        time = istep
        ifxyos = ifxyo
        ifxyo = .true.
        call prepost(.true.,'mvm')
        ifxyo = ifxyos
c       call dumpmesh('mvm')

        call fix_geom
        call mesh_metrics

      enddo

      ifield = ifield_sv

      return
      end
c-----------------------------------------------------------------------
      subroutine laplaceh
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
