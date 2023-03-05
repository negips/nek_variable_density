c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)

      implicit none  
  
      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TSTEP'
      include 'NEKUSE'

      integer e,ix,iy,iz,ieg

      real rhog         ! density of air
      real rhow         ! density of water
      real mug          ! dynamic viscosity of air
      real muw          ! dynamic viscosity of water
      real nug          ! kinematic viscosity of air
      real nuw          ! kinematic viscosity of water
      real alpha

      real setvp

      real distn, eps

!     densities      
      rhog   = uparam(1)         ! 1.2061e-3
      rhow   = uparam(2)         ! 1.0

!     viscosities      
      mug    = uparam(3)         ! 1.5052e-6
      muw    = uparam(4)         ! 8.3e-5

      nug    = mug/rhog
      nuw    = muw/rhow

      eps    = 1.0e-1

      if (ifield.eq.1) then
        utrans = setvp(rhog,rhow,-temp,eps)
        udiff  = setvp(mug,muw,-temp,eps)
      endif  

      if (ifield .eq. 2) then
        e = gllel(ieg)
        utrans = 1.0   
        udiff = param(8)
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)

      implicit none        
  
      include 'SIZE'
      include 'INPUT'
      include 'NEKUSE'
      include 'TSTEP'

      integer ix,iy,iz,ieg
      real x0,y0,mu
      real xo,yo        ! old
      real cs,sn,th
      real amp
      real omega

!     Location based on old coordinates 
      x0 = 1.0 
      y0 = 0.0
      mu = 0.1

!      pi       = 4.0*atan(1.0)
      th       = uparam(7)*pi/180.0
      cs       = cos(th)
      sn       = sin(th)

      xo = cs*x + sn*y        ! Old x
      yo = -sn*x + cs*y       ! Old y
      amp = uparam(6)
      omega = uparam(8)

      ffx = 0.0
      ffy = amp*exp(-((xo-x0)*(yo-y0)/mu/mu)**2)*sin(2*pi*omega*time)
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'

      integer ix,iy,iz,ieg

      qvol =  0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'MVGEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'PARALLEL'      ! nelgv

      include 'TEST'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      integer i,j,k,e,n,n2

      real x,y,z

      real gl2norm,glsc2,glsc3

      real rhotr2(lx2*ly2*lz2*lelv)
      real exrhotr(lx2*ly2*lz2*lelv)
      real rhotrlag(lx2*ly2*lz2*lelv,2)
      common /density_var/ rhotr2,exrhotr,rhotrlag


      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      if (istep.eq.0) then

        if (param(95).gt.0) then
          param(95) = 50        ! start of projections
        endif

        call rone(vtrans(1,1,1,1,2),n)
        call rone(vdiff(1,1,1,1,2),n)

        call phi0(t)

        ifield = 1
        call vprops
        ifield = 2
        call vprops

!        ifheat = .false.

        call frame_start

        ifto = .true.

        call outpost(v1mask,v2mask,v3mask,pr,tmask,'msk')

        call outpost(vx,vy,vz,pr,t,'   ')

12      format(A4,2x,16(E12.5,2x))


!!       Just testing
!        call rone(vx,n)
!        call rone(vy,n)
!
!        call calc_rhs_pressure()
!
!        call outpost(vx,vy,vz,rhotr2,t,'tst')
!        call exitt

        call gen_mapping_mvb()
!        call fs_get_intpos(fs_intpos)

        call exitt
      endif 

      call frame_monitor

      call chkpt_main

      ifto = .true.
!     Reset preconditioner
      if (mod(istep,iostep).eq.0) then
!       Reynolds number of the field
        do i=1,n
          t(i,1,1,1,2) = vtrans(i,1,1,1,1)*1.0*1.0/vdiff(i,1,1,1,1)
        enddo       
        call outpost2(vtrans,vdiff,t(1,1,1,1,2),exrhotr,
     $                vdiff(1,1,1,1,2),1,'vis')
      endif

      if (istep.eq.nsteps.or.lastep.eq.1) then
        call frame_end
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'
      include 'GEOM'
      include 'INPUT'

      integer ix,iy,iz,iside,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real rmid,xmid

      real mu,x0
      real cs,sn,th

      real xo,yo

      pi = 4.0*atan(1.0)
      th       = uparam(7)*pi/180.0
      cs       = cos(th)
      sn       = sin(th)

      ux   = 1.0*cs - 0.0*sn
      uy   = 1.0*sn + 0.0*cs
      uz   = 0.0

!     Temperature
      xo = cs*x + sn*y        ! Old x
      yo = -sn*x + cs*y       ! Old y
      temp = yo ! - y0        ! signed distance from interface


      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'PARALLEL'
      include 'NEKUSE'
      include 'GEOM'

      integer ix,iy,iz,ieg
      real pi

      real mu,x0
      real cs,sn,th

      real rnd

      call random_number(rnd)

      pi = 4.0*atan(1.0)

      th     = uparam(7)*pi/180.0
      cs     = cos(th)
      sn     = sin(th)

      ux   = 1.0*cs - 0.0*sn + 0.01*rnd
      uy   = 1.0*sn + 0.0*cs + 0.01*rnd
      uz   = 0.0


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices

      implicit none  
  
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'INPUT'
      include 'TSTEP'
      include 'PARALLEL'

      integer e,i,nc
      real th
      real cs,sn
      real x,y

      pi       = 4.0*atan(1.0)
      th       = uparam(7)*pi/180.0
      cs       = cos(th)
      sn       = sin(th)

      nc = 2**ndim


      do e=1,nelv
      do i=1,nc
        x       = xc(i,e)
        y       = yc(i,e)
        xc(i,e) = cs*x - sn*y
        yc(i,e) = sn*x + cs*y
      enddo
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates

      implicit none

      include 'SIZE'
      include 'INPUT'      ! cbc
      include 'PARALLEL'
      include 'GEOM'

      integer iel,ifc
      integer n


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3

      implicit none        

      include 'SIZE'
      include 'SOLN'    ! tmult
      include 'INPUT'
      include 'GEOM'

      integer i,n


      return
      end
c-----------------------------------------------------------------------
      subroutine phi0(phi)

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)
     
      integer i,n
      real x,y,xo,yo
      real r,r0
      real x0,y0

      real pi
      real cs,sn,th

      n        = lx1*ly1*lz1*nelv
 
      pi       = 4.0*atan(1.0)
      th       = uparam(7)*pi/180.0
      cs       = cos(th)
      sn       = sin(th)
     
      y0 = 0.0
      do i=1,n
        x = xm1(i,1,1,1)
        y = ym1(i,1,1,1)
        xo = cs*x + sn*y
        yo = -sn*x + cs*y
        phi(i) = yo ! - y0        ! signed distance from interface
      enddo  

      return
      end subroutine phi0
!---------------------------------------------------------------------- 

      real function heavyside(phi,eps)

      real phi,eps
      real pi

      pi = 4.0*atan(1.0) 

      if (phi.lt.-eps) then
        heavyside = 0.0
      elseif (phi.gt.eps) then
        heavyside = 1.0
      else
        heavyside = 0.5*(1.0 + phi/eps + 1.0/pi*sin(pi*phi/eps))
!        heavyside = 1.0*(1.0 + phi/eps + 0.0/pi*sin(pi*phi/eps))
      endif  

      return
      end function heavyside
!---------------------------------------------------------------------- 

      real function setvp(vg,vl,phi,eps)

      implicit none

      real vg     ! gas property
      real vl     ! liquid property
      real phi    ! disgned distance function
      real eps    
      real heavyside    ! heavyside function
      real heavy

      heavy = heavyside(phi,eps)
      setvp = vg + (vl - vg)*heavy      

      return
      end function 
!---------------------------------------------------------------------- 

      subroutine heavydist(phi,eps,off)

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)
      real phioff
      real heavyside          ! function
      
      integer i,j,k,e,n
      real hd
      real eps
      real off

      n=lx1*ly1*lz1*nelv

      do i=1,n
        phioff = phi(i)+off
        hd = heavyside(phioff,eps)
        phi(i) = hd
      enddo  

      return
      end subroutine heavydist
!---------------------------------------------------------------------- 

      subroutine outmat_formatted(a,m,n,name6,ie)

      real a(m,n)
      character*6 name6

      character*28 str

      call blank(str,28)

      write(str,7) '(I5,2x,A6,2x,',n,'(E15.8E2,2x))'
    7 format(a13,I2,A13)

      write(6,*) 
      write(6,*) ie,' matrix: ',name6,m,n,str

      do i=1,m
         write(6,str) ie,name6,(a(i,j),j=1,n)
      enddo
      write(6,*) 

      return
      end
c-----------------------------------------------------------------------

      subroutine incomprn_test (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure currection req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c     Notes  1.  up is _not_ scaled by bd(1)/dt.  This should be done
c                external to incompr().
c
c            2.  up accounts _only_ for the perturbation pressure,
c                not the current pressure derived from extrapolation.
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      include 'TEST'
c
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      parameter(nset = 1 + lbelv/lelv)
      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
      common /orthbi/ nprv(2)
      logical ifprjp

      real rhotr2(lx2*ly2*lz2*lelv),exrhotr(lx2*ly2*lz2*lelv)
      real rhotrlag(lx2*ly2*lz2*lelv,2)
      common /density_var/ rhotr2,exrhotr,rhotrlag



      ifprjp=.false.    ! Project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      if (icalld.eq.0) tpres=0.0
      icalld = icalld+1
      npres  = icalld
      etime1 = dnekclock()

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      intype = 1

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call invers2 (h2inv,h2,ntot1)

!      call opdiv   (dp,ux,uy,uz)
      call opdiv_rho (dp,ux,uy,uz)        ! prabal

      bdti = -bd(1)/dt
      call cmult   (dp,bdti,ntot2)

      call add2(dp,exrhotr,ntot2)         ! added extrapolated density transport

      call add2col2(dp,bm2,usrdiv,ntot2)  ! User-defined divergence.

      call ortho   (dp)

      i = 1 + ifield/ifldmhd
      if (ifprjp)   call setrhsp  (dp,h1,h2,h2inv,pset(1,i),nprv(i))
                    scaledt = dt/bd(1)
                    scaledi = 1./scaledt
                    call cmult(dp,scaledt,ntot2)        ! scale for tol
                    call esolver  (dp,h1,h2,h2inv,intype)
                    call cmult(dp,scaledi,ntot2)
      if (ifprjp)   call gensolnp (dp,h1,h2,h2inv,pset(1,i),nprv(i))

      call add2(up,dp,ntot2)

      call opgradt  (w1 ,w2 ,w3 ,dp)
      call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      dtb  = dt/bd(1)
      call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )

      if (ifmhd)  call chkptol	! to avoid repetition

      tpres=tpres+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------


      subroutine h2inv (out1,inp1,mask,h2)

!     Compute OUT = (H2*B)-1 * INP   (explicit)
      
      implicit none
  
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'

      include 'OPCTR'

      real out1  (1)
      real inp1  (1)
      real mask  (1)
      real h2    (1)

      real temp
      common /scrcg/ temp(lx1,ly1,lz1,lelt)

      integer n

#ifdef TIMER
      if (isclld.eq.0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'h2inv'
      endif
#endif

      n = lx1*ly1*lz1*nelv

      call copy(out1,inp1,n)
      call col2(out1,mask,n)
      call dssum (out1,lx1,ly1,lz1)

!     temp = H2*BM1      
      call col3 (temp,bm1,h2,n)
      call dssum (temp,lx1,ly1,lz1)
!     temp = (H2*BM1)^-1
      call invcol1(temp,n)

!     out1 = [(H2*BM1)^-1]*out1
      call col2(out1,temp,n)

!     print out
      if (nio.eq.0) write(6,*) '(H2^-1):: Done'      

      return
      end
c-----------------------------------------------------------------------





c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
