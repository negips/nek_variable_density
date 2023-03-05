!====================================================================== 
!     
!     Author: Prabal Negi
!     Description: Terms for variable density pressure solve      
!      
!====================================================================== 
      subroutine opdiv_rho(outfld,inpx,inpy,inpz)
     
      implicit none

      include 'SIZE'
      include 'SOLN'          ! vtrans
      include 'TSTEP'         ! ifield

      real outfld(1),inpx(1),inpy(1),inpz(1)
      real tmp1,tmp2,tmp3
      common /scrsf/ tmp1(lx1,ly1,lz1,lelt),
     $               tmp2(lx1,ly1,lz1,lelt), 
     $               tmp3(lx1,ly1,lz1,lelt) 


      call opcopy(tmp1,tmp2,tmp3,inpx,inpy,inpz)
      call opcol2(tmp1,tmp2,tmp3,
     $            vtrans(1,1,1,1,ifield),vtrans(1,1,1,1,ifield),
     $            vtrans(1,1,1,1,ifield))

      call opdiv(outfld,inpx,inpy,inpz)

      return
      end subroutine
!---------------------------------------------------------------------- 
     
      subroutine density_transport_m2(rhotr2)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'   ! ifield
      include 'MASS'
      include 'GEOM'
      include 'DXYZ'
      include 'IXYZ'
      include 'WZ'

      integer i,e,n1,n2,nxyz1,nxyz2
      integer ifld

      real ta1,ta2,ta3,ta4
      common /scruz/ ta1(lx1*ly1*lz1,lelt),
     $               ta2(lx1*ly1*lz1,lelt),    
     $               ta3(lx1*ly1*lz1,lelt), 
     $               ta4(lx1*ly1*lz1,lelt) 

      real rhotr2(lx2*ly2*lz2,nelv)
      
      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2
      n1    = nxyz1*nelv
      n2    = nxyz2*nelv

      ifld = 1

      call rzero(rhotr2,n2)
      do e=1,nelv
        call copy(ta1,vtrans(1,1,1,e,ifld),nxyz1)

!       X-direction 
!       d\rho/dx        
        call tensor3_op(ta3,ta1,lx1,ly1,lz1,
     $                  dxm12,iytm12,iztm12,lx2,ly2,lz2)
        call col2(ta3,rxm2(1,1,1,e),nxyz2)    

        call tensor3_op(ta2,ta1,lx1,ly1,lz1,
     $                  ixm12,dytm12,iztm12,lx2,ly2,lz2)
        call col2(ta2,sxm2(1,1,1,e),nxyz2)    
        call add2(ta3,ta2,nxyz2)    
            
        if (ndim.eq.3) then
          call tensor3_op(ta2,ta1,lx1,ly1,lz1,
     $                  ixm12,iytm12,dztm12,lx2,ly2,lz2)
          call col2(ta2,txm2(1,1,1,e),nxyz2)    
          call add2(ta3,ta2,nxyz2)
        endif 
         
        call copy(ta4,vx(1,1,1,e),nxyz1)
!       U: M1 -> M2        
        call tensor3_op(ta2,ta4,lx1,ly1,lz1,
     $                  ixm12,iytm12,iztm12,lx2,ly2,lz2)

!       U*d\rho/dx        
        call col2(ta3,ta2,nxyz2)

!       Wt*U*d\rho/dx        
        call col2(ta3,w3m2,nxyz2)

        call add2(rhotr2(1,e),ta3,nxyz2)

!       Y-direction 
!       d\rho/dy        
        call tensor3_op(ta3,ta1,lx1,ly1,lz1,
     $                  dxm12,iytm12,iztm12,lx2,ly2,lz2)
        call col2(ta3,rym2(1,1,1,e),nxyz2)    

        call tensor3_op(ta2,ta1,lx1,ly1,lz1,
     $                  ixm12,dytm12,iztm12,lx2,ly2,lz2)
        call col2(ta2,sym2(1,1,1,e),nxyz2)    
        call add2(ta3,ta2,nxyz2)    
            
        if (ndim.eq.3) then
          call tensor3_op(ta2,ta1,lx1,ly1,lz1,
     $                  ixm12,iytm12,dztm12,lx2,ly2,lz2)
          call col2(ta2,tym2(1,1,1,e),nxyz2)    
          call add2(ta3,ta2,nxyz2)
        endif 
         
        call copy(ta4,vy(1,1,1,e),nxyz1)
!       V: M1 -> M2        
        call tensor3_op(ta2,ta4,lx1,ly1,lz1,
     $                  ixm12,iytm12,iztm12,lx2,ly2,lz2)

!       V*d\rho/dy
        call col2(ta3,ta2,nxyz2)

!       Wt*V*d\rho/dy        
        call col2(ta3,w3m2,nxyz2)

        call add2(rhotr2(1,e),ta3,nxyz2)


        if (ndim.eq.3) then
!         Z-direction 
!         d\rho/dz
          call tensor3_op(ta3,ta1,lx1,ly1,lz1,
     $                    dxm12,iytm12,iztm12,lx2,ly2,lz2)
          call col2(ta3,rzm2(1,1,1,e),nxyz2)    

          call tensor3_op(ta2,ta1,lx1,ly1,lz1,
     $                    ixm12,dytm12,iztm12,lx2,ly2,lz2)
          call col2(ta2,szm2(1,1,1,e),nxyz2)    
          call add2(ta3,ta2,nxyz2)    
              
          call tensor3_op(ta2,ta1,lx1,ly1,lz1,
     $                  ixm12,iytm12,dztm12,lx2,ly2,lz2)
          call col2(ta2,tzm2(1,1,1,e),nxyz2)    
          call add2(ta3,ta2,nxyz2)
           
          call copy(ta4,vz(1,1,1,e),nxyz1)
!         V: M1 -> M2        
          call tensor3_op(ta2,ta4,lx1,ly1,lz1,
     $                    ixm12,iytm12,iztm12,lx2,ly2,lz2)

!         W*d\rho/dz
          call col2(ta3,ta2,nxyz2)

!         Wt*V*d\rho/dy        
          call col2(ta3,w3m2,nxyz2)

          call add2(rhotr2(1,e),ta3,nxyz2)

        endif  


      enddo


      return
      end subroutine
!---------------------------------------------------------------------- 

      subroutine extrapolate_f(exf,f,flag1,flag2,n)

C     Sum up contributions to kth order extrapolation scheme.

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'MASS'

      integer n
      real ab0,ab1,ab2

      real exf(n),f(n),flag1(n),flag2(n)

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)

!     exf = ab1*flag1 + ab2*flag2      
      call add3s2 (exf,flag1,flag2,ab1,ab2,n)
      call copy   (flag2,flag1,n)
      call copy   (flag1,f,n)

!     exf = ab0*f + ab1*flag1 + ab2*flag2
      call add2s2 (exf,f,ab0,n)

      return
      END
C
c-----------------------------------------------------------------------

      subroutine calc_rhs_pressure()

      implicit none

      include 'SIZE'
      include 'TSTEP'

      real rhotr2(lx2*ly2*lz2*lelv)
      real exrhotr(lx2*ly2*lz2*lelv)
      real rhotrlag(lx2*ly2*lz2*lelv,2)
      common /density_var/ rhotr2,exrhotr,rhotrlag

      integer n

      n = lx2*ly2*lz2*nelv

      call density_transport_m2(rhotr2)
      call extrapolate_f(exrhotr,rhotr2,rhotrlag(1,1),
     $                   rhotrlag(1,2),n)

      return
      end subroutine
!---------------------------------------------------------------------- 

