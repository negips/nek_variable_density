!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for parallel sorting
!
!----------------------------------------------------------------------
!====================================================================== 
!---------------------------------------------------------------------- 
      subroutine fs_global_sort_mp(fld1,fld2,gnum,n)

!     mp - Multiple Processor        
!     Get global numbering of entries
!     When the Single processor work memory is not large enough
!     for all entries.

!     fld1  - primary field (input)       ! scruz
!     fld2  - secondary field (input)     ! scruz
!     gnum  - Unique global numbering (output)
!     n     - Local number of unique elements        

      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'FS_ALE'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

      integer n               ! No. of entries.
      real fld1(n)            ! Input field 1
      real fld2(n)            ! Input field 2
      integer gnum(n)         ! Permutation indicies

!     Temporary arrays
      integer wki
      common /scrcg/ wki (lt*2)

!     Temporary arrays
      real wk1,wk2,wk3
      integer ind
      common /scrns/ wk1(lt*2),
     $               wk2(lt*2),
     $               wk3(lt*2),
     $               ind(lt*2) 

!     Temporary arrays
      real wk5
      integer owner0
      common /scrsf/ wk5(lt*2),
     $               owner0(lt*2)    

!     Temporary arrays
      integer owner1,owner2
      common /scrvh/ owner1 (lt*2)   ! Owner
     $ ,             owner2 (lt*2)   ! 

!     Temporary arrays
      integer ninseg
      logical ifseg
      common /screv/ ninseg (lt*2)
     $ ,             ifseg  (lt*2)

!     Temporary arrays
      real ta1,ta2
      common /scruz/ ta1    (lt*2)
     $ ,             ta2    (lt*2)

!     Temporary arrays
      integer gnum1,gnum2
      common /gnumscrch/ gnum1    (lt*2)
     $ ,                 gnum2    (lt*2)

      integer nkey
      integer own
      integer iseg,nseg
      integer fdim
      integer nglobalu


      integer alen            ! Work array length
      integer alen2           ! Half of work array
      integer iglsum
      integer igl_running_sum
      integer iglmax

      integer gl_n            ! Total number of global entries to sort
      integer gl_n0
      integer local_n_max     ! Largest local array size 
      character*128 str

      integer i,j,k,ip,il,jloop

      real tempr1,tempr2
      integer phase1_steps
      integer phase2_steps
      integer phase0_steps
      integer tempi

      integer nsum,nsum2
      integer nrecv,nsend
      integer ierr

      integer np2       ! Nearest power of 2 number of processors

!     communication
      integer ic
      integer irecv
      integer msg_id0, msg_id1, msg_id2, msg_id3  ! message id for non-blocking receive
      integer msg_id4                             ! message id for non-blocking receive
      integer srcid, dstid                        ! source and destination node id
      integer leng                                ! buffer size

!     temporary variables
      integer gntemp

!     Timing
      real*8 tim0,dnekclock

      alen = 2*lt                  ! Available array length
      local_n_max = iglmax(n,1)    ! Largest local array
!      alen2 = alen/2              ! Half of work array
      alen2 = alen - local_n_max   ! How much of the array we can use
!      alen2 = 500

      tempr1 = (alen2+0.0)/(local_n_max+0.0)

      tempi = int(log(np+0.0)/log(2.0))
      phase0_steps = iglmax(tempi,1)
      np2          = 2**(phase0_steps)
      if (np.gt.np2) phase0_steps = phase0_steps+1

!     No of steps till the array is full
      tempi = int(log(tempr1)/log(2.0))
      phase1_steps = iglmax(tempi,1)

!     No of steps in the second phase      
      tempi = np2/(2**phase1_steps)
      phase2_steps = iglmax(tempi,1)

!     We don't need phase 2 if phase 1 already
!     covers all the data/processors      
      if (phase1_steps.gt.phase0_steps) then
        phase1_steps = phase0_steps
        phase2_steps = 0
      endif  

!     Just writing out some stuff
!     Can be deleted.      
      call blank(str,128)
      write(str,'(A21,I7,A24,I7)') 
     $ 'Total Array Length = ',alen2,
     $ '; Max local Data size = ',local_n_max
      call mntr_log(fs_id,fs_log,str)

      call blank(str,128)
      write(str,'(A6,I5)') 
     $ 'NP2 = ',np2
      call mntr_log(fs_id,fs_log,str)

      call blank(str,128)
      write(str,'(A16,I5)') 
     $ 'Phase 1 Steps = ',phase1_steps
      call mntr_log(fs_id,fs_log,str)
      
      call blank(str,128)
      write(str,'(A16,I5)') 
     $ 'Phase 2 Steps = ',phase2_steps
      call mntr_log(fs_id,fs_log,str)


!     Initialize arrays      
      call rzero(wk1,alen)
      call rzero(wk2,alen)

      call rzero(ta1,alen)
      call rzero(ta2,alen)

      call copy(ta1,fld1,n)
      call copy(ta2,fld2,n)

!     Give elements an initial ownership 
      gl_n  = iglsum(n,1)           ! Total no of global elements
      gl_n0 = igl_running_sum(n)    ! Running sum
      gl_n0 = gl_n0 - n             ! entries before this nid.
      do i=1,n
        gnum2(i)  = gl_n0 + i       ! Assign a global number
        owner2(i) = nid             ! Assign current processor as owner

        gnum(i)   = gnum2(i) 
        owner0(i) = owner2(i)
      enddo  

!     Cummulative no of elements      
      nsum = n

!     First Phase      
!     Collect Data
!----------------------------------------
      tim0 = dnekclock() 
      do ic=2,phase1_steps

        il = 2**(ic-1)
        srcid = nid - il
        dstid = nid + il
        if (srcid.lt.0)  srcid = srcid + np2
        if (dstid.ge.np2) dstid = dstid - np2

!       set buffer for the number of elements to receive
        leng  = isize
        msg_id0 = irecv(0,nrecv,leng)

!       send local size of array
        nsend = nsum
        call csend(0,nsend,leng,dstid,0)

!       wait for length data
        call msgwait(msg_id0)

!       receive x1 data
        leng = wdsize*nrecv
        msg_id1 = irecv(1,ta1(nsum+1),leng)

!       receive x2 data
        leng = wdsize*nrecv
        msg_id2 = irecv(2,ta2(nsum+1),leng)

!       receive ownership data
        leng = isize*nrecv
        msg_id3 = irecv(3,owner2(nsum+1),leng)

!       receive ownership data
        leng = isize*nrecv
        msg_id4 = irecv(4,gnum2(nsum+1),leng)


!       Send data
        leng = wdsize*nsend
        call csend(1,ta1,leng,dstid,0)
        call csend(2,ta2,leng,dstid,0)
        leng = isize*nsend
        call csend(3,owner2,leng,dstid,0)
        call csend(4,gnum2,leng,dstid,0)

!       finish communication            
        call msgwait(msg_id1)
        call msgwait(msg_id2)
        call msgwait(msg_id3)
        call msgwait(msg_id4)

        nsum      = nsum + nrecv

        call nekgsync()
      enddo         ! ic

!     Error check
      i = 0 
      if (nsum.gt.alen2) then
        i = 1
      endif
      ierr = iglsum(i,1)
      if (ierr.gt.0) then   
!        if (nid.eq.0) write(6,*) 
!     $    'No of elements exceeded array size during collection in',
!     $    'fs_global_sort_mp', nsum, alen2
!        if (nid.eq.0) write(6,*) 'Exitting'

        call blank(str,128)
        write(str,'(A55,I7,A11,I7)') 
     $   'No of elements exceeded array size during collection N=',
     $      nsum,', Arr. Len=', alen2
        call mntr_log(fs_id,fs_log,str)
       
        call exitt
      endif  

!     Perform Local sort
!     Assuming data is continuous
!     Inputs:
!     ta1         - Primary field entries
!     ta2         - Secondary field entries
!     nsum        - no of entries
!     fdim        - No of fields (2/3)
!     Outputs:
!     ta1         - sorted entries (Primary field)
!     ta2         - sorted entries (Secondary field)
!     ind         - permutation indicies
!     ifseg       - segment breaks
!     nseg        - no of segments
!     ninseg      - entries in each segment
!     Work:
!     wk3         - real work array
!     wki         - integer work array
      fdim = ndim
      call fs_local_tuple_sort(ta1,ta2,ind,ifseg,nseg,
     $                               ninseg,nsum,fdim,wk3,wki)

!     Shuffle ownership and global numbering      
      call iswap_ip(owner2,ind,nsum)
      call iswap_ip(gnum2,ind,nsum)

!     Debugging      
!      do ic=0,np-1
!        if (nid.eq.ic) then
!          do i=1,nsum
!            write(6,*) nid, gnum2(i),ind(i),owner2(i)
!          enddo
!        endif
!        call nekgsync()
!      enddo  


!     Reassign ownership/Numbering. 
!     Take smallest global number.
!     And the corresponding nid.
      i = 1
      j = 1
      do iseg=1,nseg
        own = 999999999
        gntemp = 999999999
        do k=1,ninseg(iseg)
          if (gnum2(i).lt.gntemp) then
!          if (owner2(i).lt.own) then
            own = owner2(i)
            gntemp = gnum2(i)
          endif
          i = i+1
        enddo
        call ifill(owner2(j),own,ninseg(iseg))
        call ifill(gnum2(j),gntemp,ninseg(iseg))
        j = j + ninseg(iseg)
      enddo  

!     Reassign ownership/numbering of local elements
      do i=1,nsum
        j = ind(i)     ! j = Original position in array
                       ! i = Curent position in array
!       Ownership of local point        
        if (j.le.n) owner0(j) = owner2(i)     ! local point
!       Global numbering of local point        
        if (j.le.n) gnum(j)  = gnum2(i)       ! local point
      enddo  

      tim0 = dnekclock() - tim0

      call blank(str,128)
      write(str,'(A41,E14.8)') 
     $ 'Global Sorting: Phase I complete. Time = ',tim0
      call mntr_log(fs_id,fs_log,str)

!     Second Phase
!----------------------------------------
      tim0 = dnekclock()
      do ic=2,phase2_steps

        il = ic*(2**phase1_steps)
        srcid = nid - il
        dstid = nid + il
        if (srcid.lt.0)  srcid = srcid + np2
        if (dstid.ge.np2) dstid = dstid - np2

!       set buffer for the number of elements to receive
        leng  = isize
        msg_id0 = irecv(0,nrecv,leng)

!       send local size of array
        nsend = nsum
        call csend(0,nsend,leng,dstid,0)

!       wait for length data
        call msgwait(msg_id0)

!       Put local fields into work array
        call copy(wk1,fld1,n)
        call copy(wk2,fld2,n)
        call icopy(owner1,owner0,n)
        call icopy(gnum1,gnum,n)

!       receive x1 data
        leng    = wdsize*nrecv
        msg_id1 = irecv(1,wk1(n+1),leng)

!       receive x2 data
        leng    = wdsize*nrecv
        msg_id2 = irecv(2,wk2(n+1),leng)

!       receive ownership data
        leng    = isize*nrecv
        msg_id3 = irecv(3,owner1(n+1),leng)

!       receive global numbering data
        leng    = isize*nrecv
        msg_id4 = irecv(4,gnum1(n+1),leng)


!       Send data
        leng = wdsize*nsend
        call csend(1,ta1,leng,dstid,0)
        call csend(2,ta2,leng,dstid,0)
        leng = isize*nsend
        call csend(3,owner2,leng,dstid,0)
        call csend(4,gnum2,leng,dstid,0)

!       finish communication            
        call msgwait(msg_id1)
        call msgwait(msg_id2)
        call msgwait(msg_id3)
        call msgwait(msg_id4)

        nsum2 = n + nrecv

!       Perform Local sort
        fdim = ndim
        call fs_local_tuple_sort(wk1,wk2,ind,ifseg,nseg,
     $                                 ninseg,nsum2,fdim,wk3,wki)

!       Shuffle ownership and global numbering      
        call iswap_ip(owner1,ind,nsum2)        ! Ownership
        call iswap_ip(gnum1,ind,nsum2)         ! Global numbering

!       Reassign ownership/Numbering. 
!       Take smallest global number.
!       And the corresponding nid as ownership.
        i = 1
        j = 1
        do iseg=1,nseg
          own = 999999999
          gntemp = 999999999
          do k=1,ninseg(iseg)
            if (gnum1(i).lt.gntemp) then
              own = owner1(i)
              gntemp = gnum1(i)
            endif
            i = i+1
          enddo
          call ifill(owner1(j),own,ninseg(iseg))
          call ifill(gnum1(j),gntemp,ninseg(iseg))
          j = j + ninseg(iseg)
        enddo  

!       Reassign local ownership
        do i=1,nsum2
          j = ind(i)     ! j = Original position in array
!         Ownership of local point        
          if (j.le.n) owner0(j) = owner1(i)      ! local point
!         Global numbering of local point        
          if (j.le.n) gnum(j)  = gnum1(i)       ! local point
        enddo  

      enddo         ! ic=2,phase2_steps

      tim0 = dnekclock() - tim0

      call blank(str,128)
      write(str,'(A42,E14.8)') 
     $ 'Global Sorting: Phase II complete. Time = ',tim0
      call mntr_log(fs_id,fs_log,str)

      return
      end subroutine fs_global_sort_mp
!---------------------------------------------------------------------- 

