      program psf

      implicit none

      integer NrMol,ApM
      logical bnds
      parameter (ApM   = 4)
      parameter (bnds  = .true.)

      integer NrUit
      integer i,j,k,l
      read (*,*) NrMol
      NrUit=NrMol*ApM

      open (20,file='config.psf',status='new',form='formatted')
      write (20,'(a3)') 'PSF'
      write (20,*)

      write (20,'(i8,a8)') 1,'!NTITLE'
      write (20,'(a12)') 'REMARKS etc'
      write (20,*)

      write (20,'(i8,a14)') NrUit,'!NATOM: atoms'
      l=0
      do i=1,NrMol
        do j=1,ApM
          l=l+1
          write (20,9000) l,1,'UNK',"H",0,"H",0,0.0,0.0
        enddo ! j
      enddo ! i
      write (20,*)

      if (bnds) then
        write (20,'(i8,a14)') NrMol*ApM,'!NBOND: bonds'
        write (20,'(8i8)') (l*ApM+1,l*ApM+2,
     x                      l*ApM+1,l*ApM+3,
     x                      l*ApM+3,l*ApM+4,
     x                      l*ApM+2,l*ApM+4, l=0,NrMol-1)
      else
        write (20,'(i8,a14)') 0,'!NBOND: bonds'
      endif
      write (20,*)

      write (20,'(i8,a16)') 0,'!NTHETA: angles'
      write (20,*)

      write (20,'(i8,a17)') 0,'!NPHI: dihedrals'
      write (20,*)

      write (20,'(i8,a19)') 0,'!NIMPHI: impropers'
      write (20,*)

      write (20,'(i8,a14)') 0,'!NDON: donors'
      write (20,*)

      write (20,'(i8,a17)') 0,'!NACC: acceptors'
      write (20,*)

      write (20,'(i8,a28)') 0,'!NNB: non-bonded exclusions'
      write (20,*)

*     do i=1,797 ! why?
*       write (20,'(8i8)') (0, j=1,8)
*     enddo ! i
*     write (20,9000) 0
*     write (20,*)

      write (20,'(2i8,a15)') 0,0,'!NGRP: ignored'
      write (20,*)

      close (20)

 9000 format(i8,i7,a7,2(a3,i2),f12.6,f14.4)
      end
