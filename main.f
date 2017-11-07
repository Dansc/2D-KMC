c---------------------------------------------------------------------
c
c This program is used to perform monte carlo trial moves on
c 4'4-biphenyldicarboxyclic molecules on a square Copper 100
c lattice 
c Daniel Schwarz  2011
c
c---------------------------------------------------------------------
  
      program main

      implicit none

c---------------------------------------------------------------------
c
c Declaration of variables and parameters
c
c---------------------------------------------------------------------

      integer   i,j,k,l
      integer   jj
      integer   x,xx,y,yy,xi,xj,yi,yj
      integer   diffx, diffy, orient, orientA, orientB
      integer   PosLst(0:1,3,12)! (Orientation, (x,y,kind),Aindex)
      real      Earray(5,5)

      integer   TotalMol, StartMol, ActualMol, rate
      integer   TotalSize
      integer   Grids,dx
      integer   MaxPerGrid
      integer   GridSize
      integer   StepSize
      real      Size
      integer   continued
      integer*8 captured, trialsdone, Totalsteps, MC_step, AccSteps 
      integer*8 addmol

      integer   xold,yold
      integer   xnew,ynew
      integer   xtemp, ytemp, rotemp
      integer   xgrdtmp, ygrdtmp
      integer   xgrdold, ygrdold
      integer   xgrdnew, ygrdnew
      integer   posx, posy
      integer   molecule, move, possible

      real  density, startdensity, rho, test
      real  randy
    
      real  deltaEE, deltaEmax, expDeltaEmax
      real  Einteract, TotalEnergy
      real  Etemp, MolEnergy
      real  Eold, Enew, deltaE
      real  Eoverlap
      parameter (Eoverlap = 999.)
      real,   allocatable :: Energy(:,:,:,:)
      real*4, allocatable ::    xprint(:,:)
      real*4    zero4
      parameter ( zero4 = 0. )
    
      integer iseed(4)
    
      integer, allocatable :: MoList(:,:)
      integer, allocatable :: GrdLst(:,:,:)
      integer, allocatable :: PBCoff(:)
      integer, allocatable :: GrdNghbr(:)


c       parameters of the system
      parameter (deltaEmax = 20.)

c       derived variabels
      parameter (expDeltaEmax = exp(-deltaEmax))


c-----------------------------------------------------------------------
c 
c   1. Load simulation parameters and generate random start
c   condition
c
c-----------------------------------------------------------------------
 

      open (21,file='input.control', status='old', form='formatted')          
      open (98,file='output.ene', form='formatted')
      open (99,file='output.dcd', form='unformatted')
    
      read (21,*) iseed(1),iseed(2),iseed(3),iseed(4)
      read (21,*) TotalSteps
      read (21,*) Stepsize
      read (21,*) Grids
      read (21,*) Gridsize
      read (21,*) density
      read (21,*) startdensity
      read (21,*) rate
      read (21,*) MaxPerGrid
      read (21,*) Einteract
      close (21)


      TotalSize      = GridSize*Grids
      TotalMol       = ceiling(Grids*Grids*density)
      StartMol       = ceiling(Grids*Grids*startdensity)
      ActualMol      = StartMol
      rho            = startdensity
      Size           = grids*grids


      AccSteps = 0
      addmol = 0

      i = 2*GridSize+2
      allocate(GrdNghbr(0:Grids+1))
      allocate(PBCoff(0:Grids+1))
      allocate(MoList(3,TotalMol))
      allocate(GrdLst(0:MaxPerGrid, Grids, Grids))
      allocate(Energy(-i:i,-i:i,0:1,0:1))
      allocate(xprint(2,4*TotalMol))

      if (Totalsteps.lt.0) then ! -100 op input = maak config + 100 stappen
        call rstart(iseed) ! of ReadRandstate voor continuatie
       ! maak random start configuratie zoals hieronder
        Totalsteps=-Totalsteps
        trialsdone=0
        continued=0

      else

        open (20,file='input.config',status='old',form='unformatted')
        read (20) i,trialsdone
        if (i.ne.StartMol) stop 'mismatch in particle numbers'
        read (20) dx
        if (dx.ne.Grids) stop 'mismatch in box dimension'
        read (20) ( (MoList(i,k), i=1,3), k=1,StartMol )
        close (20)
        if (trialsdone.eq.0) then
            call rstart(iseed)
        else
            call ReadRandState64
        endif
        TotalSteps=trialsdone+TotalSteps
        continued = 1
c        ! nog energie berekenen voor we verder kunnen
      endif 

         ! Write out parameters of simulation to stdout
      write (*,*) "Size of system:            ", TotalSize*TotalSize
      write (*,*) "Number of molecules:       ", TotalMol
      write (*,*) "TotalSteps:            ", TotalSteps
      write (*,*) "Steps per frame:           ", StepSize
       
       ! At each step we write out parameters of the simulation
      captured = trialsdone + Stepsize
c-----------------------------------------------------------------------
c   1.1 Initialize GridList
c-----------------------------------------------------------------------
 
      do i= 1,Grids
        do j= 1,Grids
            GrdLst(0,j,i) = 0
        enddo !j
      enddo !i

c-----------------------------------------------------------------------
c   1.2 Generate list of neighbouring cells -
c   Periodic Boundary Condition!
c-----------------------------------------------------------------------

      do i=1,Grids
        PBCoff(i) = 0
      enddo !i
      
      PBCoff(0)         = -TotalSize
      PBCoff(Grids+1)   =  TotalSize

      do i=1,Grids
        GrdNghbr(i)     = i
      enddo !i
      GrdNghbr(0)       = Grids
      GrdNghbr(Grids+1) = 1

c-----------------------------------------------------------------------
c   1.3 Generate list of atom positions inside molecule for 
c   two possible orientations
c   Can be generalized!
c-----------------------------------------------------------------------
    

      do i = 1,6
        PosLst(0,1,i) = i -3
        PosLst(0,2,i) = 1
      enddo !i         
    
      do i = 7,12
        PosLst(0,1,i) = i -6 - 3
        PosLst(0,2,i) = 0
      enddo !i    

      do i = 1,6
        PosLst(1,1,i) = -1
        PosLst(1,2,i) = i-3
      enddo !i         
    
      do i = 7,12
        PosLst(1,1,i) = 0
        PosLst(1,2,i) = i -6 -3
      enddo !!
    
    
      do orient = 0,1
        PosLst(orient,3,1)  = 1
        PosLst(orient,3,2)  = 2
        PosLst(orient,3,3)  = 3
        PosLst(orient,3,4)  = 3
        PosLst(orient,3,5)  = 2
        PosLst(orient,3,6)  = 1
        PosLst(orient,3,7)  = 1
        PosLst(orient,3,8)  = 2
        PosLst(orient,3,9)  = 3
        PosLst(orient,3,10) = 3
        PosLst(orient,3,11) = 2
        PosLst(orient,3,12) = 1
      enddo !orient

      PosLst(1,3,1)       = 4
      PosLst(1,3,6)       = 4
      PosLst(1,3,7)       = 4
      PosLst(1,3,12)      = 4
    
      PosLst(1,3,3)       = 5
      PosLst(1,3,4)       = 5
      PosLst(1,3,9)       = 5
      PosLst(1,3,10)      = 5


c-----------------------------------------------------------------------
c   1.4 Generate Array with interaction energies for all atom 
c   types. Only 5-1 and 4-3 (Benzolring - Oxygen) bind!
c-----------------------------------------------------------------------

      Earray(1:4,1) = 0
      Earray(5,1)   = Einteract
    
      Earray(1:4,2) = 0
      Earray(5,2)   = 0

      Earray(1:3,3) = 0
      Earray(4,3)   = Einteract
      Earray(5,3)   = 0

      Earray(1:2,4) = 0
      Earray(3,4)   = Einteract
      Earray(4:5,4) = 0
    
      Earray(1,5)   = Einteract
      Earray(2:5,5) = 0


c-----------------------------------------------------------------------
c   1.5 Calculate the Energies for all configurations of two
c   molecules within reach of two neighboring Gridcells
c-----------------------------------------------------------------------

      do orientA = 0,1
      do orientB = 0,1
      do diffx = -16, 16
      do diffy = -16, 16
        Etemp = 0
        do i = 1,12 ! i is fixed molecule CMS is in 0,0
            do j = 1,12
                xi = PosLst(orientA,1,i)
                yi = PosLst(orientA,2,i)
                xj = PosLst(orientB,1,j) + diffx
                yj = PosLst(orientB,2,j) + diffy

                if ( (abs(xi-xj).eq.1).AND.((yi-yj).eq.0) ) then
                    Etemp = Etemp+Earray(PosLst(orientA,3,i),
     x               PosLst(orientB,3,j))
                endif

                if ( ((xi-xj).eq.0).AND.(abs(yi-yj).eq.1) ) then
                    Etemp = Etemp+Earray(PosLst(orientA,3,i),
     x               PosLst(orientB,3,j))
                endif    
                if ( (xi.eq.xj).AND.(yi.eq.yj) ) then    !Overlap!
                    Etemp = Eoverlap
                    goto 5
                endif
            enddo !j
        enddo !i
 5    continue ! end goto if Overlap
      Energy(diffx,diffy,orientA,orientB) = Etemp
c   if ( Etemp.LT.(Eoverlap+0.1*Eoverlap)) then
c   Energy(diffx,diffy,orientA,orientB) = 0
c   endif

      enddo !y
      enddo !x
      enddo !orientA
      enddo !orientB
    

c-----------------------------------------------------------------------
c   1.6 Fill the lattice up to TotalMol
c   xtemp, ytemp - Position of Mol - square molecule 2x2 with
c   Pos defined as bottom left lattice site
c      
c   rotemp  - 0: lying along x axis - 1: lying along y axis (only
c   needed for elongated molecules
c-----------------------------------------------------------------------

      Totalenergy = 0

      if (continued.eq.1) goto 11     !only do this if starting new sim


      do i = 1, StartMol

 10     xtemp = (Totalsize*randy())+1
        ytemp = (Totalsize*randy())+1
        xgrdtmp = Int((xtemp-1)/GridSize) +1
        ygrdtmp = Int((ytemp-1)/GridSize) +1
        orientA = (2*randy())
        Etemp = 0
        MolEnergy = 0
        do x = xgrdtmp-1,xgrdtmp+1
            xx = GrdNghbr(x)
        do y = ygrdtmp-1,ygrdtmp+1
            yy = GrdNghbr(y)
            do j= 1,GrdLst(0,xx,yy)
               jj      = GrdLst(j,xx,yy)
               posx    = MoList(1,jj) + PBCoff(x)
               posy    = MoList(2,jj) + PBCoff(y)
               orientB = MoList(3,jj)
               
               diffx   = posx-xtemp
               diffy   = posy-ytemp

           if (Energy(diffx,diffy,orientA,orientB).GT.Eoverlap-10) then
                goto 10
              endif

                Etemp       = Energy(diffx,diffy,orientA,orientB) 
                Totalenergy = Totalenergy + Etemp
                MolEnergy   = MolEnergy + Etemp
            enddo !j
        
        enddo !y
        enddo !x
        
        MoList(1,i) = xtemp
        MoList(2,i) = ytemp
        MoList(3,i) = orientA
        
        GrdLst(0,xgrdtmp,ygrdtmp) = GrdLst(0,xgrdtmp,ygrdtmp)+1

        if (Grdlst(0,xgrdtmp,ygrdtmp).gt.MaxPerGrid) then
            stop 'too many particles in a cell'
        endif
        
        GrdLst(GrdLst(0,xgrdtmp,ygrdtmp),xgrdtmp,ygrdtmp) = i
     
      enddo ! i
    
      do i = StartMol+1,TotalMol
        MoList(1,i) = -200
        MoList(2,i) = -200
        MoList(3,i) = 0
      enddo !i
      goto 12


c-----------------------------------------------------------------------
c   
c   Calculate Energy and fill Gridlist if continuation of old config    
c   
c-----------------------------------------------------------------------



 11   do i = 1,StartMol
        xtemp   = MoList(1,i)
        ytemp   = MoList(2,i)
        xgrdtmp = Int((xtemp-1)/GridSize) +1
        ygrdtmp = Int((ytemp-1)/GridSize) +1
        GrdLst(0,xgrdtmp,ygrdtmp) = GrdLst(0,xgrdtmp,ygrdtmp) +1

        if (GrdLst(0,xgrdtmp,ygrdtmp).gt.MaxPerGrid) then
            stop 'too many particles in a cell!'
        endif

        GrdLst(GrdLst(0,xgrdtmp,ygrdtmp),xgrdtmp,ygrdtmp) = i

      enddo !i
    

      do i = 1,StartMol
        xtemp   = MoList(1,i)
        ytemp   = MoList(2,i)
        orientA = MoList(3,i)
        xgrdtmp =  Int((xtemp-1)/GridSize) +1
        ygrdtmp =  Int((ytemp-1)/GridSize) +1
        Etemp = 0
        MolEnergy = 0
          do x = xgrdtmp-1,xgrdtmp+1
            xx = GrdNghbr(x)
              do y = ygrdtmp-1,ygrdtmp+1
                yy = GrdNghbr(y)
                  do j= 1,GrdLst(0,xx,yy)
                    jj = GrdLst(j,xx,yy)
                    if (jj.NE.i) then
                       posx    = MoList(1,jj) + PBCoff(x)
                       posy    = MoList(2,jj) + PBCoff(y)
                       orientB = MoList(3,jj)
                       
                       diffx   = posx-xtemp
                       diffy   = posy-ytemp

            if (Energy(diffx,diffy,orientA,orientB).GT.Eoverlap-1) then
                        stop 'Error in input configuration, overlap!'
                    endif
                    Etemp       = Energy(diffx,diffy,orientA,orientB) 
                    Totalenergy = Totalenergy + Etemp
                    MolEnergy   = MolEnergy + Etemp
                    endif !
                enddo !j
            enddo !y
          enddo !x
     
      enddo ! i

      do i = StartMol+1,TotalMol
        MoList(1,i) = -200
        MoList(2,i) = -200
        MoList(3,i) = 0
      enddo !i


 12   write (98,*) trialsdone, Totalenergy, (rho*ActualMol), ActualMol
      write (99) 'CORD',9999,        ! header, write once only
     x   0,1,0,0,0,0,0,0,1.0d-15,
     x   0,0,0,0,0,0,0,0,0
      write(99) 0                    ! zero lines of comments
      write(99) 4*TotalMol           ! particles in movie: 2 per BDA
      j = 0

      do i = 1, TotalMol
        orient = MoList(3,i)
        xprint(1,j+1) = real( MoList(1,i) + PosLst(orient,1,1) ) 
        xprint(1,j+2) = real( MoList(1,i) + PosLst(orient,1,6) )
        xprint(1,j+3) = real( MoList(1,i) + PosLst(orient,1,7) )
        xprint(1,j+4) = real( MoList(1,i) + PosLst(orient,1,12) )

        xprint(2,j+1) = real( MoList(2,i) + PosLst(orient,2,1) ) 
        xprint(2,j+2) = real( MoList(2,i) + PosLst(orient,2,6) )
        xprint(2,j+3) = real( MoList(2,i) + PosLst(orient,2,7) )
        xprint(2,j+4) = real( MoList(2,i) + PosLst(orient,2,12) )
        j = j + 4
      enddo !i

      write(99) (xprint(1,j), j=1, 4*TotalMol)
      write(99) (xprint(2,j), j=1, 4*TotalMol)
      write(99) (zero4, j=1,4*TotalMol)            ! z = 0.

c-----------------------------------------------------------------------
c
c   2. Do Monte Carlo move
c
c-----------------------------------------------------------------------

      do MC_step = (trialsdone+1), TotalSteps

c-----------------------------------------------------------------------
c
c   2. Deposit new molecule
c
c-----------------------------------------------------------------------


      if (AccSteps.eq.addmol) then
        if (ActualMol.eq.TotalMol) goto 41
40      xtemp = (Totalsize*randy())+1
        ytemp = (Totalsize*randy())+1
        xgrdtmp = Int((xtemp-1)/GridSize) +1
        ygrdtmp = Int((ytemp-1)/GridSize) +1
        orientA = (2*randy())
        Etemp = 0
        MolEnergy = 0
        do x = xgrdtmp-1,xgrdtmp+1
            xx = GrdNghbr(x)
            do y = ygrdtmp-1,ygrdtmp+1
                yy = GrdNghbr(y)
                    do j= 1,GrdLst(0,xx,yy)
                       jj      = GrdLst(j,xx,yy)
                       posx    = MoList(1,jj) + PBCoff(x)
                       posy    = MoList(2,jj) + PBCoff(y)
                       orientB = MoList(3,jj)
                       diffx   = posx-xtemp
                       diffy   = posy-ytemp

        if (Energy(diffx,diffy,orientA,orientB).GT.Eoverlap-1) then
                            goto 40 ! No Molecules deposited ontop of each other
                        endif
                    Etemp = Energy(diffx,diffy,orientA,orientB) 
                    Totalenergy = Totalenergy + Etemp
                    MolEnergy = MolEnergy + Etemp
                    enddo !j
            enddo !y
        enddo !x


        MoList(1,ActualMol+1) = xtemp
        MoList(2,ActualMol+1) = ytemp
        MoList(3,ActualMol+1) = orientA
        GrdLst(0,xgrdtmp,ygrdtmp) = GrdLst(0,xgrdtmp,ygrdtmp)+1

        if ( abs(MolEnergy).LT.(-2*Einteract) ) then
        rho = (exp(-1./ActualMol))*rho +(1.-exp(-1./ActualMol))*(1)
        else
            rho = (exp(-1./ActualMol))*rho
        endif

        if (Grdlst(0,xgrdtmp,ygrdtmp).gt.MaxPerGrid) then
           stop 'too many particles in a cell'
        endif
        GrdLst(GrdLst(0,xgrdtmp,ygrdtmp),xgrdtmp,ygrdtmp) = ActualMol+1
        ActualMol = ActualMol + 1
        addmol = addmol + (rho*ActualMol*rate)*randy()
        goto 20
      endif 


c-----------------------------------------------------------------------
c
c   2. Move molecule at random
c
c-----------------------------------------------------------------------

41    continue
      Eold = 0
      Enew = 0
      molecule = ActualMol*randy()+1       ! select molecule to move
    
      xold = MoList(1,molecule)
      yold = MoList(2,molecule)
      orientA = MoList(3,molecule)

      xgrdold = Int((xold -1 )/GridSize) +1
      ygrdold = Int((yold -1 )/GridSize) +1

      move = 4*randy() + 1

      goto (51,52,53,54), move

 51   xnew = xold -1
      ynew = yold

      goto 59

 52   xnew = xold +1
      ynew = yold

      goto 59

 53   ynew = yold +1
      xnew = xold

      goto 59

 54   ynew = yold -1
      xnew = xold

      goto 59

 59   continue ! end goto

      if( xnew.LE.0 ) xgrdnew = 0
      if( xnew.GT.0)  xgrdnew = Int((xnew-1)/Gridsize) +1
      xtemp = xnew - PBCoff(xgrdnew)
      xgrdnew = GrdNghbr(xgrdnew)

      if( ynew.LE.0 ) ygrdnew = 0
      if( ynew.GT.0)  ygrdnew = Int((ynew-1)/Gridsize) +1
      ytemp = ynew - PBCoff(ygrdnew)
      ygrdnew = GrdNghbr(ygrdnew)

c-----------------------------------------------------------------------
c   2.1 0 - move left, 1 - move right, 2 - move up, 3 - move down
c-----------------------------------------------------------------------
 

c-----------------------------------------------------------------------
c
c   2.2 Is move possible?
c
c-----------------------------------------------------------------------

      do x = xgrdold-1,xgrdold+1
        xx = GrdNghbr(x)
        do y = ygrdold-1,ygrdold+1
            yy = GrdNghbr(y)

            do j= 1,GrdLst(0,xx,yy)
                jj = GrdLst(j,xx,yy)
                if (jj.NE.molecule) then
                    posx    = MoList(1,jj) + PBCoff(x)
                    yi      = posy 
                    posy    = MoList(2,jj) + PBCoff(y)
                    orientB = MoList(3,jj)
                    !continue if below
c-----------------------------------------------------------------------
c   2.3 First calculate energy before move
c-----------------------------------------------------------------------

                    Eold = Eold + 
     x               Energy(posx-xold,posy-yold,orientA,orientB)

c-----------------------------------------------------------------------
c   2.4 Calculate energy after move, if overlap go to next MC step
c-----------------------------------------------------------------------

                    Etemp = Energy(posx-xnew,posy-ynew,orientA,orientB)
      
                    if ( Etemp.GT.Eoverlap-1) goto 20 !Overlap in new pos!
    
                    Enew = Enew + Etemp

                endif !! jj not equal molecule
    
            enddo !j
        enddo !y
      enddo !x

c-----------------------------------------------------------------------
c   2.5 Energy of new state lower? Do the move
c       else: roll the dice!
c-----------------------------------------------------------------------
      if ( abs(Eold).LT.(-2*Einteract) ) then    ! Update molecular density
        rho = (exp(-1./ActualMol))*rho +(1.-exp(-1./ActualMol))
      else
        rho = (exp(-1./ActualMol))*rho
      endif

      deltaE = Enew - Eold

      if (deltaE.LE.0) goto 19    ! move molecule
      deltaEE = deltaE
      do while (deltaEE.GT.deltaEmax)
        if (randy().GT.expDeltaEmax) goto 20 ! reject -> next mcstep
        deltaEE = deltaEE - deltaEmax
      enddo ! while
      if (randy().GT.exp(-deltaEE)) goto 20 ! reject

c-----------------------------------------------------------------------
c   2.6 Delete old information in GrdLst
c-----------------------------------------------------------------------        

19    continue

      xnew    = xtemp !Periodic Boundary Cond.
      ynew    = ytemp !Periodic Boundary Cond.
    
    
      if( (xgrdnew.NE.xgrdold).OR.(ygrdnew.NE.ygrdold) ) then 
! Find the molecule in the old gridlst, only done when it changed grid!
        do i = 1, Grdlst(0,xgrdold,ygrdold)
            if (Grdlst(i,xgrdold,ygrdold).EQ.molecule) jj = i
        enddo !i

        do i = jj, GrdLst(0,xgrdold,ygrdold)-1
            GrdLst(i,xgrdold,ygrdold) = GrdLst(i+1,xgrdold,ygrdold)
        enddo !i

        GrdLst(0,xgrdold,ygrdold) = GrdLst(0,xgrdold,ygrdold) -1
        GrdLst(0,xgrdnew,ygrdnew) = GrdLst(0,xgrdnew,ygrdnew) +1

c-----------------------------------------------------------------------
c   2.7 Moving the molecule
c-----------------------------------------------------------------------

        if ( GrdLst(0,xgrdnew,ygrdnew) .GT. MaxPerGrid) then 
            stop 'Too many particles in cell'
        endif !
        GrdLst(GrdLst(0,xgrdnew,ygrdnew),xgrdnew,ygrdnew) = molecule
    
      endif !
    
      MoList(1,molecule) = xnew
      MoList(2,molecule) = ynew

      Totalenergy = Totalenergy + deltaE

      AccSteps = AccSteps +1
    
 20   continue  ! continue without changing molist / grdlst

c-----------------------------------------------------------------------
c   
c   2.8 Write output every Stepsize's step
c   Step / Totalenergy into output.ene
c   Positions of all molecules in output.dcd
c   
c-----------------------------------------------------------------------

      if (MC_step.eq.captured) then
        captured = captured + stepsize   ! for next time
        write (*,*) "Step:    ",captured," out of:",TotalSteps
        write (*,*) "Totalenergy:     ",Totalenergy
        write (*,*) "2D Gas molecules: ",rho*ActualMol/(Size)
        write (*,*) "Molecules on surface: ",ActualMol

        j = 0
        do i = 1, TotalMol
            orient = MoList(3,i)
            xprint(1,j+1) = real( MoList(1,i) + PosLst(orient,1,1) )
            xprint(1,j+2) = real( MoList(1,i) + PosLst(orient,1,6) )
            xprint(1,j+3) = real( MoList(1,i) + PosLst(orient,1,7) )
            xprint(1,j+4) = real( MoList(1,i) + PosLst(orient,1,12) )

            xprint(2,j+1) = real( MoList(2,i) + PosLst(orient,2,1) )
            xprint(2,j+2) = real( MoList(2,i) + PosLst(orient,2,6) )
            xprint(2,j+3) = real( MoList(2,i) + PosLst(orient,2,7) )
            xprint(2,j+4) = real( MoList(2,i) + PosLst(orient,2,12) )
            j = j + 4
        enddo ! 1
    
        write(98,*) MC_step, Totalenergy, (rho*ActualMol), ActualMol
        write(99) (xprint(1,j), j=1, 4*TotalMol)
        write(99) (xprint(2,j), j=1, 4*TotalMol)
        write(99) (zero4, j=1,4*TotalMol)            ! z = 0.

      endif ! MC step captured
      enddo !MC_step

    
c-----------------------------------------------------------------------
c
c   3. Write Output
c
c-----------------------------------------------------------------------

c   do i = 1, 12
c       write(*,*) PosLst(0,1,i), PosLst(0,2,i) 
c   enddo !i

    
c   do i = -14, 14
c       do j = -14,14
c       write(*,*) Energy(i,j,0,0), i, j    
c       enddo !j
c   enddo !i


    
 9999 open  (90,file='output.config',form='unformatted')!status='old',
      write (90) TotalMol,MC_step
      write (90) Grids
      write (90) ( (MoList(i,k), i=1,3), k=1,TotalMol )
      close (90)
      call WriteRandState

      close(98)
      close(99)

      end program main


      real   function randy()

c     start_doc
c     Name
c       randy
c       Adapted by Johan Padding for use with real*8's.
c
Cmr   Adapted by Marten Renes for 64 bit flat distribution
c
c     Usage
c       integer i, j, k, l
c       real*8 random, randy
c       call rstart (i, j, k, l) ! seed the generator once and for all
c       random = randy()
c
c     Description
c       Universal portable random number generator.  Uniformly
c       distributed random numbers in [0, 1).
c       Ref. G Marsaglia, A Zaman and W W Tsang, Statistics and
c       Probability Letters 8, 35 (1990)
c
c       Correct performance of this implementation was tested using the
c       test driver of table 3 of the reference.
c
Cmr     Changed c,cd,cm to allow intepretation of 64bit accurate random 
Cmr     numbers
c
c     Parameters
c       i, j, k, l   integer seeds: i, j, k: 1-178 and not all 1
c                                   l:       0-168
c     Memo
c       random number generator
c
c     Imported
c
c     Author
c       fmp nov 90, hacked from above reference
c
c     end_doc

      implicit real   (a-h, o-z)
      real   u(97)
      common / crandy / u, c, cd, cm, ip, jp

      randy = u(ip) - u(jp)
      if (randy .lt. 0.d0) randy = randy + 1.d0
      u(ip) = randy

      ip = ip - 1
      if (ip .eq. 0) ip = 97
      jp = jp - 1
      if (jp .eq. 0) jp = 97

      c = c - cd
      if (c .lt. 0.0) c = c + cm

      randy = randy - c
      if (randy .lt. 0.d0) randy = randy + 1.d0

      end

c----------------------------------------------------------------------

      subroutine rstart (iseed)

      implicit real   (a-h, o-z)
      integer iseed(4)
      real   u(97)
      common / crandy / u, c, cd, cm, ip, jp

      i = iseed(1)
      j = iseed(2)
      k = iseed(3)
      l = iseed(4)

      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod (mod (i * j, 179) * k, 179)
            i = j
            j = k
            k = m
            l = mod (53 * l + 1, 169)
            if (mod (l * m, 64) .ge. 32) s = s + t
            t = 0.5 * t
 3       continue
         u(ii) = s
 2    continue

Cmr
      c  = 0. !  362436.0 / 16777216.0
      cd = 362 436 069 876.        /9 007 199 254 740 992.
           !7654321.0 / 16777216.0
      cm =   9 007 199 254 740 881./9 007 199 254 740 992.
           !16777213.0 / 16777216.0
Cmr   

      ip = 97
      jp = 33

      end

c----------------------------------------------------------------------

      subroutine ReadRandstate64

      implicit real   (a-h, o-z)
      real   u(97)
      common / crandy / u, c, cd, cm, ip, jp

      logical ex
   
      inquire(file='input.rand',exist=ex)
      if (ex) then
         open(unit=40,file='input.rand',status='old',form='unformatted')
         read(40) u,c,cd,cm,ip,jp
         close(40)
         print *,'Finished ReadRandState'
      else

        stop 'File input.rand niet gevonden'

      endif

      end

c----------------------------------------------------------------------


      subroutine WriteRandstate

      implicit real   (a-h, o-z)
      real   u(97)
      common / crandy / u, c, cd, cm, ip, jp

      open(unit=40,file='output.rand',status='unknown',
     x     form='unformatted')
      write(40) u,c,cd,cm,ip,jp
      close(40)

      end

c----------------------------------------------------------------------

      REAL   FUNCTION GAUSS ( )

C    *******************************************************************
C    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
C    **                                                               **
C    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
C    **    ADDISON-WESLEY), 1978                                      **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION randy()                                         **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    *******************************************************************

      REAL        A1, A3, A5, A7, A9
      PARAMETER ( A1 = 3.949846138d0, A3 = 0.252408784d0 )
      PARAMETER ( A5 = 0.076542912d0, A7 = 0.008355968d0 )
      PARAMETER ( A9 = 0.029899776d0                   )

      REAL        SUM, R, R2
      REAL        randy
      INTEGER     I

C    *******************************************************************

      SUM = 0.d0

      DO 10 I = 1, 12

      SUM = SUM + randy()

10    CONTINUE

      R  = ( SUM - 6.d0 ) / 4.d0
      R2 = R * R

      GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 )
     :          * R

      END


