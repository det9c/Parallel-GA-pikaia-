      program xpkaia
c======================================================================
c     Sample driver program for pikaia.f
c======================================================================
      implicit none
      integer n, seed, i, status,ndim,nparams
      parameter (n=6)
      double precision ctrl(12), x(n), f, twod
      external twod

c***********
c**********
c***********
c**********
c parallel stuff
c***********
c**********
      common /parallel/ndim1,npopsize,nparams,myrank,rankno
      integer myrank,ierr,nprocs,npopsize,ndim1
      include 'mpif.h'
      character(50) word 
      character(4) rankno(10000)


	call system('date')
     
 



      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      call mpi_comm_size(mpi_comm_world,nprocs,ierr)


c     Set control variables (use defaults)
      do 10 i=1,12
         ctrl(i) = -1
 10       continue

       npopsize=128
       ctrl(1)=npopsize
       ctrl(2)=1000
       ctrl(11)=0
       nparams=n

       call slavdrv
        if(myrank.eq.0)then
          open(unit=24,file='genwrite')
          write(24,*)0
          close(24)
        end if


       call mpi_barrier(mpi_comm_world,ierr)




c make a directory for each processor
c      call system('module load g09/d01')
      word(1:6)='mkdir '
      word(7:)=rankno(myrank+1)
      call system(word)
      word=''
      word(1:8)='./xcopy '
      word(9:)=rankno(myrank+1)
      call system(word)



c    1 write(*,'(/A$)') ' Random number seed (I*4)? '
c      read(*,*) seed
c      seed=123456
      open(unit=23,file='seed')
      read(23,*)seed
      close(23)
      call rninit(seed)
c
c     Set control variables (use defaults)
c      do 10 i=1,12
c         ctrl(i) = -1
c   10 continue

c       npopsize=128
c       ctrl(1)=npopsize
c       ctrl(2)=2000
c       ctrl(11)=0
c       nparams=n

c       call slavdrv
       call mpi_barrier(mpi_comm_world,ierr)
c     Now call pikaia
c       if(myrank.eq.0)then
      call pikaia(twod,n,ctrl,x,f,status)
c      end if
	call system('date')
      call mpi_finalize(ierr)
      end





c*********************************************************************
      function urand()
c=====================================================================
c     Return the next pseudo-random deviate from a sequence which is
c     uniformly distributed in the interval [0,1]
c
c     Uses the function ran0, the "minimal standard" random number
c     generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
c     Comm. ACM 36 No. 7, 105-110, July 1993).
c=====================================================================
      implicit none
c
c     Input - none
c
c     Output
      double precision     urand
c
c     Local
      integer  iseed
      double precision     ran0
      external ran0
c
c     Common block to make iseed visible to rninit (and to save
c     it between calls)
      common /rnseed/ iseed
c
      urand = ran0( iseed )
      return
      end
c*********************************************************************
      subroutine rninit( seed )
c=====================================================================
c     Initialize random number generator urand with given seed
c=====================================================================
      implicit none
c
c     Input
      integer seed
c
c     Output - none
c
c     Local
      integer iseed
c
c     Common block to communicate with urand
      common /rnseed/ iseed
c
c     Set the seed value
      iseed = seed
      if(iseed.le.0) iseed=123456
      return
      end
c*********************************************************************
      function ran0( seed )
c=====================================================================
c     "Minimal standard" pseudo-random number generator of Park and
c     Miller.  Returns a uniform random deviate r s.t. 0 < r < 1.0.
c     Set seed to any non-zero integer value to initialize a sequence,
c     then do not change seed between calls for successive deviates
c     in the sequence.
c
c     References:
c        Park, S. and Miller, K., "Random Number Generators: Good Ones
c           are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
c        Park, S. and Miller, K., in "Remarks on Choosing and Imple-
c           menting Random Number Generators", Comm. ACM 36 No. 7,
c           105-110 (July 1993)
c=====================================================================
c *** Declaration section ***
c
      implicit none
c
c     Input/Output:
      integer seed
c
c     Output:
      double precision ran0
c
c     Constants:
      integer A,M,Q,R
      parameter (A=48271,M=2147483647,Q=44488,R=3399)
      double precision SCALE,EPS,RNMX
      parameter (SCALE=1./M,EPS=1.2e-7,RNMX=1.-EPS)
c
c     Local:
      integer j
c
c *** Executable section ***
c
      j = seed/Q
      seed = A*(seed-j*Q)-R*j
      if (seed .lt. 0) seed = seed+M
      ran0 = min(seed*SCALE,RNMX)
      return
      end
c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      double precision      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      double precision      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end
c***********************************************************************
      subroutine pikaia(ff,n,ctrl,x,f,status)
c=======================================================================
c*******************************************

      implicit none
 
c     Input:
      integer   n
      double precision      ff
      external  ff


      common /parallel/ndim1,npopsize,nparams,myrank,rankno
      character(4) rankno(10000)
      integer myrank,ierr,nprocs,npopsize,inum,length,ndim1
      include 'mpif.h'
      character(50) process

c
c      o Integer  n  is the parameter space dimension, i.e., the number
c        of adjustable parameters. 
c
c      o Function  ff  is a user-supplied scalar function of n vari-
c        ables, which must have the calling sequence f = ff(n,x), where
c        x is a double precision parameter array of length n.  This function must
c        be written so as to bound all parameters to the interval [0,1];
c        that is, the user must determine a priori bounds for the para-
c        meter space, and ff must use these bounds to perform the appro-
c        priate scalings to recover true parameter values in the
c        a priori ranges.
c
c        By convention, ff should return higher values for more optimal
c        parameter values (i.e., individuals which are more "fit").
c        For example, in fitting a function through data points, ff
c        could return the inverse of chi**2.
c
c        In most cases initialization code will have to be written
c        (either in a driver or in a separate subroutine) which loads
c        in data values and communicates with ff via one or more labeled
c        common blocks.  An example exercise driver and fitness function
c        are provided in the accompanying file, xpkaia.f.
c
c
c      Input/Output:
       double precision ctrl(12)
c
c      o Array  ctrl  is an array of control flags and parameters, to
c        control the genetic behavior of the algorithm, and also printed
c        output.  A default value will be used for any control variable
c        which is supplied with a value less than zero.  On exit, ctrl
c        contains the actual values used as control variables.  The
c        elements of ctrl and their defaults are:
c
c           ctrl( 1) - number of individuals in a population (default
c                      is 100)
c           ctrl( 2) - number of generations over which solution is
c                      to evolve (default is 500)
c           ctrl( 3) - number of significant digits (i.e., number of
c                      genes) retained in chromosomal encoding (default
c                      is 6)  (Note: This number is limited by the
c                      machine floating point precision.  Most 32-bit
c                      floating point representations have only 6 full
c                      digits of precision.  To achieve greater preci-
c                      sion this routine could be converted to double
c                      precision, but note that this would also require
c                      a double precision random number generator, which
c                      likely would not have more than 9 digits of
c                      precision if it used 4-byte integers internally.)
c           ctrl( 4) - crossover probability; must be  <= 1.0 (default
c                      is 0.85)
c           ctrl( 5) - mutation mode; 1/2=steady/variable (default is 2)
c           ctrl( 6) - initial mutation rate; should be small (default
c                      is 0.005) (Note: the mutation rate is the proba-
c                      bility that any one gene locus will mutate in
c                      any one generation.)
c           ctrl( 7) - minimum mutation rate; must be >= 0.0 (default
c                      is 0.0005)
c           ctrl( 8) - maximum mutation rate; must be <= 1.0 (default
c                      is 0.25)
c           ctrl( 9) - relative fitness differential; range from 0
c                      (none) to 1 (maximum).  (default is 1.)
c           ctrl(10) - reproduction plan; 1/2/3=Full generational
c                      replacement/Steady-state-replace-random/Steady-
c                      state-replace-worst (default is 3)
c           ctrl(11) - elitism flag; 0/1=off/on (default is 0)
c                      (Applies only to reproduction plans 1 and 2)
c           ctrl(12) - printed output 0/1/2=None/Minimal/Verbose
c                      (default is 0)
c
c
c     Output:
      double precision      x(n), f
      integer   status
c
c      o Array  x(1:n)  is the "fittest" (optimal) solution found,
c         i.e., the solution which maximizes fitness function ff
c
c      o Scalar  f  is the value of the fitness function at x
c
c      o Integer  status  is an indicator of the success or failure
c         of the call to pikaia (0=success; non-zero=failure)
c
c
c     Constants
      integer   NMAX, PMAX, DMAX
      parameter (NMAX = 60, PMAX = 1280, DMAX = 10)
c
c      o NMAX is the maximum number of adjustable parameters
c        (n <= NMAX)
c
c      o PMAX is the maximum population (ctrl(1) <= PMAX)
c
c      o DMAX is the maximum number of Genes (digits) per Chromosome
c        segement (parameter) (ctrl(3) <= DMAX)
c
c
c     Local variables

      integer        np, nd, ngen, imut, irep, ielite, ivrb, k, ip, ig,
     +               ip1, ip2, new, newtot,ndim,nparams
      double precision           pcross, pmut, pmutmn, pmutmx, fdif
c
      double precision           ph(NMAX,2), oldph(NMAX,PMAX), 
     $   newph(NMAX,PMAX),tmp(n,npopsize),fitmp(npopsize)
c
      integer        gn1(NMAX*DMAX), gn2(NMAX*DMAX)
      integer        ifit(PMAX), jfit(PMAX)
      double precision           fitns(PMAX)
c
c     User-supplied uniform random number generator
      double precision           urand
      external       urand
c
c     Function urand should not take any arguments.  If the user wishes
c     to be able to initialize urand, so that the same sequence of
c     random numbers can be repeated, this capability could be imple-
c     mented with a separate subroutine, and called from the user's
c     driver program.  An example urand function (and initialization
c     subroutine) which uses the function ran0 (the "minimal standard"
c     random number generator of Park and Miller [Comm. ACM 31, 1192-
c     1201, Oct 1988; Comm. ACM 36 No. 7, 105-110, July 1993]) is
c     provided.
c
c
c     Set control variables from input and defaults






      call setctl
     +   (ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut,
     +    fdif,irep,ielite,ivrb,status)
      if (status .ne. 0) then
         write(*,*) ' Control vector (ctrl) argument(s) invalid'
         return
      endif
 
c     Make sure locally-dimensioned arrays are big enough
      if (n.gt.NMAX .or. np.gt.PMAX .or. nd.gt.DMAX) then
         write(*,*)
     +      ' Number of parameters, population, or genes too large'
         status = -1
         return
      endif
 
c     Compute initial (random but bounded) phenotypes
      if(myrank.eq.0)then
         do 1 ip=1,np
             do 2 k=1,n
               oldph(k,ip)=urand()
               tmp(k,ip)=oldph(k,ip)
 2            continue
c       fitns(ip) = ff(n,oldph(1,ip))
    1 continue
      end if

c      if(myrank.eq.0)call matprt(tmp,n,np,n,np)
cparallel evaluate fitness of population with several processors

c        call mpi_barrier(mpi_comm_world,ierr)
        call passout(tmp,fitmp)
        if(myrank.eq.0)then
        call dcopy(npopsize,fitmp,1,fitns,1)
        end if



 
c     Rank initial population by fitness order
      call rnkpop(np,fitns,ifit,jfit)
 
c     Main Generation Loop
      do 10 ig=1,ngen

         print*,'generation',ig
          if(myrank.eq.0)then
          call system('date')
          open(unit=24,file='genwrite')
          write(24,*)ig
          close(24)
          end if

          call mpi_barrier(mpi_comm_world,ierr)

c        Main Population Loop
c         if(myrank.eq.0)then
 
         newtot=0
         do 20 ip=1,np/2
            
     
c           1. pick two parents
            call select(np,jfit,fdif,ip1)
   21       call select(np,jfit,fdif,ip2)
            if (ip1.eq.ip2) goto 21
 
c           2. encode parent phenotypes
            call encode(n,nd,oldph(1,ip1),gn1)
            call encode(n,nd,oldph(1,ip2),gn2)
 
c           3. breed
            call cross(n,nd,pcross,gn1,gn2)
            call mutate(n,nd,pmut,gn1)
            call mutate(n,nd,pmut,gn2)
 
c           4. decode offspring genotypes
            call decode(n,nd,gn1,ph(1,1))
            call decode(n,nd,gn2,ph(1,2))
 
c           5. insert into population
            if (irep.eq.1) then
               call genrep(NMAX,n,np,ip,ph,newph)
            else
               call stdrep(ff,NMAX,n,np,irep,ielite,
     +                     ph,oldph,fitns,ifit,jfit,new)
               newtot = newtot+new
            endif
 
c        End of Main Population Loop
   20    continue
c         end if

         call mpi_barrier(mpi_comm_world,ierr)
c         open(unit=13,file='generation')
c      write(13,*)'generation # ',ig,' out of',ngen,' total'
c      close(13)
 
c        if running full generational replacement: swap populations
         if (irep.eq.1)
     +      call newpop(ff,ielite,NMAX,n,np,oldph,newph,
     +         ifit,jfit,fitns,newtot)
 
c        adjust mutation rate?
         if (imut.eq.2) call adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
c
c        print generation report to standard output?
         if (ivrb.gt.0) call report
     +      (ivrb,NMAX,n,np,nd,oldph,fitns,ifit,pmut,ig,newtot)
 
c     End of Main Generation Loop
   10 continue
c
c     Return best phenotype and its fitness
      do 30 k=1,n
         x(k) = oldph(k,ifit(np))
   30 continue
      f = fitns(ifit(np))
c


      end
c********************************************************************
      subroutine setctl
     +   (ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut,
     +    fdif,irep,ielite,ivrb,status)
c===================================================================
c     Set control variables and flags from input and defaults
c===================================================================
      implicit none
c
c     Input
      integer  n
c
c     Input/Output
      double precision     ctrl(12)
c
c     Output
      integer  np, ngen, nd, imut, irep, ielite, ivrb, status
      double precision     pcross, pmutmn, pmutmx, pmut, fdif
c
c     Local
      integer  i
      double precision     DFAULT(12)
      save     DFAULT
      data     DFAULT /100,500,5,.85,2,.005,.0005,.25,1,1,1,0/
c
      do 1 i=1,12
         if (ctrl(i).lt.0.) ctrl(i)=DFAULT(i)
    1 continue
 
      np = ctrl(1)
      ngen = ctrl(2)
      nd = ctrl(3)
      pcross = ctrl(4)
      imut = ctrl(5)
      pmut = ctrl(6)
      pmutmn = ctrl(7)
      pmutmx = ctrl(8)
      fdif = ctrl(9)
      irep = ctrl(10)
      ielite = ctrl(11)
      ivrb = ctrl(12)
      status = 0
c
c     Print a header
      if (ivrb.gt.0) then
 
         write(*,2) ngen,np,n,nd,pcross,pmut,pmutmn,pmutmx,fdif
    2    format(/1x,60('*'),/,
     +      ' *',13x,'PIKAIA Genetic Algorithm Report ',13x,'*',/,
     +           1x,60('*'),//,
     +      '   Number of Generations evolving: ',i4,/,
     +      '       Individuals per generation: ',i4,/,
     +      '    Number of Chromosome segments: ',i4,/,
     +      '    Length of Chromosome segments: ',i4,/,
     +      '            Crossover probability: ',f9.4,/,
     +      '            Initial mutation rate: ',f9.4,/,
     +      '            Minimum mutation rate: ',f9.4,/,
     +      '            Maximum mutation rate: ',f9.4,/,
     +      '    Relative fitness differential: ',f9.4)
         if (imut.eq.1) write(*,3) 'Constant'
         if (imut.eq.2) write(*,3) 'Variable'
    3    format(
     +      '                    Mutation Mode: ',A)
         if (irep.eq.1) write(*,4) 'Full generational replacement'
         if (irep.eq.2) write(*,4) 'Steady-state-replace-random'
         if (irep.eq.3) write(*,4) 'Steady-state-replace-worst'
    4    format(
     +      '                Reproduction Plan: ',A)
      endif
 
c     Check some control values
      if (imut.ne.1 .and. imut.ne.2) then
         write(*,10)
         status = 5
      endif
   10 format(' ERROR: illegal value for imut (ctrl(5))')
 
      if (fdif.gt.1.) then
         write(*,11)
         status = 9
      endif
   11 format(' ERROR: illegal value for fdif (ctrl(9))')
 
      if (irep.ne.1 .and. irep.ne.2 .and. irep.ne.3) then
         write(*,12)
         status = 10
      endif
   12 format(' ERROR: illegal value for irep (ctrl(10))')
 
      if (pcross.gt.1.0 .or. pcross.lt.0.) then
         write(*,13)
         status = 4
      endif
   13 format(' ERROR: illegal value for pcross (ctrl(4))')
 
      if (ielite.ne.0 .and. ielite.ne.1) then
         write(*,14)
         status = 11
      endif
   14 format(' ERROR: illegal value for ielite (ctrl(11))')
 
      if (irep.eq.1 .and. imut.eq.1 .and. pmut.gt.0.5 .and.
     +    ielite.eq.0) then
         write(*,15)
      endif
   15 format(' WARNING: dangerously high value for pmut (ctrl(6));',
     +      /' (Should enforce elitism with ctrl(11)=1.)')
 
      if (irep.eq.1 .and. imut.eq.2 .and. pmutmx.gt.0.5 .and.
     +    ielite.eq.0) then
         write(*,16)
      endif
   16 format(' WARNING: dangerously high value for pmutmx (ctrl(8));',
     +      /' (Should enforce elitism with ctrl(11)=1.)')
 
      if (fdif.lt.0.33 .and. irep.ne.3) then
         write(*,17)
      endif
   17 format(' WARNING: dangerously low value of fdif (ctrl(9))')
 
      if (mod(np,2).gt.0) then
         np=np-1
         write(*,18) np
      endif
   18 format(' WARNING: decreasing population size (ctrl(1)) to np=',i4)
 
      return
      end
c********************************************************************
      subroutine report
     +   (ivrb,ndim,n,np,nd,oldph,fitns,ifit,pmut,ig,nnew)
c
c     Write generation report to standard output
c
      implicit none
 
c     Input:
      integer np,ifit(np),ivrb,ndim,n,nd,ig,nnew
      double precision oldph(ndim,np),fitns(np),pmut
c
c     Output: none
c
c     Local
      double precision bestft,pmutpv
      save bestft,pmutpv
      integer ndpwr,k
      logical rpt
      data bestft,pmutpv /0,0/
c
      rpt=.false.
 
      if (pmut.ne.pmutpv) then
         pmutpv=pmut
         rpt=.true.
      endif
 
      if (fitns(ifit(np)).ne.bestft) then
         bestft=fitns(ifit(np))
         rpt=.true.
      endif
 
      if (rpt .or. ivrb.ge.2) then
 
c        Power of 10 to make integer genotypes for display
         ndpwr = nint(10.**nd)
 
         write(*,'(/i6,i6,f10.6,4f10.6)') ig,nnew,pmut,
     +      fitns(ifit(np)), fitns(ifit(np-1)), fitns(ifit(np/2))
         do 15 k=1,n
            write(*,'(22x,3i10)')
     +         nint(ndpwr*oldph(k,ifit(np  ))),
     +         nint(ndpwr*oldph(k,ifit(np-1))),
     +         nint(ndpwr*oldph(k,ifit(np/2)))
   15    continue
 
      endif
      end

c**********************************************************************
c                         GENETICS MODULE
c**********************************************************************
c
c     ENCODE:    encodes phenotype into genotype
c                called by: PIKAIA
c
c     DECODE:    decodes genotype into phenotype
c                called by: PIKAIA
c
c     CROSS:     Breeds two offspring from two parents
c                called by: PIKAIA
c
c     MUTATE:    Introduces random mutation in a genotype
c                called by: PIKAIA
c
c     ADJMUT:    Implements variable mutation rate
c                called by: PIKAIA
c
c**********************************************************************
      subroutine encode(n,nd,ph,gn)
c======================================================================
c     encode phenotype parameters into integer genotype
c     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
c======================================================================
c
      implicit none
c
c     Inputs:
      integer   n, nd
      double precision      ph(n)
c
c     Output:
      integer   gn(n*nd)
c
c     Local:
      integer   ip, i, j, ii
      double precision      z
c
      z=10.**nd
      ii=0
      do 1 i=1,n
         ip=int(ph(i)*z)
         do 2 j=nd,1,-1
            gn(ii+j)=mod(ip,10)
            ip=ip/10
    2   continue
        ii=ii+nd
    1 continue
 
      return
      end
 
c**********************************************************************
      subroutine decode(n,nd,gn,ph)
c======================================================================
c     decode genotype into phenotype parameters
c     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
c======================================================================
c
      implicit none
c
c     Inputs:
      integer   n, nd, gn(n*nd)
c
c     Output:
      double precision      ph(n)
c
c     Local:
      integer   ip, i, j, ii
      double precision      z
c
      z=10.**(-nd)
      ii=0
      do 1 i=1,n
         ip=0
         do 2 j=1,nd
            ip=10*ip+gn(ii+j)
    2    continue
         ph(i)=ip*z
         ii=ii+nd
    1 continue
 
      return
      end
 
c**********************************************************************
      subroutine cross(n,nd,pcross,gn1,gn2)
c======================================================================
c     breeds two parent chromosomes into two offspring chromosomes
c     breeding occurs through crossover starting at position ispl
c======================================================================
c     USES: urand
      implicit none
c
c     Inputs:
      integer        n, nd
      double precision           pcross
c
c     Input/Output:
      integer        gn1(n*nd), gn2(n*nd)
c
c     Local:
      integer        i, ispl, t
c
c     Function
      double precision           urand
      external       urand
 
 
c     Use crossover probability to decide whether a crossover occurs
      if (urand().lt.pcross) then
 
c        Compute crossover point
         ispl=int(urand()*n*nd)+1
 
c        Swap genes at ispl and above
         do 10 i=ispl,n*nd
            t=gn2(i)
            gn2(i)=gn1(i)
            gn1(i)=t
   10    continue
      endif
 
      return
      end
 
c**********************************************************************
      subroutine mutate(n,nd,pmut,gn)
c======================================================================
c     Mutations occur at rate pmut at all gene loci
c======================================================================
c     USES: urand
      implicit none
c
c     Input:
      integer        n, nd
      double precision           pmut
c
c     Input/Output:
      integer        gn(n*nd)
c
c     Local:
      integer        i
c
c     Function:
      double precision           urand
      external       urand
c
c     Subject each locus to mutation at the rate pmut
      do 10 i=1,n*nd
         if (urand().lt.pmut) then
            gn(i)=int(urand()*10.)
         endif
   10 continue
 
      return
      end
 
c**********************************************************************
      subroutine adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
c======================================================================
c     dynamical adjustment of mutation rate; criterion is relative
c     difference in absolute fitnesses of best and median individuals
c======================================================================
c
      implicit none
c
c     Input:
      integer        np, ifit(np)
      double precision           fitns(np), pmutmn, pmutmx
c
c     Input/Output:
      double precision           pmut
c
c     Local:
      double precision           rdif, rdiflo, rdifhi, delta
      parameter      (rdiflo=0.05, rdifhi=0.25, delta=1.5)
 
      rdif=abs(fitns(ifit(np))-fitns(ifit(np/2)))/
     +        (fitns(ifit(np))+fitns(ifit(np/2)))
      if(rdif.le.rdiflo)then
         pmut=min(pmutmx,pmut*delta)
      else if(rdif.ge.rdifhi)then
         pmut=max(pmutmn,pmut/delta)
      endif
 
      return
      end

c**********************************************************************
c                       REPRODUCTION MODULE
c**********************************************************************
c
c     SELECT:   Parent selection by roulette wheel algorithm
c               called by: PIKAIA
c
c     RNKPOP:   Ranks initial population
c               called by: PIKAIA, NEWPOP
c
c     GENREP:   Inserts offspring into population, for full
c               generational replacement
c               called by: PIKAIA
c
c     STDREP:   Inserts offspring into population, for steady-state
c               reproduction
c               called by: PIKAIA
c               calls:     FF
c
c     NEWPOP:   Replaces old generation with new generation
c               called by: PIKAIA
c               calls:     FF, RNKPOP
c
c**********************************************************************
      subroutine select(np,jfit,fdif,idad)
c======================================================================
c     Selects a parent from the population, using roulette wheel
c     algorithm with the relative fitnesses of the phenotypes as
c     the "hit" probabilities [see Davis 1991, chap. 1].
c======================================================================
c     USES: urand
      implicit none
c
c     Input:
      integer        np, jfit(np)
      double precision           fdif
c
c     Output:
      integer        idad
c
c     Local:
      integer        np1, i
      double precision           dice, rtfit
c
c     Function:
      double precision           urand
      external       urand
c
c
      np1 = np+1
      dice = urand()*np*np1
      rtfit = 0.
      do 1 i=1,np
         rtfit = rtfit+np1+fdif*(np1-2*jfit(i))
         if (rtfit.ge.dice) then
            idad=i
            goto 2
         endif
    1 continue
c     Assert: loop will never exit by falling through
 
    2 return
      end
 
c**********************************************************************
      subroutine rnkpop(n,arrin,indx,rank)
c======================================================================
c     Calls external sort routine to produce key index and rank order
c     of input array arrin (which is not altered).
c======================================================================
c     USES: rqsort
      implicit none
c
c     Input
      integer    n
      double precision       arrin(n)
c
c     Output
      integer    indx(n),rank(n)
c
c     Local
      integer    i
c
c     External sort subroutine
      external rqsort
c
c
c     Compute the key index
      call rqsort(n,arrin,indx)
c
c     ...and the rank order
      do 1 i=1,n
         rank(indx(i)) = n-i+1
    1 continue
      return
      end
 
c***********************************************************************
      subroutine genrep(ndim,n,np,ip,ph,newph)
c=======================================================================
c     full generational replacement: accumulate offspring into new
c     population array
c=======================================================================
c
      implicit none
 
c     Input:
      integer        ndim, n, np, ip
      double precision           ph(ndim,2)
c
c     Output:
      double precision           newph(ndim,np)
c
c     Local:
      integer        i1, i2, k
c
c
c     Insert one offspring pair into new population
      i1=2*ip-1
      i2=i1+1
      do 1 k=1,n
         newph(k,i1)=ph(k,1)
         newph(k,i2)=ph(k,2)
    1 continue
 
      return
      end
 
c**********************************************************************
      subroutine stdrep
     +   (ff,ndim,n,np,irep,ielite,ph,oldph,fitns,ifit,jfit,nnew)
c======================================================================
c     steady-state reproduction: insert offspring pair into population
c     only if they are fit enough (replace-random if irep=2 or
c     replace-worst if irep=3).
c======================================================================
c     USES: ff, urand
      implicit none
c
c     Input:
      integer        ndim, n, np, irep, ielite
      double precision           ff, ph(ndim,2)
      external       ff
c
c     Input/Output:
      double precision           oldph(ndim,np), fitns(np)
      integer        ifit(np), jfit(np)
c
c     Output:
      integer        nnew
 
c     Local:
      integer        i, j, k, i1, if1
      double precision           fit
c
c     External function
      double precision           urand
      external       urand
c
c
      nnew = 0
      do 1 j=1,2
 
c        1. compute offspring fitness (with caller's fitness function)
          fit=ff(n,ph(1,j))
 
c        2. if fit enough, insert in population
         do 20 i=np,1,-1
            if (fit.gt.fitns(ifit(i))) then
 
c              make sure the phenotype is not already in the population
               if (i.lt.np) then
                  do 5 k=1,n
                     if (oldph(k,ifit(i+1)).ne.ph(k,j)) goto 6
    5             continue
                  goto 1
    6             continue
               endif
 
c              offspring is fit enough for insertion, and is unique
 
c              (i) insert phenotype at appropriate place in population
               if (irep.eq.3) then
                  i1=1
               else if (ielite.eq.0 .or. i.eq.np) then
                  i1=int(urand()*np)+1
               else
                  i1=int(urand()*(np-1))+1
               endif
               if1 = ifit(i1)
               fitns(if1)=fit
               do 21 k=1,n
                  oldph(k,if1)=ph(k,j)
   21          continue
 
c              (ii) shift and update ranking arrays
               if (i.lt.i1) then
 
c                 shift up
                  jfit(if1)=np-i
                  do 22 k=i1-1,i+1,-1
                     jfit(ifit(k))=jfit(ifit(k))-1
                     ifit(k+1)=ifit(k)
   22             continue
                  ifit(i+1)=if1
               else
 
c                 shift down
                  jfit(if1)=np-i+1
                  do 23 k=i1+1,i
                     jfit(ifit(k))=jfit(ifit(k))+1
                     ifit(k-1)=ifit(k)
   23             continue
                  ifit(i)=if1
               endif
               nnew = nnew+1
               goto 1
            endif
   20    continue
 
    1 continue
 
      return
      end
 
c**********************************************************************
      subroutine newpop
     +   (ff,ielite,ndim,n,np,oldph,newph,ifit,jfit,fitns,nnew)
c======================================================================
c     replaces old population by new; recomputes fitnesses & ranks
c======================================================================
c     USES: ff, rnkpop
      implicit none
c
c     Input:
      integer        ndim, np, n, ielite
      double precision           ff
      external       ff
c
c     Input/Output:
      double precision           oldph(ndim,np), newph(ndim,np)
c
c     Output:
      common /parallel/ndim1,npopsize,nparams,myrank,rankno
      character(4) rankno(10000)
      integer        ifit(np), jfit(np), nnew,npopsize,nparams,myrank
     $,ndim1
      double precision      fitns(np),fitmp(np),tmp(nparams,npopsize)
c
c     Local:
      integer        i, k
c
c
      nnew = np
 
c     if using elitism, introduce in new population fittest of old
c     population (if greater than fitness of the individual it is
c     to replace)
      if (ielite.eq.1 .and. ff(n,newph(1,1)).lt.fitns(ifit(np))) then
         do 1 k=1,n
            newph(k,1)=oldph(k,ifit(np))
    1    continue
         nnew = nnew-1
      endif
 
c     replace population
      do 2 i=1,np
c         open(unit=14,file='individual')
c         write(14,*)'working on person number',i
c         close(14)
         do 3 k=1,n
            oldph(k,i)=newph(k,i)
            tmp(k,i)=newph(k,i)
    3    continue
 
c        get fitness using caller's fitness function
c        fitns(i)=ff(n,oldph(1,i))
    2 continue

      
        call passout(tmp,fitmp)
        if(myrank.eq.0)then
        call dcopy(np,fitmp,1,fitns,1)
        end if

 
c     compute new population fitness rank order
      call rnkpop(np,fitns,ifit,jfit)
 
      return
      end

c*********************************************************************
      subroutine twod2(n,x,value,y,ndata,septs,ab)
      implicit double precision (a-h,o-z)
      common /parallel/ndim1,npopsize,nparams,myrank,rankno
      integer n,i,npts
      double precision x(n),y(n),ab(ndata),septs(ndata),s(n)
      double precision twod,factor
      character(50) word
      character(4) rankno(10000)
      dimension namea(1000),nameb(1000),axyz(3,1000),bxyz(3,1000)
     $     ,xyz_rot(3,2000),nameab(2000)
      integer freeze(6),na,nb
      character*10 cjunk
      character(5) ctrash,ch,ch3
      data itimes /0/

      value=0.0d0
      itimes=itimes+1

           icount=0
      open(unit=10,file='input')
      do 100 i=1,6
         icount=icount+1
c         print*,'reading'
         read(10,*)s(i),freeze(icount),cjunk
c         print*,x(icount),freeze(icount),cjunk
 100                      continue
      read(10,*)force,cjunk
      read(10,*)na,cjunk
      do i=1,na
      read(10,*)namea(i),axyz(1,i),axyz(2,i),axyz(3,i)
c      print*,na,namea(i),axyz(1,i),axyz(2,i),axyz(3,i)
      end do
      read(10,*)nb,cjunk
      do i=1,nb
      read(10,*)nameb(i),bxyz(1,i),bxyz(2,i),bxyz(3,i)
      end do
      close(10)

c      rfac=0.9
c       do i=1,n
c       y(i)=s(i)+(2.*rfac*(x(i)-0.5d0)*s(i))
c       end do

      y(1)= 20.0d0*x(1)
      y(2)= 180.0d0*x(2)
      y(3)= 360.0d0*x(3)
      y(4)= 360.0d0*x(4)
      y(5)= 360.0d0*x(5)
      y(6)= 360.0d0*x(6)





        word=''
       word(1:4)='xyz.'
       word(5:)=rankno(myrank+1)
       open(unit=1,file=word)
       rewind 1
        write(1,30)y(1),y(2),y(3),y(4),y(5),y(6)
        write(1,*)na
        do i=1,na
        write(1,*)namea(i),axyz(1,i),axyz(2,i),axyz(3,i)
        end do
        write(1,*)nb
        do i=1,nb
        write(1,*)nameb(i),bxyz(1,i),bxyz(2,i),bxyz(3,i)
        end do
        close(1)
c        call system('./gen.exe')

 30                    format(f15.10,f15.10,f15.10,f15.10,f15.10,f15.10)

          word=''
          word(1:11)='./xrunprog '
          word(12:)=rankno(myrank+1)
          call system(word)
          print*,myrank

          return


 3                       format(I4,f20.5,f20.5,f20.5)


         word=''
         word(1:5)='gout.'
         word(6:)=rankno(myrank+1)
         open(unit=1,file=word)
         do i=1,1000000
         read(1,*,iostat=io)ch3
          if(io<0)then
          value=1.0d-10
          print*,'Gaussian did not run'
          close(1)
          return
          end if
          if(ch3 =="SCF")then
             BACKSPACE (unit=1)
             exit
          end if
       end do
       read(1,*)ch3,ch3,ch3,ch3,ener,ch3,ch3,isteps,ch3
       close(1)
       if(isteps .gt. 50)then
           value=1.0d-10
          print*,'SCF did not converge'
          return
       end if


         word=''
         word(1:6)='gouta.'
         word(7:)=rankno(myrank+1)
         open(unit=1,file=word)
         do i=1,1000000
         read(1,*,iostat=io)ch3
          if(io<0)then
          value=1.0d-10
          print*,'Gaussian did not run'
          close(1)
          return
          end if
          if(ch3 =="SCF")then
             BACKSPACE (unit=1)
             exit
          end if
       end do
       read(1,*)ch3,ch3,ch3,ch3,enera,ch3,ch3,isteps,ch3
       close(1)
       if(isteps .gt. 50)then
           value=1.0d-10
          print*,'SCF did not converge'
          return
       end if


        word=''
         word(1:6)='goutb.'
         word(7:)=rankno(myrank+1)
         open(unit=1,file=word)
         do i=1,1000000
         read(1,*,iostat=io)ch3
          if(io<0)then
          value=1.0d-10
          print*,'Gaussian did not run'
          close(1)
          return
          end if
          if(ch3 =="SCF")then
             BACKSPACE (unit=1)
             exit
          end if
       end do
       read(1,*)ch3,ch3,ch3,ch3,enerb,ch3,ch3,isteps,ch3
       close(1)
       if(isteps .gt. 50)then
           value=1.0d-10
          print*,'SCF did not converge'
          return
       end if








c                    value=627.51*ener/27.2113838565563d0

c              value2=value- (-71.76252768506d0*627.51) 

           value=ener-(enera+enerb)

                    value=627.51*value

              value2=value



c        print*,'value',value,value2
c        call evalf(value,y,n,istop,value2)

c        twod=-value2
        print*,'value on node ',myrank, 'is: ',value2
c        print*,myrank,value2,value,268.4823039*627.51
        value=value2

      RETURN
      END







      function twod(n,x)
      twod=0
      return
      end














