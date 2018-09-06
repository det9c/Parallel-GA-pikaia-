              subroutine evalf(f,Y,n,ab,septs,npts,istop,ilow)
      implicit double precision (a-h,o-z)
      common /parallel/ndim1,npopsize,nparams,myrank,rankno
      character(3) rankno(128)
      character(50) word
       include 'mpif.h'

      DOUBLE PRECISION f,fold,Y(n),ab(npts),septs(npts),Y2(108)
      integer npts
            data fold /1.0D20/
      if(f.lt.fold)then
      open(unit=11,file='lowest-so-far')
      rewind 11
      write(11,*)'lowest f so far is',f,y
      close(11)

      fold=f

               open(unit=24,file='genwrite')
         read(24,*)ig
         close(24)



        open(unit=12,file='history',access='append')
        write(12,*)f,ig
        close(12)




!      call system('cp INPUT INPUT.lowest')
!      call system('cp out out.lowest')
!      call system('cp md_history.xyz md_history.lowest.xyz')
      
         word(1:8)='xgetlow '
         word(9:)=rankno(ilow)

       call system(word)
       print*,'lowest is in direcoty',rankno(ilow)



      end if
 
c      if(f.lt.1.00D-3)then
c      print*,'stoppin cuz f is small enough in evalf',f
c      istop=1
c              OPEN(UNIT=40,FILE='stop')
c      WRITE(40,*)2
c      CLOSE(40)
c      end if

      return
      end


