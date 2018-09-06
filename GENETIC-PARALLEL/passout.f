      subroutine passout(xtrial,fitout)
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      parameter(ndata=1)
      common /parallel/ndim1,npopsize,nparams,myrank,rankno
      character(2) rankno(128)
      dimension xtrial(nparams,npopsize),xnode(nparams,ndim1)
     $,fitnode(ndim1),fitout(npopsize)
     $,ynode(nparams,ndim1),yscale(nparams,npopsize)
     $,ab(ndata),rlocpts(ndata,ndim1)
     $,allpts(ndata,npopsize)
c      include 'mpif.h'
      character(20) process

       call mpi_barrier(mpi_comm_world,ierr)
      nsend=nparams*npopsize

       call mpi_bcast(xtrial,nsend,mpi_double_precision
     $     ,0,mpi_comm_world,ierr)

      call mpi_barrier(mpi_comm_world,ierr)


      icol=myrank*ndim1+1
      nsend=nparams*ndim1
      call dcopy(nsend,xtrial(1,icol),1,xnode,1)
c      call matprt(xnode,nparams,ndim,nparams,ndim)


      do i=1,ndim1
      call twod2(nparams,xnode(1,i),value,ynode(1,i),ndata
     $,rlocpts(1,i),ab)
      fitnode(i)=-value
      end do

      call mpi_barrier(mpi_comm_world,ierr)
      

       call mpi_gather(fitnode,ndim1,mpi_double_precision
     $,fitout,ndim1
     $,mpi_double_precision,0,mpi_comm_world,ierr)

       nsend=nparams*ndim1
      call mpi_gather(ynode,nsend,mpi_double_precision
     $,yscale,nsend
     $,mpi_double_precision,0,mpi_comm_world,ierr)


      nsend=ndata*ndim1
      call mpi_gather(rlocpts,nsend,mpi_double_precision
     $,allpts,nsend
     $,mpi_double_precision,0,mpi_comm_world,ierr)



      ilow=1
      flow=-fitout(1)
       do i=2,npopsize
        ftmp=-fitout(i)
           if(ftmp.lt.flow)then
           flow=ftmp
           ilow=i
           end if
          end do

          istop=0
          if(myrank.eq.0)then
       call evalf(flow,yscale(1,ilow),nparams,ab,allpts(1,ilow)
     $,ndata,istop,ilow)

c       istop=1
       end if

       call mpi_barrier(mpi_comm_world,ierr)


        call mpi_bcast(istop,1,mpi_integer,0,mpi_comm_world,ierr)

        istop=1
         if(istop.eq.1)then
           call mpi_finalize(ierr)
	call system('date')
           stop
           end if


       
      return
      end


