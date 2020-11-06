c
c *******************************************
c
      subroutine wrapdlacpy(m,n,A,ldA,B,ldB)
      external dlacpy
      call dlacpy('G',m,n,A,ldA,B,ldB)
      end
c
c *******************************************
c
      subroutine wrapdpotrf(upper,n,A,ldA,info)
      integer   upper
      character uplo
      external dpotrf
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dpotrf(uplo,n,A,ldA,info)
      end
c
c *******************************************
c
      subroutine wrapdpstrf(upper,n,A,ldA,piv,rank,
     $                      tol,work,info)
      integer   upper
      character uplo
      external dpstrf
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dpstrf(uplo,n,A,ldA,piv,rank,tol,work,info)
      end
c
c *******************************************
c
      subroutine wrapdpotri(upper,n,A,ldA,info)
      integer   upper
      character uplo
      external dpotri
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dpotri(uplo,n,A,ldA,info)
      end
c
c *******************************************
c
      subroutine wrapdtrtri(upper,n,A,ldA,info)
      integer   upper
      character uplo
      external dpotri
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dtrtri(uplo,'N',n,A,ldA,info)
      end
c
c *******************************************
c
      subroutine wrapdlauum(upper,n,A,ldA,info)
      integer   upper
      character uplo
      external dlauum
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dlauum(uplo,n,A,ldA,info)
      end
c
c *******************************************
c
      subroutine wrapdgetrf(m,n,A,ldA,ipiv,info)
      external dgetrf
      call dgetrf(m,n,A,ldA,ipiv,info);
      end
c
c *******************************************
c
      subroutine wrapdgesv(n,nrhs,A,ldA,ipiv,B,ldB,info)
      external dgesv
      call dgesv(n,nrhs,A,ldA,ipiv,B,ldB,info);
      end
c
c *******************************************
c
      subroutine wrapdsgesv(n,nrhs,A,ldA,ipiv,B,ldB,X,ldX,
     $                      work,swork,iter,info)
      external dsgesv
      call dsgesv(n,nrhs,A,ldA,ipiv,B,ldB,X,ldX,
     $            work,swork,iter,info);
      end
c
c *******************************************
c
      subroutine wrapdgetri(n,A,ldA,ipiv,work,lwork,info)
      external dgetri
      call dgetri(n,A,ldA,ipiv,work,lwork,info);
      end
c
c *******************************************
c
      subroutine wrapdgels(trans,m,n,nrhs,A,ldA,B,ldB,work,lwork,info)
      integer   trans
      character tr
      external dgels
      if(trans.eq.0) then
        tr = 'N'
      else
        tr = 'T'
      endif
      call dgels(tr,m,n,nrhs,A,ldA,B,ldB,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdgeqrf(m,n,A,ldA,tau,work,lwork,info)
      external dgeqrf
      call dgeqrf(m,n,A,ldA,tau,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdormqr(left,trans,m,n,k,A,ldA,tau,
     $                      C,ldC,work,lwork,info)
      external dormqr
      integer   left, trans
      character side,  tr
      if(left.eq.0) then
        side = 'R'
      else
        side = 'L'
      endif
      if(trans.eq.0) then
        tr = 'N'
      else
        tr = 'T'
      endif
      call dormqr(side,tr,m,n,k,A,ldA,tau,C,ldC,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdorgqr(m,n,k,A,ldA,tau,work,lwork,info)
      external dgeqrf
      call dorgqr(m,n,k,A,ldA,tau,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdpbtrf(upper,n,kd,A,ldA,info)
      integer   upper
      character uplo
      external dpbtrf
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dpbtrf(uplo,n,kd,A,ldA,info)
      end
c
c *******************************************
c
      subroutine wrapdpbsv(upper,n,kd,nrhs,A,ldA,B,ldB,info)
      integer   upper
      character uplo
      external dpbtrf
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dpbsv(uplo,n,kd,nrhs,A,ldA,B,ldB,info)
      end
c
c *******************************************
c
      subroutine wrapdgbtrf(n,m,kl,ku,A,ldA,ipiv,info)
      external dgbtrf
      call dgbtrf(n,m,kl,ku,A,ldA,ipiv,info)
      end
c
c *******************************************
c
      subroutine wrapdtbtrs(upper,trans,diag,n,kd,nrhs,A,ldA,B,ldB,info)
      external dtrsm
      integer   upper, trans, diag
      character uplo,  tr,    diago
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(trans.eq.0) then
        tr = 'N'
      else
        tr = 'T'
      endif
      if(diag.eq.0) then
        diago = 'N'
      else
        diago = 'U'
      endif
      call dtbtrs(uplo,tr,diago,n,kd,nrhs,A,ldA,B,ldB,info)
      end
c
c *******************************************
c
      subroutine wrapdgbtrs(trans,n,kl,ku,nrhs,A,ldA,ipiv,B,ldB,info)
      external  dgbtrs
      integer   trans
      character tr
      if(trans.eq.0) then
        tr = 'N'
      else
        tr = 'T'
      endif
      call dgbtrs(tr,n,kl,ku,nrhs,A,ldA,ipiv,B,ldB,info)
      end
c
c *******************************************
c
      subroutine wrapdsyev(job,upper,n,A,ldA,W,work,lwork,info)
      external dsyev
      integer   job, upper
      character jobz,uplo
      if(upper.eq.0) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(job.eq.0) then
        jobz = 'N'
      else
        jobz = 'V'
      endif
      call dsyev(jobz,uplo,n,A,ldA,W,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdgeev(jobvl,jobvr,n,A,ldA,WR,WI,VL,ldVL,VR,
     $                     ldVR,work,lwork,info)
      external dgeev
      integer   jobvl, jobvr
      character jobvlc, jobvrc
      if(jobvl.eq.0) then
        jobvlc = 'N'
      else
        jobvlc = 'V'
      endif
      if(jobvr.eq.0) then
        jobvrc = 'N'
      else
        jobvrc = 'V'
      endif
      call dgeev(jobvlc,jobvrc,n,A,ldA,WR,WI,VL,ldVL,VR,
     $           ldVR,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdgesvd(jobz,jobvt,m,n,A,ldA,S,
     $                      U,ldU,Vt,ldVt,work,lwork,info)
      external dsyev
      integer   jobz,  jobvt
      character cjobz,cjobvt
      if(jobz.eq.0) then
        cjobz = 'N'
      else
        cjobz = 'S'
      endif
      if(jobvt.eq.0) then
        cjobvt = 'N'
      else
        cjobvt = 'S'
      endif
      call dgesvd(cjobz,cjobvt,m,n,A,ldA,
     $            S,U,ldU,Vt,ldVt,work,lwork,info)
      end
c
c *******************************************
c
      subroutine wrapdgesdd(jobz,m,n,A,ldA,S,
     $                      U,ldU,Vt,ldVt,work,lwork,iwork,info)
      external dsyev
      integer   jobz
      character cjobz
      if(jobz.eq.0) then
        cjobz = 'N'
      else
        cjobz = 'S'
      endif
      call dgesdd(cjobz,m,n,A,ldA,
     $            S,U,ldU,Vt,ldVt,work,lwork,iwork,info)
      end
c
c *******************************************
c
