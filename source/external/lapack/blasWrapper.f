c
c *******************************************
c
      subroutine wrapdgemv(trans,m,n,alpha,A,ldA,x,incx,beta,y,incy)
      logical*1 trans
      character tr
      external dgemv
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      call dgemv(tr,m,n,alpha,A,ldA,x,incx,beta,y,incy)
      end
c
c *******************************************
c
      subroutine wrapdsymv(upper,n,alpha,A,ldA,x,incx,beta,y,incy)
      logical*1 upper
      character uplo
      external dsymv
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dsymv(uplo,n,alpha,A,ldA,x,incx,beta,y,incy)
      end
c
c *******************************************
c
      subroutine wrapdtrmv(upper,trans,diag,n,A,ldA,x,incx)
      external dtrmv
      logical*1 upper, trans, diag
      character uplo,  tr,    diago
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      if(.not.diag) then
        diago = 'N'
      else
        diago = 'U'
      endif
      call dtrmv(uplo,tr,diago,n,A,ldA,x,incx)
      end
c
c *******************************************
c
      subroutine wrapdtrsv(upper,trans,diag,n,A,ldA,x,incx)
      external dtrsv
      logical*1 upper, trans, diag
      character uplo,  tr,    diago
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      if(.not.diag) then
        diago = 'N'
      else
        diago = 'U'
      endif
      call dtrsv(uplo,tr,diago,n,A,ldA,x,incx);
      end
c
c *******************************************
c
      subroutine wrapdgemm(transA,transB,m,n,k,alpha,
     $                     A,ldA,B,ldB,beta,C,ldC)
      external dgemm
      logical*1 transA, transB
      character trA, trB
      if(.not.transA) then
        trA = 'N'
      else
        trA = 'T'
      endif
      if(.not.transB) then
        trB = 'N'
      else
        trB = 'T'
      endif
      call dgemm(trA,trB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC)
      end
c
c *******************************************
c
      subroutine wrapdsymm(left,upper,m,n,alpha,A,ldA,B,ldB,beta,C,ldC)
      external dsymm
      logical*1 left, upper
      character side, uplo
      if(.not.left) then
        side = 'R'
      else
        side = 'L'
      endif
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      call dsymm(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC)
      end
c
c *******************************************
c
      subroutine wrapdtrmm(left,upper,trans,diag,m,n,alpha,A,ldA,B,ldB)
      external dtrmm
      logical*1 left, upper, trans, diag
      character side, uplo,  tr,    diago
      if(.not.left) then
        side = 'R'
      else
        side = 'L'
      endif
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      if(.not.diag) then
        diago = 'N'
      else
        diago = 'U'
      endif
      call dtrmm(side,uplo,tr,diago,m,n,alpha,A,ldA,B,ldB)
      end
c
c *******************************************
c
      subroutine wrapdtrsm(left,upper,trans,diag,m,n,alpha,A,ldA,B,ldB)
      external dtrsm
      logical*1 left, upper, trans, diag
      character side, uplo,  tr,    diago
      if(.not.left) then
        side = 'R'
      else
        side = 'L'
      endif
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      if(.not.diag) then
        diago = 'N'
      else
        diago = 'U'
      endif
      call dtrsm(side,uplo,tr,diago,m,n,alpha,A,ldA,B,ldB)
      end
c
c *******************************************
c
      subroutine wrapdsyrk(upper,trans,n,k,alpha,A,ldA,beta,C,ldC)
      external dsyrk
      logical*1 upper, trans
      character uplo, tr
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      call dsyrk(uplo,tr,n,k,alpha,A,ldA,beta,C,ldC)
      end
c
c *******************************************
c
      subroutine wrapdsyr2k(upper,trans,n,k,alpha,A,ldA,
     $                      B,ldB,beta,C,ldC)
      external dsyr2k
      logical*1 upper, trans
      character uplo, tr
      if(.not.upper) then
        uplo = 'L'
      else
        uplo = 'U'
      endif
      if(.not.trans) then
        tr = 'N'
      else
        tr = 'T'
      endif
      call dsyr2k(uplo,tr,n,k,alpha,A,ldA,B,ldB,beta,C,ldC)
      end
c
c *******************************************
c
