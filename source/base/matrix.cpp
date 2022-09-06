/***********************************************/
/**
* @file matrix.cpp
*
* @brief matrix computations.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @author Sebastian Strasser
* @date 2001-06-16
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "external/lapack/blas.h"
#include "external/lapack/lapack.h"
#include "base/matrix.h"

/***********************************************/
/***** MatrixBase ******************************/
/***********************************************/

MatrixBase::MatrixBase(UInt size) : _size(size)
{
  try
  {
    ptr = std::shared_ptr<Double>(new Double[_size], std::default_delete<Double[]>());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** const_MatrixSlice ***********************/
/***********************************************/

const_MatrixSlice::const_MatrixSlice(UInt rows, UInt columns, Matrix::Type type, Matrix::Uplo uplo, Double fill) :
  _rows(rows),
  _columns(columns),
  _start(0),
  _ld(rows),
  _type(type),
  _uplo(uplo),
  _rowMajorOrder(FALSE)
{
  if(size())
  {
    base = std::make_shared<MatrixBase>(size());
    std::fill_n(base->field(), base->size(), fill);
  }
}

/***********************************************/

const_MatrixSlice const_MatrixSlice::slice(UInt startRow, UInt startColumn, UInt height, UInt width) const
{
  try
  {
    if((height==0) || (width==0))
      return Matrix(height, width);
    if((startRow+height>rows()) || (startColumn+width>columns()))
      throw(Exception("Dimension error: ("+rows()%"%i x "s+columns()%"%i).slice("s+startRow%"%i, "s+startColumn%"%i, "s+height%"%i, "s+width%"%i)"s));

    const_MatrixSlice x(*this);
    x._rows    = height;
    x._columns = width;
    x._start  += (x._rowMajorOrder) ? (startColumn + startRow*_ld) : (startRow + startColumn*_ld);
    if((startRow!=startColumn) || (height!=width))
      x._type = GENERAL;
    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** MatrixSlice *****************************/
/***********************************************/

MatrixSlice MatrixSlice::slice(UInt startRow, UInt startColumn, UInt height, UInt width) const
{
  try
  {
    if((height==0) || (width==0))
      return Matrix(height, width);
    if((startRow+height>rows()) || (startColumn+width>columns()))
      throw(Exception("Dimension error: ("+rows()%"%i x "s+columns()%"%i).slice("s+startRow%"%i, "s+startColumn%"%i, "s+height%"%i, "s+width%"%i)"s));

    MatrixSlice x(*this);
    x._rows    = height;
    x._columns = width;
    x._start  += (x._rowMajorOrder) ? (startColumn + startRow*_ld) : (startRow + startColumn*_ld);
    if((startRow!=startColumn) || (height!=width))
      x._type = GENERAL;
    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

const MatrixSlice &MatrixSlice::fill(Double f) const
{
  if(!size())
    return *this;
  if(_rowMajorOrder)
    for(UInt i=0; i<_rows; i++)
      std::fill_n(base->field()+(_start+i*_ld), _columns, f);
  else
    for(UInt i=0; i<_columns; i++)
      std::fill_n(base->field()+(_start+i*_ld), _rows, f);
  return *this;
}

/***********************************************/

MatrixSliceRef MatrixSlice::operator+= (const const_MatrixSlice &x) const
{
  try
  {
    axpy(1., x, *this);
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixSliceRef MatrixSlice::operator-= (const const_MatrixSlice &x) const
{
  try
  {
    axpy(-1., x, *this);
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixSliceRef MatrixSlice::operator += (Double c) const
{
  for(UInt i=0; i<rows(); i++)
    for(UInt j=0; j<columns(); j++)
      (*this)(i,j) += c;
  return *this;
}

/***********************************************/

MatrixSliceRef MatrixSlice::operator -= (Double c) const
{
  return (*this += -c);
}

/***********************************************/

MatrixSliceRef MatrixSlice::operator *= (Double c) const
{
  if(!size())
    return *this;
  const UInt rows    = (_rowMajorOrder) ? _columns : _rows;
  const UInt columns = (_rowMajorOrder) ? _rows : _columns;
  if(rows == _ld)
    blas_dscal(static_cast<F77Int>(rows*columns), c, field(), 1);
  else
    for(UInt i=0; i<columns; i++)
      blas_dscal(static_cast<F77Int>(rows), c, field()+i*_ld, 1);
  return *this;
}

/***********************************************/

MatrixSliceRef MatrixSlice::operator /= (Double c) const
{
  return (*this *= 1./c);
}

/***********************************************/
/***** Matrix **********************************/
/***********************************************/

void Matrix::assignment(const const_MatrixSlice &x)
{
  if(!size())
  {
    base = nullptr;
  }
  else if((_ld!=_rows) || _rowMajorOrder)
  {
    _start         = 0;
    _ld            = _rows;
    _rowMajorOrder = FALSE;
    base = std::make_shared<MatrixBase>(_rows*_columns);
    copy(x, *this);
  }
  else
  {
    base = std::make_shared<MatrixBase>(*x.base);
  }
}

/***********************************************/

Matrix::Matrix(const Matrix &x) : MatrixSlice(x)
{
  try
  {
    assignment(x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix::Matrix(const const_MatrixSlice &x) : MatrixSlice(x)
{
  try
  {
    assignment(x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix::Matrix(std::initializer_list<std::initializer_list<Double>> list) : MatrixSlice(list.size(), (list.size() ? list.begin()->size() : 0), GENERAL, UPPER)
{
  try
  {
    for(const auto &l : list)
      if(l.size() != columns())
        throw(Exception("All rows of matrix must have same number of columns when constructing from vectors. Expected "+columns()%"%i columns but got "s+l.size()%"%i."s));

    for(UInt i=0; i<rows(); i++)
      for(UInt k=0; k<columns(); k++)
        (*this)(i, k) = *((list.begin()+i)->begin()+k);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix &Matrix::operator=(const const_MatrixSlice &x)
{
  try
  {
    _rows          = x._rows;
    _columns       = x._columns;
    _start         = x._start;
    _ld            = x._ld;
    _type          = x._type;
    _uplo          = x._uplo;
    _rowMajorOrder = x._rowMajorOrder;
    assignment(x);
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix &Matrix::operator=(const Matrix &x)
{
  return operator=(static_cast<const const_MatrixSlice &>(x));
}

/***********************************************/
/***** Vector **********************************/
/***********************************************/

Vector::Vector(const const_MatrixSlice &x) : Matrix(x)
{
  try
  {
    if(x.columns() > 1)
      throw(Exception("Matrix assignment to Vector with more than 1 column"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector::Vector(std::initializer_list<Double> list) : Matrix(list.size(), 1)
{
  try
  {
    for(UInt i=0; i<list.size(); i++)
      (*this)(i) = *(list.begin()+i);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector::Vector(const std::vector<Double> &x) : Matrix(x.size(), 1)
{
  try
  {
    for(UInt i=0; i<x.size(); i++)
      (*this)(i) = x[i];
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector &Vector::operator=(const const_MatrixSlice &x)
{
  try
  {
    if(x.columns() > 1)
      throw(Exception("Matrix assignment to Vector with more than 1 column"));
    Matrix::operator=(x);
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector::operator std::vector<Double> () const
{
  try
  {
    std::vector<Double> x(rows());
    for(UInt i=0; i<x.size(); i++)
      x[i] = operator()(i);
    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** FUNCTIONS *******************************/
/***********************************************/

Matrix identityMatrix(UInt size, Matrix::Type type)
{
  Matrix I(size, type);
  for(UInt i=0; i<size; i++)
    I(i,i) = 1.0;
  return I;
}

/***********************************************/

void fillSymmetric(MatrixSliceRef A)
{
  try
  {
    if(A.getType()!=Matrix::SYMMETRIC)
      throw(Exception("Matrix must be SYMMETRIC"));
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));

    if(A.isUpper())
      for(UInt i=0; i<A.rows(); i++)
        copy(A.slice(i, i+1, 1, A.columns()-i-1).trans(), A.slice(i+1, i, A.rows()-i-1, 1));
    else
      for(UInt i=0; i<A.columns(); i++)
        copy(A.slice(i+1, i, A.rows()-i-1, 1).trans(), A.slice(i, i+1, 1, A.columns()-i-1));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void zeroUnusedTriangle(MatrixSliceRef A)
{
  try
  {
    if(A.getType()!=Matrix::TRIANGULAR)
      throw(Exception("Matrix must be TRIANGULAR"));
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));

    if(A.isUpper()) // set lower triangle to zero
      for(UInt k=0; k<A.columns()-1; k++)
        A.slice(k+1, k, A.rows()-k-1, 1).setNull();
    else  // set upper triangle to zero
      for(UInt k=1; k<A.columns(); k++)
        A.slice(k, k+1, 1, A.columns()-k-1).setNull();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/
//***********************************************/

void copy(const_MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if((A.size()==0) && (B.size()==0))
      return;
    if((A.rows()!=B.rows()) || (A.columns()!=B.columns()))
      throw(Exception("Dimension error: ("+A.rows()%"%i x "s+A.columns()%"%i) and ("s+B.rows()%"%i x "s+B.columns()%"%i)"s));

    const UInt rows    = (B.isRowMajorOrder()) ? B.columns() : B.rows();
    const UInt columns = (B.isRowMajorOrder()) ? B.rows() : B.columns();
    if(A.isRowMajorOrder() == B.isRowMajorOrder())
      lapack_dlacpy(rows, columns, A.field(),A.ld(), B.field(),B.ld());
    else
      for(UInt i=0; i<columns; i++)
        blas_dcopy(static_cast<F77Int>(rows), A.field()+i, static_cast<F77Int>(A.ld()), B.field()+i*B.ld(), 1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void swap(MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if((A.size()==0) && (B.size()==0))
      return;
    if((A.rows()!=B.rows()) || (A.columns()!=B.columns()))
      throw(Exception("Dimension error: ("+A.rows()%"%i x "s+A.columns()%"%i) and ("s+B.rows()%"%i x "s+B.columns()%"%i)"s));

    const UInt ldA  = (A.isRowMajorOrder()) ? 1 : A.ld();
    const UInt incA = (A.isRowMajorOrder()) ? A.ld() : 1;
    const UInt ldB  = (B.isRowMajorOrder()) ? 1 : B.ld();
    const UInt incB = (B.isRowMajorOrder()) ? B.ld() : 1;
    if(A.rows()>=A.columns())
      for(UInt i=0; i<A.columns(); i++)
        blas_dswap(static_cast<F77Int>(A.rows()), A.field()+i*ldA, static_cast<F77Int>(incA), B.field()+i*ldB, static_cast<F77Int>(incB));
    else
      for(UInt i=0; i<A.rows(); i++)
        blas_dswap(static_cast<F77Int>(A.columns()), A.field()+i*incA, static_cast<F77Int>(ldA), B.field()+i*incB, static_cast<F77Int>(ldB));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

//***********************************************/

void reshape(const_MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if(A.size() != B.size())
      throw(Exception("Dimension error: ("+A.rows()%"%i x "s+A.columns()%"%i) can't be reshaped to ("s+B.rows()%"%i x "s+B.columns()%"%i)"s));
    if(A.size()==0)
      return;

    if(A.rows() == B.rows())
      copy(A, B);
    else
      for(UInt i=0; i<A.size(); i++)
        B(i%B.rows(), i/B.rows()) = A(i%A.rows(), i/A.rows());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

//***********************************************/

Matrix reshape(const_MatrixSliceRef A, UInt rows, UInt columns)
{
  try
  {
    if(!rows && !columns)
      throw(Exception("Only one dimension can be automatically determined."));
    Matrix B(rows ? rows : A.size()/columns, columns ? columns : A.size()/rows);
    reshape(A, B);
    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

Matrix reorder(const_MatrixSliceRef A, const std::vector<UInt> &rowIndex)
{
  try
  {
    Matrix B(rowIndex.size(), A.columns());
    for(UInt i=0; i<B.rows(); i++)
      if(rowIndex[i] != NULLINDEX)
        copy(A.row(rowIndex.at(i)), B.row(i));
    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool isStrictlyZero(const_MatrixSliceRef A)
{
  try
  {
    for(UInt c = 0; c<A.columns(); c++)
      for(UInt r = 0; r<A.rows(); r++)
        if(A(r, c) != 0.0)
          return FALSE;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }

}

/***********************************************/

Double inner(const_MatrixSliceRef A, const_MatrixSliceRef B)
{
  try
  {
    if((A.size()==0) && (B.size()==0))
      return 0.;
    if((A.rows()!=B.rows()) || (A.columns()!=B.columns()))
      throw(Exception("Dimension error: ("+A.rows()%"%i x "s+A.columns()%"%i)  and ("s+B.rows()%"%i x "s+B.columns()%"%i)"s));

    const UInt rows    = (B.isRowMajorOrder()) ? B.columns() : B.rows();
    const UInt columns = (B.isRowMajorOrder()) ? B.rows() : B.columns();
    const Double *ptr1 = A.field();
    const Double *ptr2 = B.field();
    Double sum  = 0;

    if(A.isRowMajorOrder() == B.isRowMajorOrder())
      for(UInt i=0; i<columns; i++)
        sum += blas_ddot(static_cast<F77Int>(rows), ptr1+i*A.ld(), 1, ptr2+i*B.ld(), 1);
    else
      for(UInt i=0; i<columns; i++)
        sum += blas_ddot(static_cast<F77Int>(rows), ptr1+i, static_cast<F77Int>(A.ld()), ptr2+i*B.ld(), 1);
    return sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double sum(const_MatrixSliceRef A)
{
  try
  {
    Double sum = 0;
    for(UInt i=0; i<A.rows(); i++)
      for(UInt k=0; k<A.columns(); k++)
        sum += A(i,k);
    return sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double trace(const_MatrixSliceRef A)
{
  try
  {
    if(A.rows() != A.columns())
      throw(Exception("Dimension error: A("+A.rows()%"%i x "s+A.columns()%"%i) Matrix must be square!"s));

    Double tr = 0;
    for(UInt i=0; i<A.rows(); i++)
      tr += A(i,i);
    return tr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double determinant(const_MatrixSliceRef A)
{
  Double sign = 0.0;
  Double logdet = logdeterminant(A, sign);

  return sign * std::exp(logdet);
}

/***********************************************/

Double logdeterminant(const_MatrixSliceRef A, Double& sign)
{
  try
  {
    if(A.rows() != A.columns())
      throw(Exception("Determinant is only defined for square matrices: A("+A.rows()%"%i x "s+A.columns()%"%i)."s));

    sign = 0.0;
    Double logdet = -std::numeric_limits<Double>::infinity();

    if(A.getType() == Matrix::TRIANGULAR)
    {
      UInt countNegative = 0;
      logdet = 0.0;

      // eigenvalues of a triangular matrix are its diagonal elements
      for(UInt k = 0; k<A.rows(); k++)
      {
        logdet += std::log(std::abs(A(k, k)));
        if(A(k, k)<0)
          countNegative++;
      }
      sign = (countNegative % 2) == 0 ? 1.0 : -1.0;
    }
    else // GENERAL and SYMMETRIC matrices from QR decomposition: det(A) = det(Q)*det(R) = det(R)*(-1)^s
    {
      Matrix R = A; // copy for decomposition
      if(R.getType() == Matrix::SYMMETRIC)
      {
        fillSymmetric(R);
        R.setType(Matrix::GENERAL);
      }

      Vector tau = QR_decomposition(R); // A = QR
      Double signR = 0.0;
      logdet = logdeterminant(R, signR); // recursive call for det(R)
      UInt houseHolderCount = R.columns()-1; // det(Q) depends on number of reflectors: number of columns-1

      sign = (houseHolderCount % 2) == 0 ? signR : -signR;  // flip sign if necessary
    }

    return logdet;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double maxabs(const_MatrixSliceRef A)
{
  Double d = 0;
  for(UInt i=0; i<A.rows(); i++)
    for(UInt k=0; k<A.columns(); k++)
      d = std::max(d, std::fabs(A(i,k)));
  return A.size() ? d : NAN_EXPR;
}

/***********************************************/

Double min(const_MatrixSliceRef A)
{
  Double d = std::numeric_limits<double>::max();
  for(UInt i=0; i<A.rows(); i++)
    for(UInt k=0; k<A.columns(); k++)
      d = std::min(d, A(i,k));
  return A.size() ? d : NAN_EXPR;
}

/***********************************************/

Double max(const_MatrixSliceRef A)
{
  Double d = std::numeric_limits<double>::min();
  for(UInt i=0; i<A.rows(); i++)
    for(UInt k=0; k<A.columns(); k++)
      d = std::max(d, A(i,k));
  return A.size() ? d : NAN_EXPR;
}

/***********************************************/

Double median(const_MatrixSliceRef A)
{
  try
  {
    if(!A.size())
      return NAN_EXPR;
    std::vector<Double> data = flatten(A);
    std::partial_sort(data.begin(), data.begin()+data.size()/2+1, data.end());
    return (data.size()%2) ? data.at(data.size()/2) : (0.5*(data.at(data.size()/2-1)+data.at(data.size()/2)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }

}

/***********************************************/

void axpy(Double c, const_MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if((A.size()==0) && (B.size()==0))
      return;
    if((A.rows()!=B.rows()) || (A.columns()!=B.columns()))
      throw(Exception("Dimension error: ("+A.rows()%"%i x "s+A.columns()%"%i)  and ("s+B.rows()%"%i x "s+B.columns()%"%i)"s));

    const UInt rows    = (B.isRowMajorOrder()) ? B.columns() : B.rows();
    const UInt columns = (B.isRowMajorOrder()) ? B.rows() : B.columns();
    const Double *ptr1 = A.field();
          Double *ptr2 = B.field();

    if(A.isRowMajorOrder() == B.isRowMajorOrder())
      for(UInt i=0; i<columns; i++)
        blas_daxpy(static_cast<F77Int>(rows), c, ptr1+i*A.ld(), 1, ptr2+i*B.ld(), 1);
    else
      for(UInt i=0; i<columns; i++)
        blas_daxpy(static_cast<F77Int>(rows), c, ptr1+i, static_cast<F77Int>(A.ld()), ptr2+i*B.ld(), 1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Matrix-Matrix multiplication: C += c * A * B.
void matMult(Double c, const_MatrixSliceRef A, const_MatrixSliceRef B, MatrixSliceRef C)
{
  try
  {
    if((A.columns()!=B.rows())||(A.rows()!=C.rows())||(B.columns()!=C.columns()))
      throw(Exception("Dimension error"));
    if(C.size()==0)
      return;
    if(C.getType()!=Matrix::GENERAL)
      throw(Exception("Matrix C must be GENERAL"));

    if((A.getType()==Matrix::TRIANGULAR) && (B.getType()==Matrix::GENERAL))
    {
      Matrix B2 = B;
      triangularMult(c,A,B2);
      axpy(1., B2, C);
      return;
    }

    if((A.getType()==Matrix::GENERAL) && (B.getType()==Matrix::TRIANGULAR))
    {
      Matrix A2 = A;
      triangularMult(c, B.trans(), A2.trans());
      axpy(1., A2, C);
      return;
    }

    if((A.getType()==Matrix::SYMMETRIC) && (B.getType()==Matrix::GENERAL) && (C.isRowMajorOrder() == B.isRowMajorOrder()))
    {
      const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
      if(!C.isRowMajorOrder())
        blas_dsymm(TRUE, isUpper, static_cast<F77Int>(C.rows()), static_cast<F77Int>(C.columns()),
                   c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()), 1.0, C.field(), static_cast<F77Int>(C.ld()));
      else
        blas_dsymm(FALSE, isUpper, static_cast<F77Int>(C.columns()), static_cast<F77Int>(C.rows()),
                   c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()), 1.0, C.field(), static_cast<F77Int>(C.ld()));
      return;
    }

    if((A.getType()==Matrix::GENERAL) && (B.getType()==Matrix::SYMMETRIC) && (C.isRowMajorOrder() == A.isRowMajorOrder()))
    {
      const Bool isUpper = (B.isRowMajorOrder()) ? (!B.isUpper()) : B.isUpper();
      if(!C.isRowMajorOrder())
        blas_dsymm(FALSE, isUpper, static_cast<F77Int>(C.rows()), static_cast<F77Int>(C.columns()),
                   c, B.field(), static_cast<F77Int>(B.ld()), A.field(), static_cast<F77Int>(A.ld()), 1.0, C.field(), static_cast<F77Int>(C.ld()));
      else
        blas_dsymm(TRUE, isUpper, static_cast<F77Int>(C.columns()), static_cast<F77Int>(C.rows()),
                   c, B.field(), static_cast<F77Int>(B.ld()), A.field(), static_cast<F77Int>(A.ld()), 1.0, C.field(), static_cast<F77Int>(C.ld()));
      return;
    }

    if((A.getType()!=Matrix::GENERAL) || (B.getType()!=Matrix::GENERAL))
      throw(Exception("Combination of Matrix types not implemented"));

    if(!C.isRowMajorOrder())
      blas_dgemm(A.isRowMajorOrder(), B.isRowMajorOrder(), static_cast<F77Int>(C.rows()), static_cast<F77Int>(C.columns()), static_cast<F77Int>(A.columns()),
                 c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()), 1.0, C.field(), static_cast<F77Int>(C.ld()));
    else
      blas_dgemm(!B.isRowMajorOrder(), !A.isRowMajorOrder(), static_cast<F77Int>(C.columns()), static_cast<F77Int>(C.rows()), static_cast<F77Int>(A.columns()),
                 c, B.field(), static_cast<F77Int>(B.ld()), A.field(), static_cast<F77Int>(A.ld()), 1.0, C.field(), static_cast<F77Int>(C.ld()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("C("s+C.rows()%"%i x "s+C.columns()%"%i) = A("s+A.rows()%"%i x "s+A.columns()%"%i) * B("s+B.rows()%"%i x "s+B.columns()%"%i)"s, e)
  }
}

/***********************************************/

void rankKUpdate(Double c, const_MatrixSliceRef A, MatrixSliceRef N)
{
  try
  {
    if((N.rows()!=N.columns())||(A.columns()!=N.rows()))
      throw(Exception("Dimension error"));
    if((N.size()==0) || (A.size()==0))
      return;
    if((N.getType()!=Matrix::SYMMETRIC)||(A.getType()!=Matrix::GENERAL))
      throw(Exception("Matrix A must be GENERAL and Matrix N must be SYMMETRIC"));

    const Bool isUpper = (N.isRowMajorOrder()) ? (!N.isUpper()) : N.isUpper();
    blas_dsyrk(isUpper, !A.isRowMajorOrder(), static_cast<F77Int>(N.rows()), static_cast<F77Int>(A.rows()),
               c, A.field(), static_cast<F77Int>(A.ld()), 1.0, N.field(), static_cast<F77Int>(N.ld()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("N = A'A mit A = ("s+A.rows()%"%i x "s+A.columns()%"%i) and N = ("s+N.rows()%"%i x "s+N.columns()%"%i)"s, e)
  }
}

/***********************************************/

void rank2KUpdate(Double c, const_MatrixSliceRef A, const_MatrixSliceRef B, MatrixSliceRef N)
{
  try
  {
    if((N.rows()!=N.columns())||(A.columns()!=N.rows())||(B.columns()!=N.rows()))
      throw(Exception("Dimension error"));
    if(N.size()==0)
      return;
    if((N.getType()!=Matrix::SYMMETRIC)||(A.getType()!=Matrix::GENERAL)||(B.getType()!=Matrix::GENERAL))
      throw(Exception("Matrix A and B must be GENERAL and Matrix N must be SYMMETRIC"));
    if(A.isRowMajorOrder() != B.isRowMajorOrder())
      throw(Exception("not possible"));

    const Bool isUpper = (N.isRowMajorOrder()) ? (!N.isUpper()) : N.isUpper();
    blas_dsyr2k(isUpper, !A.isRowMajorOrder(), static_cast<F77Int>(N.rows()), static_cast<F77Int>(A.rows()),
                c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()), 1.0, N.field(), static_cast<F77Int>(N.ld()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("N = A'B+B'A mit A = ("s+A.rows()%"%i x "s+A.columns()%"%i), B = ("s+B.rows()%"%i x "s+B.columns()%"%i), N = ("s+N.rows()%"%i x "s+N.columns()%"%i)"s, e)
  }
}

/***********************************************/

void inverse(MatrixSliceRef A)
{
  try
  {
    if(A.size()==0)
      return;
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));

    Int info = 0;
    if(A.getType()==Matrix::TRIANGULAR)
    {
      const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
      info = lapack_dtrtri(isUpper, A.rows(),A.field(),A.ld());
    }
    else if(A.getType()==Matrix::SYMMETRIC)
    {
      cholesky(A);
      const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
      info = lapack_dpotri(isUpper, A.rows(),A.field(),A.ld());
      const_cast<MatrixSlice&>(A).setType(Matrix::SYMMETRIC);
    }
    else
    {
      Int *ipiv = new Int[A.rows()];
      info = lapack_dgetrf(A.rows(),A.columns(), A.field(),A.ld(),ipiv);
      if(info==0)
        info = lapack_dgetri(A.rows(),A.field(),A.ld(),ipiv);
      delete[] ipiv;
    }

    if(info!=0)
      throw(Exception("can not compute inverse, error = "+info%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Vector eigenValueDecomposition(MatrixSliceRef A, Bool computeEigenVectors)
{
  try
  {
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));
    if(A.getType()!=Matrix::SYMMETRIC)
      throw(Exception("Matrix must be SYMMETRIC"));

    Vector eigen(A.rows());
    const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    const Int info = lapack_dsyev(computeEigenVectors, isUpper, A.rows(),A.field(),A.ld(), eigen.field());
    const_cast<MatrixSlice&>(A).setType(Matrix::GENERAL);
    if(info!=0)
      throw(Exception("cannot compute eigen values, error = "+info%"%i"s));
    return eigen;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Matrix eigenValueDecomposition(MatrixSliceRef A, Matrix &VL, Matrix &VR, Bool computeEigenVectors)
{
  try
  {
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));
    if(A.getType()!=Matrix::GENERAL)
      throw(Exception("Matrix must be GENERAL"));

    Matrix eigenValues(A.rows(), 2);
    Int    info;
    if(computeEigenVectors)
    {
      VL = Matrix(A.rows(), A.rows());
      VR = Matrix(A.rows(), A.rows());
      info = lapack_dgeev(TRUE, TRUE, A.rows(), A.field(), A.ld(),   // Input matrix field
                          &eigenValues(0,0), &eigenValues(0,1),      // Real and complex eigenvalue vectors
                          VL.field(), VL.ld(), VR.field(), VR.ld()); // Left and right eigenvectors
    }
    else
    {
      info = lapack_dgeev(FALSE, FALSE, A.rows(), A.field(), A.ld(), // Input matrix field
                          &eigenValues(0,0), &eigenValues(0,1),      // Real and complex eigenvalue vectors
                          nullptr, 1, nullptr, 1);
    }

    if(info!=0)
      throw(Exception("cannot compute eigen values, error = "+info%"%i"s));
    return eigenValues;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Vector singularValueDecomposition(MatrixSliceRef A, Matrix &U, Matrix &Vt, Bool computeSingularVectors)
{
  try
  {
    if(A.isRowMajorOrder())
      throw(Exception("Matrix must not be transposed"));

    const UInt m = A.rows();
    const UInt n = A.columns();
    const UInt s = std::min(m,n);
    Vector singular(s);
    Int    info;
    if(computeSingularVectors)
    {
      U  = Matrix(m,s);
      Vt = Matrix(s,n);
      info = lapack_dgesdd(TRUE, m, n, A.field(), A.ld(), singular.field(), U.field(), U.ld(), Vt.field(), Vt.ld());
    }
    else
      info = lapack_dgesdd(FALSE, m, n, A.field(), A.ld(), singular.field(), nullptr, 1, nullptr, 1);

    if(info!=0)
      throw(Exception("cannot compute singular values, error = "+info%"%i"s));
    return singular;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Matrix pseudoInverse(const_MatrixSliceRef A, Double rcond)
{
  try
  {
    // symmetric positive semidefinite
    if(A.getType() == Matrix::SYMMETRIC)
    {
      Matrix U = A; // copy contents of A
      Vector s = eigenValueDecomposition(U);

      const Double maxS = maxabs(s);
      for(UInt k=0; k<s.rows(); k++)
        U.column(k) *= (s(k) >= (rcond*maxS)) ? 1./sqrt(s(k)) : 0.;

      Matrix Inv(A.rows(), Matrix::SYMMETRIC);
      rankKUpdate(1., U.trans(), Inv);
      return Inv;
    }

    // general inverse
    Matrix U, Vt;
    Vector s = singularValueDecomposition(Matrix(A), U, Vt);

    const Double maxS = maxabs(s);
    for(UInt k=0; k<s.rows(); k++)
      U.column(k) *= (std::fabs(s(k)) >= (rcond*maxS)) ? 1./s(k) : 0.;

    return Vt.trans() * U.trans();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Matrix matrixSquareRoot(const_MatrixSliceRef A, Double rcond)
{
  try
  {
    if(A.getType() != Matrix::SYMMETRIC)
      throw(Exception("Input matrix must be SYMMETRIC."));

    Matrix Q = A; // copy contents of A
    Matrix U, Vt;
    Vector e = singularValueDecomposition(Q, U, Vt, TRUE);

    const Double maxE = maxabs(e);
    for(UInt k=0; k<e.rows(); k++)
    {
      if(e(k) < -(rcond*maxE))
        throw(Exception("Input matrix must be positive semidefinite"));
      U.column(k) *= (e(k) >= (rcond*maxE)) ? std::sqrt(e(k)) : 0.;
    }

    return U*Vt; // Q*sqrt(E)*Q.trans()
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Matrix matrixSquareRootInverse(const_MatrixSliceRef A, Double rcond)
{
  try
  {
    if(A.getType() != Matrix::SYMMETRIC)
      throw(Exception("Input matrix must be SYMMETRIC."));

    Matrix Q = A; // copy contents of A
    Matrix U, Vt;
    Vector e = singularValueDecomposition(Q, U, Vt, TRUE);

    const Double maxE = maxabs(e);
    for(UInt k=0; k<e.rows(); k++)
    {
      if(e(k) < -(rcond*maxE))
        throw(Exception("Input matrix must be positive semidefinite"));
      U.column(k) *= (e(k) >= (rcond*maxE)) ? 1.0/std::sqrt(e(k)) : 0.;
    }

    return U*Vt; // Q*E^(-0.5)*Q.trans()
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

Matrix kron(const_MatrixSliceRef A, const_MatrixSlice B)
{
  try
  {
    Matrix C(A.rows()*B.rows(), A.columns()*B.columns());
    for(UInt i=0; i<A.rows(); i++)
      for(UInt k=0; k<A.columns(); k++)
        axpy(A(i, k), B, C.slice(i*B.rows(), k*B.columns(), B.rows(), B.columns()));
    return C;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix toeplitz(const Vector &a)
{
  try
  {
    Matrix A(a.rows(), Matrix::SYMMETRIC);
    for(UInt i=0; i<A.columns(); i++)
      copy(a.row(0, a.rows()-i), A.slice(i, i, A.rows()-i, 1));
    fillSymmetric(A);
    return A;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void solveInPlace(MatrixSliceRef N, MatrixSliceRef B)
{
  try
  {
    if(N.rows()!=N.columns()||(N.columns()!=B.rows()))
      throw(Exception("Dimension error"));
    if(B.size()==0)
      return;
    if(B.getType()!=Matrix::GENERAL)
      throw(Exception("Matrix B must be GENERAL"));

    if(N.getType()==Matrix::TRIANGULAR)
    {
      triangularSolve(1.0, N, B);
    }
    else if(N.getType()==Matrix::SYMMETRIC)
    {
      cholesky(N);
      triangularSolve(1.0, N.trans(), B);
      triangularSolve(1.0, N, B);
    }
    else if(N.getType()==Matrix::GENERAL)
    {
      if(N.isRowMajorOrder() || B.isRowMajorOrder())
        throw(Exception("Matricies must not be transposed, if N == GENERAL"));
      Int *ipiv = new Int[N.rows()];
      const Int info = lapack_dgesv(B.rows(),B.columns(), N.field(),N.ld(), ipiv, B.field(),B.ld());
      delete[] ipiv;
      if(info!=0)
        throw(Exception("cannot solve system, error = "+info%"%i"s));
    }
    else
      throw(Exception("Matrix Type not implemented"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("B := N("s+N.rows()%"%i x "s+N.columns()%"%i)^-1 * B = ("s+B.rows()%"%i x "s+B.columns()%"%i)"s, e)
  }
}

/***********************************************/

Matrix leastSquares(MatrixSliceRef A, MatrixSliceRef l)
{
  try
  {
    if((A.getType()!=Matrix::GENERAL)||(l.getType()!=Matrix::GENERAL))
      throw(Exception("Matrix A and l must be GENERAL"));
    if((A.rows()<A.columns())||(A.rows()!=l.rows()))
      throw(Exception("Dimension error"));
    if(l.isRowMajorOrder())
      throw(Exception("l must not be transposed"));

    Vector tau = QR_decomposition(A);
    QTransMult(A, tau, l);
    Matrix x = l.row(0, tau.rows());
    triangularSolve(1., A.row(0, tau.rows()), x);
    l.row(0, tau.rows()).setNull();
    QMult(A, tau, l);

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("min[ A("s+A.rows()%"%i x "s+A.columns()%"%i)*x - l("s+l.rows()%"%i x "s+l.columns()%"%i)"s, e)
  }
}

/***********************************************/

void reduceLeastSquaresFit(const_MatrixSliceRef A, MatrixSliceRef l)
{
  try
  {
    if((A.getType()!=Matrix::GENERAL)||(l.getType()!=Matrix::GENERAL))
      throw(Exception("Matrix A and l must be GENERAL"));
    if((A.rows()<A.columns())||(A.rows()!=l.rows()))
      throw(Exception("Dimension error"));

    Matrix B = A;
    Vector tau = QR_decomposition(B);
    QTransMult(B, tau, l);
    l.row(0,tau.rows()).setNull();
    QMult(B, tau, l);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("min[ A("s+A.rows()%"%i x "s+A.columns()%"%i)*x - l("s+l.rows()%"%i x "s+l.columns()%"%i)"s, e)
  }
}

/***********************************************/

void eliminationParameter(MatrixSliceRef B, const std::vector<std::reference_wrapper<Matrix>> &listA)
{
  try
  {
    const Vector tau = QR_decomposition(B);
    for(Matrix &A : listA)
    {
      if(A.rows() != B.rows())
        throw(Exception("Dimension error: B("s+B.rows()%"%i x "s+B.columns()%"%i), A("s+A.rows()%"%i x "s+A.columns()%"%i)"s));
      QTransMult(B, tau, A);
      A = A.row(B.columns(), A.rows()-B.columns());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

// Triangular matrix
// -----------------
void cholesky(MatrixSliceRef A)
{
  try
  {
    if(A.getType()!=Matrix::SYMMETRIC)
      throw(Exception("Matrix must be SYMMETRIC"));
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));

    const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    const Int info = lapack_dpotrf(isUpper, A.rows(), A.field(),A.ld());
    const_cast<MatrixSlice&>(A).setType(Matrix::TRIANGULAR, Matrix::UPPER);
    if(isUpper)
      const_cast<MatrixSlice&>(A)._rowMajorOrder = FALSE;
    else
      const_cast<MatrixSlice&>(A)._rowMajorOrder = TRUE;
    if(info!=0)
      throw(Exception("cannot compute decomposition, error = "+info%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("W'W = A("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

std::vector<UInt> choleskyPivoting(MatrixSliceRef A, UInt &rank, Double tolerance)
{
  try
  {
    if(A.getType()!=Matrix::SYMMETRIC)
      throw(Exception("Matrix must be SYMMETRIC"));
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));

    F77Int *ipiv = new F77Int[A.rows()];
    F77Int  irank;
    const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    const Int info = lapack_dpstrf(isUpper, A.rows(), A.field(), A.ld(), ipiv, irank, tolerance);
    const_cast<MatrixSlice&>(A).setType(Matrix::TRIANGULAR, Matrix::UPPER);
    if(isUpper)
      const_cast<MatrixSlice&>(A)._rowMajorOrder = FALSE;
    else
      const_cast<MatrixSlice&>(A)._rowMajorOrder = TRUE;

    if(info<0)
      throw(Exception("cannot compute decomposition, error = "+info%"%i"s));

    rank = irank;
    std::vector<UInt> piv(A.rows());
    for(UInt i=0; i<piv.size(); i++)
      piv[i] = static_cast<UInt>(ipiv[i]-1);
    delete[] ipiv;
    return piv;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("W'W = A("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void cholesky2Inverse(MatrixSliceRef A)
{
  try
  {
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));
    if(A.size()==0)
      return;
    if(A.getType()!=Matrix::TRIANGULAR)
      throw(Exception("Matrix must be TRIANGULAR"));

    const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    const Int info = lapack_dpotri(isUpper, A.rows(),A.field(),A.ld());
    const_cast<MatrixSlice&>(A).setType(Matrix::SYMMETRIC);
    if(info!=0)
      throw(Exception("can not compute inverse, error = "+info%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void choleskyProduct(MatrixSliceRef A)
{
  try
  {
    if(A.rows()!=A.columns())
      throw(Exception("Dimension error"));
    if(A.getType()!=Matrix::TRIANGULAR)
      throw(Exception("Matrix must be TRIANGULAR"));
    if(A.size()==0)
      return;

    const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    const Int info = lapack_dlauum(isUpper, A.rows(), A.field(),A.ld());
    const_cast<MatrixSlice&>(A).setType(Matrix::SYMMETRIC);
    if(info!=0)
      throw(Exception("cannot compute W'W, error = "+info%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void triangularMult(Double c, const_MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if((A.getType()!=Matrix::TRIANGULAR)||(B.getType()!=Matrix::GENERAL))
      throw(Exception("Matrix A must be TRIANGULAR and Matrix B must be GENERAL"));
    if((A.rows()!=A.columns())||(A.columns()!=B.rows()))
      throw(Exception("Dimension error"));

    Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    if(!B.isRowMajorOrder())
      blas_dtrmm(TRUE,  isUpper, A.isRowMajorOrder(), FALSE, static_cast<F77Int>(B.rows()), static_cast<F77Int>(B.columns()),
                 c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()));
    else
      blas_dtrmm(FALSE, isUpper, !A.isRowMajorOrder(), FALSE, static_cast<F77Int>(B.columns()), static_cast<F77Int>(B.rows()),
                 c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("B = A("s+A.rows()%"%i x "s+A.columns()%"%i) * B("s+B.rows()%"%i x "s+B.columns()%"%i)"s, e)
  }
}

/***********************************************/

void triangularSolve(Double c, const_MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if((A.getType()!=Matrix::TRIANGULAR)||(B.getType()!=Matrix::GENERAL))
      throw(Exception("Matrix A must be TRIANGULAR and Matrix B must be GENERAL"));
    if((A.rows()!=A.columns())||(A.columns()!=B.rows()))
      throw(Exception("Dimension error"));

    const Bool isUpper = (A.isRowMajorOrder()) ? (!A.isUpper()) : A.isUpper();
    if(!B.isRowMajorOrder())
      blas_dtrsm(TRUE,  isUpper,  A.isRowMajorOrder(), FALSE, static_cast<F77Int>(B.rows()), static_cast<F77Int>(B.columns()),
                c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()));
    else
      blas_dtrsm(FALSE, isUpper, !A.isRowMajorOrder(), FALSE, static_cast<F77Int>(B.columns()), static_cast<F77Int>(B.rows()),
                 c, A.field(), static_cast<F77Int>(A.ld()), B.field(), static_cast<F77Int>(B.ld()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("B = A("s+A.rows()%"%i x "s+A.columns()%"%i)^-1 * B("s+B.rows()%"%i x "s+B.columns()%"%i)"s, e)
  }
}

/***********************************************/
/***********************************************/

Vector QR_decomposition(MatrixSliceRef A)
{
  try
  {
    if(A.isRowMajorOrder())
      throw(Exception("A.trans not implemented"));

    Vector tau (A.columns());
    const Int info = lapack_dgeqrf(A.rows(), A.columns(), A.field(), A.ld(), tau.field());
    if(info!=0)
      throw(Exception("can not compute QR decomposition, error = "+info%"%i"s));
    const_cast<MatrixSlice&>(A).setType(Matrix::TRIANGULAR, Matrix::UPPER);
    return tau;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void QMult(const_MatrixSliceRef B, const Vector &tau, MatrixSliceRef A)
{
  try
  {
    if((B.rows()!=A.rows())||(tau.rows()!=B.columns()))
      throw(Exception("Dimension error"));
    if(A.getType()!=Matrix::GENERAL)
      throw(Exception("Matrix A must be GENERAL"));

    if(!A.isRowMajorOrder())
      lapack_dormqr(TRUE, FALSE, A.rows(),A.columns(),B.columns(), B.field(),B.ld(),tau.field(), A.field(),A.ld());
    else
      lapack_dormqr(FALSE, TRUE, A.columns(),A.rows(),B.columns(), B.field(),B.ld(),tau.field(), A.field(),A.ld());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Q("s+B.rows()%"%i x "s+B.columns()%"%i) * A("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void QTransMult(const_MatrixSliceRef B, const Vector &tau, MatrixSliceRef A)
{
  try
  {
    if((B.rows()!=A.rows())||(tau.rows()!=B.columns()))
      throw(Exception("Dimension error"));
    if(A.getType()!=Matrix::GENERAL)
      throw(Exception("Matrix C must be GENERAL"));

    if(!A.isRowMajorOrder())
      lapack_dormqr(TRUE, TRUE, A.rows(),A.columns(),B.columns(), B.field(),B.ld(),tau.field(), A.field(),A.ld());
    else
      lapack_dormqr(FALSE, FALSE, A.columns(),A.rows(),B.columns(), B.field(),B.ld(),tau.field(), A.field(),A.ld());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Q("s+B.rows()%"%i x "s+B.columns()%"%i) * A("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}

/***********************************************/

void generateQ(MatrixSliceRef A, const Vector &tau)
{
  try
  {
    const Int info = lapack_dorgqr(A.rows(), A.columns(), tau.rows(), A.field(), A.ld(), tau.field());
    if(info!=0)
      throw(Exception("can not compute Q matrix, error = "+info%"%i"s));
    const_cast<MatrixSlice&>(A).setType(Matrix::GENERAL);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("A = ("s+A.rows()%"%i x "s+A.columns()%"%i)"s, e)
  }
}
