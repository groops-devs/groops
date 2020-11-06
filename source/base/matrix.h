/***********************************************/
/**
* @file matrix.h
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

#ifndef __GROOPS_MATRIX__
#define __GROOPS_MATRIX__

#include "base/importStd.h"
#include "base/constants.h"

/** @defgroup matrixGroup Matrix and Vector
* @brief User friendly Classes for BLAS (Basic Linear Algebra Subroutines) and LAPACK (Linear Algebra Package).
* @ingroup base */
/// @{

/***** DEFINES *********************************/

class MatrixBase;
class const_MatrixSlice;
class MatrixSlice;
class Matrix;
class Vector;

typedef const const_MatrixSlice &const_MatrixSliceRef; //!< reason is given in @ref const_MatrixSlice
typedef const MatrixSlice       &MatrixSliceRef;       //!< reason is given in @ref MatrixSlice

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/** @brief Readonly view of a slice of a Matrix.
* This is a reference to a memory area of a Matrix. It is intentend to be used as parameter of a function,
* which operates on parts or all elements of the matrix without changing the size of the matrix.
* Examples are copy, matMult or inverse. To store a matrix in a variable use the class Matrix instead.
*
* The class MatrixSlice itself should be a const reference in a function parameter definition.
* Thus the function accepts temporally created const_MatrixSlice or Matrix as parameter.
* As the following definition is a little bit confusing:
* @code void f(const const_MatrixSlice &x); @endcode
* a shortcut is defined with the same meaning:
* @code void f(const_MatrixSliceRef x); @endcode */
class const_MatrixSlice
{
public:
  /// Matrix types
  enum Type {GENERAL, SYMMETRIC, TRIANGULAR};

  /** @brief upper or lower triangle.
  * For SYMMETRIC or TRIANGULAR matrices.
  * The other triangular is not acessed. */
  enum Uplo {UPPER, LOWER};

  const_MatrixSlice()                          = delete;  //!< Disallow default constructor.
  const_MatrixSlice(const const_MatrixSlice &) = default; //!< Copy constructor.
  const_MatrixSlice &operator=(const const_MatrixSlice &) = delete; //!< Disallow copying.

  UInt rows()    const {return _rows;}          //!< row count.
  UInt columns() const {return _columns;}       //!< column count.
  UInt size()    const {return _rows*_columns;} //!< size = rows*columns.

  inline const Double &operator () (UInt row, UInt column) const; /// matrix element.

  /** @brief Transposed matrix.
  * The transpose points to the same memory. */
  inline const_MatrixSlice trans() const;

  /** @brief Readonly reference to a submatrix.
  * The same memory is used. No memory is copied. */
  const_MatrixSlice slice(UInt startRow, UInt startColumn, UInt height, UInt width) const;

  /** @brief Readonly reference to one or more columns.
  * The same memory is used. No memory is copied.
  * @relates slice */
  inline const_MatrixSlice column(UInt column, UInt len=1) const {return slice(0,column,rows(),len);}

  /** @brief Readonly reference to one or more rows.
  * The same memory is used. No memory is copied.
  * @relates slice */
  inline const_MatrixSlice row(UInt row, UInt len=1) const {return slice(row,0,len,columns());}

  /** @brief Typ der Matrix.
  * GENERAL, SYMMETRIC, TRIANGLUAR. */
  Type getType() const {return _type;}

  /** @brief Which triangle is used?
  * For @a getType()==GENERAL the result is meaningless.
  * If TRUE the upper triangle is used, the lower otherwise. */
  Bool isUpper() const {return _uplo==UPPER;}

  /** @brief Pointer to the memory.
  * Acess to the element `A(row,column)` with
  * `ptr = (A.isRowMajorOrder()) ? (A.field() + (column + row*A.ld())) : (A.field() + (row + column*A.ld()));` */
  inline const Double *field() const;

  /** @brief Leading Dimension.
  * Number of elements needed to reach the next column or row,
  * depending wether `!isRowMajorOrder()` or `isRowMajorOrder()`.
  * @see field() */
  UInt ld() const {return _ld;}

  /** @brief Transposed view of the matrix?
  * If `isRowMajorOrder()==TRUE` the matrix is stored rowwise, columnwise otherwise.
  * @see field() */
  Bool isRowMajorOrder() const {return _rowMajorOrder;}

  inline Matrix operator+(Double c)                   const;
  inline Matrix operator-(Double c)                   const;
  inline Matrix operator*(Double c)                   const;
  inline Matrix operator/(Double c)                   const;
  inline Matrix operator-()                           const;
  inline Matrix operator+(const const_MatrixSlice &x) const;
  inline Matrix operator-(const const_MatrixSlice &x) const;

protected:
  UInt _rows,  _columns;
  UInt _start, _ld;
  Type _type;
  Uplo _uplo;
  Bool _rowMajorOrder;
  std::shared_ptr<MatrixBase> base;

  const_MatrixSlice(UInt rows, UInt columns, const_MatrixSlice::Type type, const_MatrixSlice::Uplo uplo, Double fill=0);

  friend class MatrixSlice;
  friend class Matrix;
  friend class Vector;
  friend void  cholesky(MatrixSliceRef A);
  friend std::vector<UInt> choleskyPivoting(MatrixSliceRef A, UInt& rank, Double tolerance);
}; // class MatrixSlice

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/** @brief Writable view of a slice of a Matrix.
* This is a reference to a memory area of a Matrix. It is intentend to be used as parameter of a function,
* which operates on parts or all elements of the matrix without changing the size of the matrix.
* Examples are copy, matMult or inverse. To store a matrix in a variable use the class Matrix instead.
*
* The class MatrixSlice itself should be a const reference in a function parameter definition,
* although the memory is still changeable (the elements). Thus the function accepts temporally created
* MatrixSlice or Matrix as parameter. As the following definition is confusing:
* @code void f(const MatrixSlice &x); // the elements of x are still writable @endcode
* a shortcut is defined with the same meaning:
* @code void f(MatrixSliceRef x); @endcode */
class MatrixSlice : public const_MatrixSlice
{
public:
  MatrixSlice() = delete;                                           //!< Disallow default constructor.
  MatrixSlice(const MatrixSlice &) = default;                       //!< Copy constructor.
  MatrixSlice(const const_MatrixSlice &x) : const_MatrixSlice(x) {} //!< Copy constructor.
  MatrixSlice &operator=(MatrixSlice &) = delete;                   //!< Disallow copying

  /// writable matrix element.
  inline Double &operator () (UInt row, UInt column) const;

  /** @brief Transposed matrix.
  * The transpose points to the same memory. */
  inline MatrixSlice trans() const;

  /** @brief Writable reference to a submatrix.
  * The same memory is used. No memory is copied. */
  MatrixSlice slice(UInt startRow, UInt startColumn, UInt height, UInt width) const;

  /** @brief Writable reference to one or more columns.
  * The same memory is used. No memory is copied.
  * @relates slice */
  MatrixSlice column(UInt column, UInt len=1) const {return slice(0,column,rows(),len);}

  /** @brief Writable reference to one or more rows.
  * The same memory is used. No memory is copied.
  * @relates slice */
  MatrixSlice row(UInt row, UInt len=1) const {return slice(row,0,len,columns());}

  /** @brief Set type of the matrix.
  * @param type GENERAL, SYMMETRIC, TRIANGLUAR. */
  void setType(Type type) {_type=type;}

  /** @brief Set type of the matrix.
  * @param type GENERAL, SYMMETRIC, TRIANGULAR.
  * @param uplo if type is not GENERAL, only the UPPER or LOWER part of the matrix is used. */
  void setType(Type type, Uplo uplo) {_type=type; _uplo=uplo;}

  /** @brief Pointer to the memory.
  * Acess to the element `A(row,column)` with
  * `ptr = (A.isRowMajorOrder()) ? (A.field() + (column + row*A.ld())) : (A.field() + (row + column*A.ld()));` */
  inline Double *field() const;

  /// fill all matrix elements.
  const MatrixSlice &fill(Double f) const;

  /// reset all matrix elements.
  const MatrixSlice &setNull() const {return fill(0);}

  const MatrixSlice &operator+=(const const_MatrixSlice &x) const;
  const MatrixSlice &operator-=(const const_MatrixSlice &x) const;
  const MatrixSlice &operator+=(Double c) const;
  const MatrixSlice &operator-=(Double c) const;
  const MatrixSlice &operator*=(Double c) const;
  const MatrixSlice &operator/=(Double c) const;

protected:
  MatrixSlice(UInt rows, UInt columns, MatrixSlice::Type type, MatrixSlice::Uplo uplo, Double fill=0) : const_MatrixSlice(rows, columns, type, uplo, fill) {}
}; // class MatrixSlice

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/** @brief Matrix representation.
* With smart memmory management. The memory is only copied if needed (copy on write, late copy).
* For parameter references use MatrixSlice or const_MatrixSlice. */
class Matrix : public MatrixSlice
{
public:
  /// Default Constructor.
  explicit Matrix(UInt rows=0, UInt columns=0, Double fill=0.) : MatrixSlice(rows, columns, GENERAL, UPPER, fill) {}

  /** @brief Constructor for quadratic matrices.
  * @param rows Number of rows and columns.
  * @param type GENERAL, SYMMETRIC, TRIANGULAR
  * @param uplo Only the UPPER or LOWER part of the matrix is used */
  Matrix(UInt rows, Type type, Uplo uplo=Matrix::UPPER) : MatrixSlice(rows, rows, type, uplo) {}

  Matrix(const Matrix &x);                                           //!< Copy Constructor.
  Matrix(const const_MatrixSlice &x);                                //!< Copy Constructor.
  Matrix(std::initializer_list<std::initializer_list<Double>> list); //!< List Constructor.
  Matrix &operator=(const Matrix &x);                                //!< Assignment.
  Matrix &operator=(const const_MatrixSlice &x);                     //!< Assignment.

  /// matrix element.
  Double operator()(UInt row, UInt column) const {return const_MatrixSlice::operator()(row,column);}

  /// writable matrix element.
  Double &operator()(UInt row, UInt column) {return MatrixSlice::operator()(row,column);}

  /** @brief Readonly transposed matrix.
  * The transpose points to the same memory. */
  const const_MatrixSlice trans() const {return const_MatrixSlice::trans();}

  /** @brief Writable transposed matrix.
  * The transpose points to the same memory. */
  MatrixSlice trans() {return MatrixSlice::trans();}

  /** @brief Readonly reference to a submatrix.
  * The same memory is used. No memory is copied. */
  const_MatrixSlice slice(UInt startRow, UInt startColumn, UInt height, UInt width) const {return const_MatrixSlice::slice(startRow, startColumn, height, width);}

  /** @brief Writable reference to a submatrix.
  * The same memory is used. No memory is copied. */
  MatrixSlice slice(UInt startRow, UInt startColumn, UInt height, UInt width) {return MatrixSlice::slice(startRow, startColumn, height, width);}

  /** @brief Readonly reference to one or more columns.
  * The same memory is used. No memory is copied.
  * @relates slice */
  const_MatrixSlice column(UInt column, UInt len=1) const {return slice(0,column,rows(),len);}

  /** @brief Writable reference to one or more columns.
  * The same memory is used. No memory is copied.
  * @relates slice */
  MatrixSlice column(UInt column, UInt len=1) {return slice(0,column,rows(),len);}

  /** @brief Readonly reference to one or more rows.
  * The same memory is used. No memory is copied.
  * @relates slice */
  const_MatrixSlice row(UInt row, UInt len=1) const {return slice(row,0,len,columns());}

  /** @brief Writable reference to one or more rows.
  * The same memory is used. No memory is copied.
  * @relates slice */
  MatrixSlice row(UInt row, UInt len=1) {return slice(row,0,len,columns());}

private:
  void assignment(const const_MatrixSlice &x);
};

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/** @brief Vector representation.
* A Vector is a Matrix with only one column. */
class Vector : public Matrix
{
public:
  /// Default Constructor (vector with zeros).
  explicit Vector(UInt rows=0, Double fill=0.) : Matrix(rows, 1, fill) {}

  Vector(const Vector &) = default;              //!< Copy Constructor.
  Vector(const const_MatrixSlice &x);            //!< Copy Constructor.
  Vector(std::initializer_list<Double> list);    //!< List Constructor.
  Vector(const std::vector<Double> &x);          //!< Cast operator.
//   explicit Vector(const std::vector<Time> &x);   //!< Cast operator.
  Vector &operator=(const Vector &) = default;   //!< Assignment.
  Vector &operator=(const const_MatrixSlice &x); //!< Assignment.

  inline Double operator()(UInt row) const {return const_MatrixSlice::operator()(row,0);} //!< vector element.
  inline Double operator[](UInt row) const {return const_MatrixSlice::operator()(row,0);} //!< vector element.
  inline Double at(UInt row)         const {return const_MatrixSlice::operator()(row,0);} //!< vector element.

  inline Double &operator()(UInt row) {return MatrixSlice::operator()(row,0);} //!< writable vector element.
  inline Double &operator[](UInt row) {return MatrixSlice::operator()(row,0);} //!< writable vector element.
  inline Double &at(UInt row)         {return MatrixSlice::operator()(row,0);} //!< writable vector element.

  /** @brief Writable reference to a submatrix.
  * The same memory is used. No memory is copied. */
  inline const_MatrixSlice slice(UInt startRow, UInt len) const {return const_MatrixSlice::slice(startRow,0,len,columns());}

  /** @brief Writable reference to a submatrix.
  * The same memory is used. No memory is copied. */
  inline MatrixSlice slice(UInt startRow, UInt len) {return MatrixSlice::slice(startRow,0,len,columns());}

  /// Cast operator to std::vector.
  operator std::vector<Double> () const;
}; // class Vector

/***********************************************/
/***** FUNCTIONS *******************************/
/***********************************************/

/** @brief Element-wise scalar matrix addition.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Matrix operator+(Double c, const_MatrixSliceRef A) {return A+c;}

/** @brief Element-wise scalar matrix subtraction.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Matrix operator-(Double c, const_MatrixSliceRef A);

/** @brief Element-wise scalar matrix multiplication.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Matrix operator*(Double c, const_MatrixSliceRef A) {return A*c;}

/** @brief Element-wise scalar matrix division.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Matrix operator/(Double c, const_MatrixSliceRef A);

/** @brief Matrix matrix multiplication.
* The @a type of matrices are considered. */
inline Matrix operator*(const_MatrixSliceRef A, const_MatrixSliceRef B);

/** @brief identity matrix. */
Matrix identityMatrix(UInt size, Matrix::Type=Matrix::GENERAL);

/** @brief Fill the not used triangle, the result is symmetric.
* A must be quadratic. */
void fillSymmetric(MatrixSliceRef A);

/** @brief Fill unused triangle with zeros.
* A must be quadratic and triangular, otherwise an exception is thrown. */
void zeroUnusedTriangle(MatrixSliceRef A);

/** @brief Copy the content of Matrix @a A in Matrix @a B.
* @warning Undefined for overlapping memory.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
void copy(const_MatrixSliceRef A, MatrixSliceRef B);

/** @brief Swap the content of Matrix @a A and Matrix @a B.
* @warning Undefined for overlapping memory.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
void swap(MatrixSliceRef A, MatrixSliceRef B);

/** @brief Reshape the content of Matrix @a A into Matrix @a B.
* Reshaping is done by column-wise order.
* @warning Undefined for overlapping memory.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
void reshape(const_MatrixSliceRef A, MatrixSliceRef B);

/** @brief Return a reshaped copy of Matrix @a A with new dimensions (@a rows x @a columns).
* Reshaping is done by column-wise order. Either @a rows or @a columns can be set to zero
* to auto-determine the dimension based on the provided one.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Matrix reshape(const_MatrixSliceRef A, UInt rows, UInt columns);

/** @brief Return a reordered copy of Matrix @a A.
* The parameters @a rowIndex and @a columnIndex contain the indices of the input matrix elements in the reordered matrix.
* Indices with NULLINDEX are inserted as zero elements into the new matrix. If either @a rowIndex or @a columnIndex are empty (size == 0),
* the ordering of the input matrix along the dimension is preserved.
* The @a type of the return matrix is the same as of the input matrix.
* (operation is always applied to both triangles). */
Matrix reorder(const_MatrixSliceRef A, const std::vector<UInt> &rowIndex, const std::vector<UInt> &columnIndex={});

/** @brief Return the columns of the matrix as column vector.
* The columns of the matrix are copied into the resulting vector.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Vector flatten(const_MatrixSliceRef A) {return reshape(A, A.size(), 1);}

/** @brief Return whether the matrix is contains only strict zeros.
* (operation is always applied to both triangles). */
Bool isStrictlyZero(const_MatrixSliceRef A);

/** @brief Scalar product.
* Sum of elementwise products.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Double inner(const_MatrixSliceRef A, const_MatrixSliceRef B);

/** @brief Sum of matrix elements. Returns 0 if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Double sum(const_MatrixSliceRef A);

/** @brief Mean of matrix elements. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Double mean(const_MatrixSliceRef A) {return A.size() ? sum(A)/A.size() : NAN_EXPR;}

/** @brief Square sum of matrix elements. Returns 0 if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Double quadsum(const_MatrixSliceRef A) {return inner(A,A);}

/** @brief Square root of squared sum of matrix elements. Returns 0 if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Double norm(const_MatrixSliceRef A) {return sqrt(quadsum(A));}

/** @brief Root mean square (RMS) of matrix elements. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Double rootMeanSquare(const_MatrixSliceRef A) {return A.size() ? std::sqrt(quadsum(A)/A.size()) : NAN_EXPR;}

/** @brief Standard deviation of matrix elements. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
inline Double standardDeviation(const_MatrixSliceRef A) {return A.size()>1 ? std::sqrt((quadsum(A)-std::pow(sum(A),2)/A.size())/(A.size()-1)) : NAN_EXPR;}

/** @brief Trace of matrix. Returns 0 if @p A is empty. */
Double trace(const_MatrixSliceRef A);

/** @brief Matrix determinant. For non-triangular matrices, the determinant is computed from the QR factorization of a copy of @p A. */
Double determinant(const_MatrixSliceRef A);

/** @brief Logarithm of matrix determinant and its sign. For non-triangular matrices, the determinant is computed from the QR factorization of a copy of @p A. */
Double logdeterminant(const_MatrixSliceRef A, Double &sign);

/** @brief Maximum absolute element of matrix. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Double maxabs(const_MatrixSliceRef A);

/** @brief Minimum element of matrix. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Double min(const_MatrixSliceRef A);

/** @brief Maximum element of matrix. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Double max(const_MatrixSliceRef A);

/** @brief Median element of matrix. Returns NAN if @p A is empty.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
Double median(const_MatrixSliceRef A);

/** @brief axpy: B += c * A.
* The @a type of matrix is @b not considered
* (operation is always applied to both triangles). */
void axpy(Double c, const_MatrixSliceRef A, MatrixSliceRef B);

/** @brief Matrix matrix multiplication: C += c * A * B.
* The @a type of matrices are considered. */
void matMult(Double c, const_MatrixSliceRef A, const_MatrixSliceRef B, MatrixSliceRef C);

/** @brief rank k update (accumulate normal equations).
* @a N must be SYMMETRIC.
* @f[ N += c \cdot (A^TA) @f] */
void rankKUpdate(Double c, const_MatrixSliceRef A, MatrixSliceRef N);

/** @brief rank 2k update.
* @a N must be SYMMETRIC.
* @f[ N += c \cdot (A^TB+B^TA) @f] */
void rank2KUpdate(Double c, const_MatrixSliceRef A, const_MatrixSliceRef B, MatrixSliceRef N);

/** @brief Inverte a matrix. */
void inverse(MatrixSliceRef A);

/** @brief Eigen value decomposition.
* @f[ N = Q\Lambda Q^T @f]
* @param[in,out]  A input: symmetric Matrix, output matrix of eigen vectors
* @param computeEigenVectors if true eigen vectors are stored in A
* @return eigen values in ascending order */
Vector eigenValueDecomposition(MatrixSliceRef A, Bool computeEigenVectors=TRUE);

/** @brief Eigen value decomposition.
* @param[in]  A input: general N-by-N nonsymmetric Matrix
* @param[out] VL N-by-N, left eigenvectors of A, one per column
* @param[out] VR N-by-N, right eigenvectors of A, one per column
* @param computeEigenVectors if true VL and VR are computed, untouched otherwise
* @return eigen values in ascending order (first column: real part, second column: imaginary part) */
Matrix eigenValueDecomposition(MatrixSliceRef A, Matrix &VL, Matrix &VR, Bool computeEigenVectors=TRUE);

/** @brief singular value decomposition (SVD).
* A = U*diag(S)*Vt.
* @param[in,out]  A input: (n x m) Matrix, output: content is destroyed.
* @param[out] U the left singular vectors (n x s)
* @param[out] Vt the right singular vectors (s x n) with s=min(n,m)
* @param computeSingularVectors if true U and Vt are computed, untouched otherwise
* @return singular values in descending order */
Vector singularValueDecomposition(MatrixSliceRef A, Matrix &U, Matrix &Vt, Bool computeSingularVectors=TRUE);

/** @brief pseudo inverse of matrix A using SVD/EVD
* @param[in] A input: Matrix content is copied to temporary matrix
* @param[in] rcond singular values with an absolute value smaller than max(S)*rcond are assumed to be zero
* @return Moore-Penrose inverse of A */
Matrix pseudoInverse(const_MatrixSliceRef A, Double rcond = 1e-15);

/** @brief matrix square root of a positive semi-definite matrix which fulfills B*B = A and B^T*B = A
* @param[in] A input: Matrix content is copied to temporary matrix
* @param[in] rcond Eigenvalues with an absolute value smaller than max(S)*rcond are assumed to be zero
* @return matrix square root of A */
Matrix matrixSquareRoot(const_MatrixSliceRef A, Double rcond = 1e-15);

/** @brief inverse of matrix square root of a positive semi-definite matrix which fulfills B*B = inv(A) and B^T*B = inv(A)
* @param[in] A input: Matrix content is copied to temporary matrix
* @param[in] rcond Eigenvalues with an absolute value smaller than max(S)*rcond are assumed to be zero
* @return matrix square root of A */
Matrix matrixSquareRootInverse(const_MatrixSliceRef A, Double rcond = 1e-15);

/***********************************************/

/** @brief Solve system of equations for mutliple right hand sides.
* @f[ B := N^{-1} B @f]
* - if @a N is GENERAL content is destroyed.
* - if @a N is SYMMETRIC content will be replaced by the cholesky decomposition N=W'W.
* - if @a N is TRIANGULAR content is unchanged. */
void solveInPlace(MatrixSliceRef N, MatrixSliceRef B);

/** @brief Solve system of equations for mutliple right hand sides.
* @f[ B := N^{-1} B @f]
* - if @a N is GENERAL content is destroyed.
* - if @a N is SYMMETRIC content will be replaced by the cholesky decomposition N=W'W.
* - if @a N is TRIANGULAR content is unchanged. */
Matrix solve(MatrixSliceRef N, const_MatrixSliceRef B);

/** @brief Least squares adjustment with multiple right hand sides.
* Using QR decomposition. Content of A will be destroyed. Residuals are returned in l.
* @f[ x = (A^TA)^{-1}A^T l @f] */
Matrix leastSquares(MatrixSliceRef A, MatrixSliceRef l);

/** @brief Observations l reduced by least suqares fit.
* Using QR decomposition.
* @f[ l := (I-A*(A^TA)^{-1}A^T) l @f] */
void reduceLeastSquaresFit(const_MatrixSliceRef A, MatrixSliceRef l);

/** @brief Eliminate additional parameters from observation equations.
* Additional parameters are defined by the desigmatrix @a B.
* Using @a B=QR decomposition.
* \f[ \bar{A}^T \bar{A} = A^T (I - B (B^T B)^{-1} B^T) A \f] */
void eliminationParameter(MatrixSliceRef B, const std::vector<std::reference_wrapper<Matrix>> &listA);

/** @brief Eliminate additional parameters from observation equations.
* Additional parameters are defined by the desigmatrix @a B.
* Using QR decomposition.
* \f[ \bar{A}^T \bar{A} = A^T (I - B (B^T B)^{-1} B^T) A \f]
* \f[ \bar{A}^T \bar{l} = A^T (I - B (B^T B)^{-1} B^T) l \f] */
inline void eliminationParameter(MatrixSliceRef B, Matrix &A, Matrix &l) {eliminationParameter(B, {A, l});}

/***********************************************/

/** @brief Cholesky decomposition of a SYMMETRIC matrix.
* @param[in,out] A Input: positive definite matrix, Output: upper triangular matrix W (@f$ A = W^TW @f$). */
void cholesky(MatrixSliceRef A);

/** @brief Cholesky decomposition with pivoting of a SYMMETRIC matrix.
* @param[in,out] A Input: positive semidefinite matrix, Output: upper triangular matrix W (@f$ P^TAP = W^TW @f$).
* @param[out] rank Estimated rank of the matrix.
* @param tolerance The algorithm terminates at the (K-1)st step if the pivot <= tolerance.
* @return pivoting vector (@f$ P(piv[i],i)=1 @f$) */
std::vector<UInt> choleskyPivoting(MatrixSliceRef A, UInt &rank, Double tolerance=0);

/** @brief Inverte a matrix, if the cholesky decomposition of the matrix is given. */
void cholesky2Inverse(MatrixSliceRef W);

/** @brief Computes the product WW'. */
void choleskyProduct(MatrixSliceRef W);

/** @brief Triangular matrix - matrix multiplication: C := c * W * C.
* W must be a TRIANGULAR matrix. */
void triangularMult(Double c, const_MatrixSliceRef W, MatrixSliceRef C);

/** @brief Inverse triangular matrix - matrix multiplication: C := c * W^-1 * C.
* W must be a TRIANGULAR matrix. */
void triangularSolve(Double c, const_MatrixSliceRef W, MatrixSliceRef C);

/***********************************************/

/** @brief QR decomposition.
* Computes a QR factorization of a real M-by-N matrix A.
* The matrix Q is represented as a product of elementary reflectors
* Q = H(0) H(1) . . . H(k), where k = min(m,n)-1.
* Each H(i) has the form H(i) = I - tau(i) * v(i) * v(i)'
* where tau is a real scalar, and v is stored on exit in A(i+1:m,i).
* @param[in,out] A On exit, the elements on and above the diagonal of the array contain
* the min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular if m >= n);
* the elements below the diagonal, with the vector tau, represent the orthogonal
* matrix Q as a product of min(m,n) elementary reflectors.
* @return  tau: The scalar factors of the elementary reflectors. */
Vector QR_decomposition(MatrixSliceRef A);

/** @brief Q Multiplication.
* Computes A := Q*A.
* The matrix Q is represented as a product of elementary reflectors
* Q = H(0) H(1) . . . H(k), where k = min(m,n)-1.
* Each H(i) has the form H(i) = I - tau(i) * v(i) * v(i)'
* where tau is a real scalar, and v is stored in B(i+1:m,i).
* @see QR_decomposition */
void QMult(const_MatrixSliceRef B, const Vector &tau, MatrixSliceRef A);

/** @brief Q' Multiplication.
* Computes A := Q'*A.
* The matrix Q is represented as a product of elementary reflectors
* Q = H(0) H(1) . . . H(k), where k = min(m,n)-1.
* Each H(i) has the form H(i) = I - tau(i) * v(i) * v(i)'
* where tau is a real scalar, and v is stored in B(i+1:m,i).
* @see QR_decomposition */
void QTransMult(const_MatrixSliceRef B, const Vector &tau, MatrixSliceRef A);

/** @brief Generate Q matrix.
* Generates  real matrix Q with orthonormal columns,
* The matrix Q is represented at input as a product of elementary reflectors
* Q = H(0) H(1) . . . H(k), where k = min(m,n)-1.
* Each H(i) has the form H(i) = I - tau(i) * v(i) * v(i)'
* where tau is a real scalar, and v is stored in A(i+1:m,i).
* @see QR_decomposition */
void generateQ(MatrixSliceRef A, const Vector &tau);

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/** @brief Internal class.
* Used for memory management of matrices.
* Memory is only copied, if needed (copy on write).
* This means multiple matrices share the same memory,
* as long as the elements are not changed. */
class MatrixBase
{
public:
  explicit MatrixBase(UInt size);  //!< Constructor

  /// count of elements in field.
  UInt size() const {return _size;}

  /// Readonly access to field.
  inline const Double *const_field() const {return ptr.get();}

  /** Writable access to field.
  * Maybe memory must be copied. */
  inline Double *field();

private:
  UInt _size;
  std::shared_ptr<Double> ptr;
};

/// @} group matrix

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Double *MatrixBase::field()
{
  if(ptr.use_count()>1) // not threat save
  {
    auto ptr2 = std::shared_ptr<Double>(new Double[size()], std::default_delete<Double[]>());
    std::copy_n(ptr.get(), size(), ptr2.get());
    std::swap(ptr, ptr2);
  }
  return ptr.get();
}

/***********************************************/

inline Matrix const_MatrixSlice::operator+(Double c)                   const {return Matrix(*this) +=  c;}
inline Matrix const_MatrixSlice::operator-(Double c)                   const {return Matrix(*this) -=  c;}
inline Matrix const_MatrixSlice::operator*(Double c)                   const {return Matrix(*this) *=  c;}
inline Matrix const_MatrixSlice::operator/(Double c)                   const {return Matrix(*this) /=  c;}
inline Matrix const_MatrixSlice::operator-()                           const {return Matrix(*this) *= -1;}
inline Matrix const_MatrixSlice::operator+(const const_MatrixSlice &x) const {return Matrix(*this) +=  x;}
inline Matrix const_MatrixSlice::operator-(const const_MatrixSlice &x) const {return Matrix(*this) -=  x;}

/***********************************************/

inline const Double *const_MatrixSlice::field() const
{
  if(!base)
    throw(Exception("const_MatrixSlice.field: Null-Pointer"));
  return base->const_field()+_start;
}

/***********************************************/

inline Double *MatrixSlice::field() const
{
  if(!base)
    throw(Exception("MatrixSlice.field: Null-Pointer"));
  return base->field()+_start;
}

/***********************************************/

// readonly
inline const Double &const_MatrixSlice::operator () (UInt row, UInt  column) const
{
  if((row>=rows()) || (column>=columns()))
    throw(Exception("Access to matrix element ("s+row%"%i x "s+column%"%i) out of range, matrix size ("s+rows()%"%i x "s+columns()%"%i)"s));
  return *(base->const_field() + (_start + (_rowMajorOrder ? (column+row*_ld) : (row+column*_ld))));
}

/***********************************************/

// writable
inline Double &MatrixSlice::operator () (UInt row, UInt  column)  const
{
  if((row>=rows()) || (column>=columns()))
    throw(Exception("Access to matrix element ("s+row%"%i x "s+column%"%i) out of range, matrix size ("s+rows()%"%i x "s+columns()%"%i)"s));
  return *(base->field() + (_start + (_rowMajorOrder ? (column+row*_ld) : (row+column*_ld))));
}

/***********************************************/

inline const_MatrixSlice const_MatrixSlice::trans() const
{
  const_MatrixSlice x(*this);
  x._rows          = _columns;
  x._columns       = _rows;
  x._uplo          = (_uplo==UPPER) ? LOWER : UPPER;
  x._rowMajorOrder = !_rowMajorOrder;
  return x;
}

/***********************************************/

inline MatrixSlice MatrixSlice::trans() const
{
  MatrixSlice x(*this);
  x._rows          = _columns;
  x._columns       = _rows;
  x._uplo          = (_uplo==UPPER) ? LOWER : UPPER;
  x._rowMajorOrder = !_rowMajorOrder;
  return x;
}

/***********************************************/

inline Matrix operator- (Double c, const_MatrixSliceRef A)
{
  Matrix C(A);
  for(UInt i=0; i<C.rows(); i++)
    for(UInt j=0; j<C.columns(); j++)
      C(i,j) = c-C(i,j);
  return C;
}

/***********************************************/

inline Matrix operator/ (Double c, const_MatrixSliceRef A)
{
  Matrix C(A);
  for(UInt i=0; i<C.rows(); i++)
    for(UInt j=0; j<C.columns(); j++)
      C(i,j) = c/C(i,j);
  return C;
}

/***********************************************/

inline Matrix operator* (const_MatrixSliceRef A, const_MatrixSliceRef B)
{
  Matrix C(A.rows(),B.columns());
  matMult(1.0, A, B, C);
  return C;
}

//***********************************************/

inline Matrix solve(MatrixSliceRef N, const_MatrixSliceRef B)
{
  Matrix X = B;
  solveInPlace(N,X);
  return X;
}

//***********************************************/

#endif
