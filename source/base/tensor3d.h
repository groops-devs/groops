/***********************************************/
/**
* @file tensor3d.h
*
* @brief tensor in 3d space.
*
* @author Torsten Mayer-Guerr
* @date 2001-01-05
*
*/
/***********************************************/

#ifndef __GROOPS_TENSOR3D__
#define __GROOPS_TENSOR3D__

#include "base/importStd.h"
#include "base/matrix.h"

/** @addtogroup vector3dGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Tensor in 3d space.
* This is the result of a dyadic product of two @a Vector3d.
* It is a symmetric matrix in cartesian coordinates.
* For representation of e.g. gravity gradients or covariances of positions.
* Can be rotated with @a Rotary3d.
* Can be added, subtracted and scaled with a Double. */
class Tensor3d
{
    std::array<Double,6> field;  // xx,yy,zz,xy,xz,yz

public:
  Tensor3d() : field{0.,0.,0.,0.,0.,0.} {} //!< Default constructor (zero matrix).

  Tensor3d(const_MatrixSliceRef x); //!< Default constructor from symmetric 3x3 matrix

  /** @brief Read/write cartesian coordinates.
  * Read/write cartesian coordinates in the form:
  * @code
  *  Tensor3d t;
  *  t.xx() = 0.5;
  *  Double xx = t.xx();
  * @endcode */
  Double xx() const {return field[0];}
  Double yy() const {return field[1];}  //!< @copydoc xx
  Double zz() const {return field[2];}  //!< @copydoc xx
  Double xy() const {return field[3];}  //!< @copydoc xx
  Double xz() const {return field[4];}  //!< @copydoc xx
  Double yz() const {return field[5];}  //!< @copydoc xx

  Double &xx() {return field[0];} //!< @copydoc xx
  Double &yy() {return field[1];} //!< @copydoc xx
  Double &zz() {return field[2];} //!< @copydoc xx
  Double &xy() {return field[3];} //!< @copydoc xx
  Double &xz() {return field[4];} //!< @copydoc xx
  Double &yz() {return field[5];} //!< @copydoc xx

  /// Cast to matrix.
  Matrix matrix() const;

  Tensor3d &operator+= (Tensor3d const &b);
  Tensor3d &operator-= (Tensor3d const &b);
  Tensor3d &operator*= (Double  c);
};

/***********************************************/

inline Tensor3d operator- (const Tensor3d &t)                      {return Tensor3d(t)  *= -1;}
inline Tensor3d operator+ (const Tensor3d &t1, const Tensor3d &t2) {return Tensor3d(t1) += t2;}
inline Tensor3d operator- (const Tensor3d &t1, const Tensor3d &t2) {return Tensor3d(t1) -= t2;}
inline Tensor3d operator* (Double c, const Tensor3d &t)            {return Tensor3d(t)  *= c;}
inline Tensor3d operator* (const Tensor3d &t, Double c)            {return Tensor3d(t)  *= c;}

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Tensor3d &Tensor3d::operator+= (Tensor3d const &b)  {for(UInt i=0; i<6; i++) field[i] += b.field[i]; return *this;}
inline Tensor3d &Tensor3d::operator-= (Tensor3d const &b)  {for(UInt i=0; i<6; i++) field[i] -= b.field[i]; return *this;}
inline Tensor3d &Tensor3d::operator*= (Double  c)          {for(UInt i=0; i<6; i++) field[i] *= c;          return *this;}

/***********************************************/

inline Tensor3d::Tensor3d(const_MatrixSliceRef T)
{
  try
  {
    if((T.getType()!=Matrix::SYMMETRIC) || (T.rows()!=3) || !T.isUpper())
      throw(Exception("Matrix T("s+T.rows()%"%i x "s+T.columns()%"%i) must be 3x3 upper SYMMETRIC"s));
    field[0] = T(0,0);
    field[1] = T(1,1);
    field[2] = T(2,2);
    field[3] = T(0,1);
    field[4] = T(0,2);
    field[5] = T(1,2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix Tensor3d::matrix() const
{
  Matrix T(3, Matrix::SYMMETRIC);
  T(0,0) = field[0];
  T(1,1) = field[1];
  T(2,2) = field[2];
  T(0,1) = field[3];
  T(0,2) = field[4];
  T(1,2) = field[5];
  fillSymmetric(T);
  return T;
}

/*************************************************/

#endif /* __GROOPS_TENSOR3D__ */


