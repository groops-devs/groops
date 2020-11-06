/***********************************************/
/**
* @file kalmanProcessing.cpp
*
* @brief Miscellaneous functions for Kalman filter applications
*
* @author Andreas Kvas
* @date 2016-12-26
*
*/
/***********************************************/

#define DOCSTRING_AutoregressiveModelSequence

#include "kalmanProcessing.h"
#include "config/config.h"
#include "config/configRegister.h"
#include "files/fileMatrix.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(AutoregressiveModelSequence, "autoregressiveModelSequenceType")
GROOPS_READCONFIG_CLASS(AutoregressiveModelSequence, "autoregressiveModelSequenceType")

/***********************************************/

void AutoregressiveModel::computeNormalEquation()
{
  Matrix N_epoch(_model.columns(), Matrix::SYMMETRIC, Matrix::UPPER);
  rankKUpdate(1.0, _model, N_epoch);
  fillSymmetric(N_epoch);

  const UInt order = this->order();
  const UInt dim = this->dimension();

  N.resize(order+1, std::vector<Matrix>(order+1));
  for(UInt r = 0; r < N.size(); r++)
    for(UInt c = r; c < N.size(); c++)
      N.at(r).at(c) = N_epoch.slice(r*dim, c*dim, dim, dim);

  hasNormals = TRUE;
}

/***********************************************/

void AutoregressiveModel::distributedNormalsBlock(UInt epochCount, UInt row, UInt column, Matrix &X)
{
  if(row > column)
    throw(Exception("Process normals only implemented for upper triangular matrices!"));

  if(!hasNormals)
    this->computeNormalEquation();

  for(UInt l = 0; l<epochCount-this->order(); l++)
    if( (row>=l) && (column<=this->order()+l) )
      axpy(1.0, N.at(row-l).at(column-l), X);
}

/***********************************************/

Matrix AutoregressiveModel::distributedNormalsBlock(UInt epochCount, UInt row, UInt column)
{
  if(row > column)
    throw(Exception("Process normals only implemented for upper triangular matrices!"));

  Matrix X;
  if(row == column)
    X = Matrix(this->dimension(), Matrix::SYMMETRIC, Matrix::UPPER);
  else
    X = Matrix(this->dimension(), this->dimension());

  this->distributedNormalsBlock(epochCount, row, column, X);

  return X;
}

/***********************************************/

std::vector<Matrix> AutoregressiveModel::coefficients() const
{
  Matrix W = pseudoInverse(_model.column(order()*dimension(), dimension()));

  std::vector<Matrix> B;
  for(UInt k = 0; k<order(); k++)
    B.insert(B.begin(), -1.0*W*_model.column(k*dimension(), dimension())); // B_k = -sqrt(Q)*A_k

  return B;
}

/***********************************************/

Matrix AutoregressiveModel::whiteNoiseCovariance() const
{
  Matrix W = pseudoInverse(_model.column(order()*dimension(), dimension()));

  Matrix Q(dimension(), Matrix::SYMMETRIC);
  rankKUpdate(1.0, W.trans(), Q); // Q = W W^T

  return Q;
}

/***********************************************/

void AutoregressiveModel::orderOneRepresentation(Matrix &B, Matrix &Q) const
{
  std::vector<Matrix> B_tmp = this->coefficients();
  Matrix Q_tmp = this->whiteNoiseCovariance();

  if(order() == 1)
  {
    B = B_tmp.front();
    Q = Q_tmp;
  }
  else
  {
    B = Matrix(order()*dimension(), order()*dimension());
    for(UInt k = 0; k<order(); k++)
      copy(B_tmp.at(k), B.slice(0, k*dimension(), dimension(), dimension()));

    for(UInt i = 0; i<dimension()*(order()-1); i++)
      B(i+dimension(), i) = 1.0;

    Q = Matrix(order()*dimension(), Matrix::SYMMETRIC);
    copy(Q_tmp, Q.slice(0,0, dimension(), dimension()));
  }
}

/***********************************************/

std::vector< std::pair<UInt, UInt> > AutoregressiveModel::distributedNormalsBlockIndex(UInt blockCount) const
{
  std::vector< std::pair<UInt, UInt> > index;
  for(UInt row = 0; row<blockCount; row++)
  {
    for(UInt column = row; column<std::min(row+this->order()+1, blockCount); column++)
      index.push_back( std::pair<UInt, UInt>(row, column));
  }

  return index;
}

/***********************************************/

void AutoregressiveModel::updatePredictionErrorCovariance(const Matrix& Q_new)
{
  if(Q_new.rows() != this->dimension())
    throw(Exception("Size of new prediction error covariance does not match."));

  Matrix W = pseudoInverse(_model.column(order()*dimension(), dimension()));
  Matrix W_new = matrixSquareRootInverse(Q_new);

  _model = (W_new*W)*_model;
}

/***********************************************/

Vector AutoregressiveModel::decorrelate(const Vector& x) const
{
  if( (x.rows() % dimension()) != 0 )
    throw(Exception("Size of solution vector does not match model dimension."));

  const UInt epochCount = x.rows()/dimension();

  Matrix X((order()+1)*dimension(), epochCount-order());
  for(UInt k = 0; k<epochCount-order(); k++)
    copy(x.row(k*dimension(), X.rows()), X.column(k));

  Matrix U(dimension(), epochCount-order());
  matMult(-1.0, _model, X, U);

  return flatten(U);
}

/***********************************************/

void AutoregressiveModel::rotate(const_MatrixSliceRef F)
{
  if(F.rows() != dimension())
    throw(Exception("Dimensions of rotation matrix and model do not match (<" + F.rows()%"%i"s + "> vs. <"+dimension()%"%i"s+">)."));

  Matrix G = pseudoInverse(_model.column(order()*dimension(), dimension())*F); // G = F^+ * W^0.5

  Matrix Q(F.columns(), Matrix::SYMMETRIC);
  rankKUpdate(1.0, G.trans(), Q);
  Matrix V = matrixSquareRootInverse(Q);

  G = V*G; // apply decorrelation to G, G := Q^-0.5 * F^+ * W^0.5

  Matrix B(F.columns(), F.columns()*(order()+1)); // new model

  for(UInt k = 0; k<order(); k++)
    copy(G*_model.column(k*dimension(), dimension())*F, B.column(k*F.columns(), F.columns()));
  copy(V, B.column(order()*F.columns(), F.columns()));

  _model = B;
}

/***********************************************/
/***********************************************/

AutoregressiveModelSequencePtr AutoregressiveModelSequence::create(Config &config, const std::string &name)
{
  try
  {
    AutoregressiveModelSequencePtr arSequence;
    readConfigSequence(config, name, Config::MUSTSET, "", "");
    arSequence = AutoregressiveModelSequencePtr(new AutoregressiveModelSequence(config));
    endSequence(config);
    return arSequence;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

AutoregressiveModelSequence::AutoregressiveModelSequence(Config &config)
{
  std::vector<FileName> fileNameArModel;
  Double sigma0 = 1.0;

  readConfig(config, "inputfileAutoregressiveModel", fileNameArModel, Config::MUSTSET,  "",    "matrix file containing an AR model");
  readConfig(config, "sigma0",                       sigma0,          Config::DEFAULT,  "1.0", "a-priori sigma for white noise covariance");
  if(isCreateSchema(config)) return;

  std::set<UInt> dimensions;
  for(auto &fileName : fileNameArModel)
  {
    Matrix tmp;
    readFileMatrix(fileName, tmp);
    tmp*=(1.0/sigma0);
    _arModels.push_back(AutoregressiveModel(tmp));

    dimensions.insert(_arModels.back().dimension());
  }

  if(dimensions.size() != 1)
    throw(Exception("Autoregressive models must have same dimension."));

  for(UInt order = 0; order<_arModels.size(); order++)
    if(order != _arModels.at(order).order())
      throw(Exception("Autoregressive models must be given in increasing and consecutive order."));
}

/***********************************************/

void AutoregressiveModelSequence::distributedNormalsBlock(UInt epochCount, UInt row, UInt column, Matrix &X)
{
  _arModels.back().distributedNormalsBlock(epochCount, row, column, X);  // last model constrains all parameters

  for(UInt order = 0; order<maximumOrder(); order++) // other models constrain the first maxOrder-1 blocks
    if(row<=order && column<=order)
      axpy(1.0, _arModels.at(order).normalEquation().at(row).at(column), X);
}

/***********************************************/

Matrix AutoregressiveModelSequence::distributedNormalsBlock(UInt epochCount, UInt row, UInt column)
{
  if(row > column)
    throw(Exception("Process normals only implemented for upper triangular matrices!"));

  Matrix X;
  if(row == column)
    X = Matrix(this->dimension(), Matrix::SYMMETRIC, Matrix::UPPER);
  else
    X = Matrix(this->dimension(), this->dimension());

  this->distributedNormalsBlock(epochCount, row, column, X);

  return X;
}


/***********************************************/

std::vector< std::vector< std::vector<Matrix> > > AutoregressiveModelSequence::normalEquationSequence()
{
  std::vector< std::vector< std::vector<Matrix> > > N;
  for(AutoregressiveModel &model : _arModels)
    N.push_back(model.normalEquation());

  return N;
}

/***********************************************/
