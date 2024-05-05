/***********************************************/
/**
* @file gnssLambda.cpp
*
* @brief LAMBDA - Least-squares AMBiguity Decorrelation Adjustment.
*
* @author Torsten Mayer-Guerr
* @date 2013-06-24
*
*/
/***********************************************/

#include "base/import.h"
#include "inputOutput/logging.h"
#include "gnssLambda.h"

/***********************************************/

Matrix GnssLambda::phaseDecorrelation(const std::vector<GnssType> &types, Double wavelengthFactor)
{
  try
  {
    if(types.size() == 0)
      return Matrix();
    if(types.size() == 1)
      return Matrix(1, 1, wavelengthFactor*types.at(0).wavelength());

    // design matrix
    // assume for every phase observation an additional range observation
    // parameters: range, TEC, ambiguities
    const UInt dim = types.size();
    Matrix A(2*dim, 2+dim);
    for(UInt i=0; i<dim; i++)
    {
      // phase observations:
      A(i, 0)   = 1.;                                                      // range
      A(i, 1)   = types.at(i).ionosphericFactor();                         // TEC
      A(i, 2+i) = wavelengthFactor*types.at(i).wavelength(); // ambiguity
      // range observations (100 times less accurate):
      A(i+dim, 0) = 1./100.;       // range
      A(i+dim, 1) = -A(i, 1)/100.; // TEC
    }

    // solve & deccorelate
    QR_decomposition(A);
    inverse(A.slice(2, 2, dim, dim));
    Transformation Z(dim);
    choleskyTransform(A.slice(2, 2, dim, dim), Z, FALSE/*timing*/);

    // Transformation matrix cycles -> meter
    Matrix T = Z.transformBack(identityMatrix(types.size()));
    for(UInt idType=0; idType<types.size(); idType++)
      T.row(idType) *= wavelengthFactor*types.at(idType).wavelength();

    return T;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssLambda::choleskyReversePivot(Matrix &N, Transformation &Z, UInt index0Z, Bool timing)
{
  try
  {
    Vector tmp(N.rows());
    Log::Timer timer(N.rows(), 1, timing);
    for(UInt i=0; i<N.rows(); i++)
    {
      timer.loopStep(i);

      // find minimium
      UInt   k    = i;
      Double minN = N(i,i)-tmp(i);
      for(UInt j=i+1; j<N.rows(); j++)
        if(N(j,j)-tmp(j) < minN)
        {
          k    = j;
          minN = N(k,k)-tmp(k);
        }
      // swap
      if(i != k)
      {
        Z.swap(index0Z+i, index0Z+k);
        std::swap(tmp(i), tmp(k));
        std::swap(N(i,i), N(k,k));
        if(i>0)          swap(N.slice(0, i,i, 1),               N.slice(0,   k,   i,     1));
        if(k<N.rows()-1) swap(N.slice(i, k+1, 1, N.rows()-k-1), N.slice(k,   k+1, 1,     N.rows()-k-1));
        if(k>i+1)        swap(N.slice(i, i+1, 1, k-i-1),        N.slice(i+1, k,   k-i-1, 1).trans());
      }
      // cholesky
      N(i,i) -= tmp(i);
      N(i,i)  = std::sqrt(N(i,i));
      if(i+1<N.rows())
      {
        if(i > 0)
          matMult(-1., N.slice(0, i, i, 1).trans(), N.slice(0, i+1, i, N.rows()-1-i), N.slice(i, i+1, 1, N.rows()-1-i));
        N.slice(i,i+1,1,N.rows()-1-i) *= 1./N(i,i);
        for(UInt k=i+1; k<N.rows(); k++)
          tmp(k) += std::pow(N(i,k), 2);
      }
    }
    timer.loopEnd();
    N.setType(Matrix::TRIANGULAR);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// reduce element (i, k) (reduce row k from i)
inline Bool GnssLambda::choleskyReduce(UInt i, UInt k, MatrixSliceRef W, Transformation &Z)
{
  try
  {
    const Double alpha = std::round(W(i,k));
    if(alpha == 0.)
      return FALSE;
    axpy(-alpha, W.slice(k,k,1,W.rows()-k), W.slice(i,k,1,W.rows()-k));
    Z.reduce(alpha, k, i);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssLambda::choleskyTransform(MatrixSliceRef W, Transformation &Z, Bool timing)
{
  try
  {
    const UInt dim = W.rows();

    // Covariance = W D W^T
    // D: diagonal matrix, W: upper triangular
    std::vector<Double> d(dim);
    for(UInt i=0; i<dim; i++)
      d[i] = W(i,i)*W(i,i);
    for(UInt i=0; i<dim; i++)
      W.column(i) *= 1./W(i,i);

    // decorrelate
    for(UInt i=dim; i-->0;)
      for(UInt k=i+1; k<dim; k++)
        choleskyReduce(i, k, W, Z);

    std::vector<Double> delta(dim-1);
    for(UInt i=0; i<dim-1; i++)
      delta[i] = d[i] + std::pow(W(i,i+1), 2) * d[i+1];

    std::multimap<Double, UInt> ratio;
    for(UInt i=0; i<dim-1; i++)
      if(delta[i] < d[i+1])
        ratio.insert(std::make_pair(delta[i]/d[i+1], i));

    Log::Timer timer(1, 1, timing);
    UInt iter = 0;
    while(!ratio.empty())
    {
      timer.loopStep(iter);
      iter++;

      // find maximum change
      const UInt k = ratio.begin()->second;

      // erase ratios that will be changed
      ratio.erase(ratio.begin());
      if((k > 0) && (delta[k-1] < d[k]))
      {
        auto iterpair = ratio.equal_range(delta[k-1]/d[k]);
        ratio.erase(std::find_if(iterpair.first, iterpair.second, [k](const std::pair<Double, UInt> &x){return x.second == k-1;}));
      }
      if((k+1 < dim-1) && (delta[k+1] < d[k+2]))
      {
        auto iterpair = ratio.equal_range(delta[k+1]/d[k+2]);
        ratio.erase(std::find_if(iterpair.first, iterpair.second, [k](const std::pair<Double, UInt> &x){return x.second == k+1;}));
      }

      // swap k and k+1
      const Double eta    = d[k]/delta[k];
      const Double lambda = d[k+1]*W(k,k+1)/delta[k];
      d[k]   = eta*d[k+1];
      d[k+1] = delta[k];

      if(k>0)
        copy(W.slice(0,k,k,2)*Matrix({{-W(k,k+1),  eta}, {1., lambda}}), W.slice(0,k,k,2)); // upper columns
      W(k,k+1) = lambda;
      if(k+2<W.rows())
        swap(W.slice(k,k+2,1,dim-2-k), W.slice(k+1,k+2,1,dim-2-k)); // swap row k and k+1
      Z.swap(k, k+1); // store applied transformation

      // decorrelate, update delta and ratio
      for(UInt i=std::min(k+2,dim-1); i-->std::max(k,UInt(1))-1;)
      {
        // mathematically only choleskyReduce(i, i+1, W, Z) is needed
        if(choleskyReduce(i, i+1, W, Z))
          for(UInt k=i+2; k<dim; k++)
            choleskyReduce(i, k, W, Z);
        delta[i] = d[i] + std::pow(W(i,i+1), 2) * d[i+1];
        if(delta[i] < d[i+1])
          ratio.insert(std::make_pair(delta[i]/d[i+1], i));
      }
    }
    timer.loopEnd();

    // decorrelate rest of the triangle
    for(UInt i=dim; i-->0;)
      for(UInt k=i+1; k<dim; k++)
        choleskyReduce(i, k, W, Z);

    return d;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssLambda::searchInteger(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d,
                               UInt maxSearchSteps, Vector &solution, Double &minNorm)
{
  try
  {
    const UInt dim = W.rows();
    std::vector<Double> xInt(dim, 0);
    std::vector<Double> norm(dim, 0);
    std::vector<Double> step(dim, 0);
    std::vector<Double> dx  (dim, 0);

    // store xBar in lower triangle of W
    // xBar vector of step i starts at W(dim-1-i, dim-1-i), reverse order
    MatrixSlice xBar(W);
    copy(xFloat, xBar.column(0));

    UInt i  = dim-1;
    xInt[i] = std::round(xBar(dim-1, dim-1-i));
    dx[i]   = xInt[i] - xBar(dim-1, dim-1-i);
    step[i] = (dx[i] < 0) ? 1 : -1;
    minNorm = 1e99;
    Double newNorm = dx[i]*dx[i]/d(i,0);

    UInt iter=0;
    for(; iter<maxSearchSteps; iter++)
    {
      // move down
      const UInt lastValidXBar = i;
      while((newNorm < minNorm) && (i-- > 0))
      {
        // compute new xBar
        for(UInt k=lastValidXBar; k-->i;)
          xBar(i+dim-1-k, dim-1-k) = xBar(i+dim-2-k, dim-2-k) + dx[k+1] * W(i,k+1);

        xInt[i]  = std::round(xBar(dim-1, dim-1-i));
        dx[i]    = xInt[i] - xBar(dim-1, dim-1-i);
        step[i]  = (dx[i] < 0) ? 1 : -1;
        norm[i]  = newNorm;
        newNorm += dx[i]*dx[i]/d(i,0);
      }

      // new solution?
      if(newNorm < minNorm)
      {
        i = 0;
        solution = Vector(xInt);
        minNorm  = newNorm;
        newNorm  = 2e99;
      }

      // move up
      while((newNorm >= minNorm) && (i++ < dim-1))
      {
        xInt[i] += step[i];
        dx[i]   += step[i];
        step[i]  = (step[i]>0) ? (-step[i]-1) : (-step[i]+1); // zig-zag search
        newNorm  = norm[i] + dx[i]*dx[i]/d(i,0);
      }

      if(newNorm >= minNorm)
        break;
    } // for(iter)

    // restore diagonal
    for(UInt i=0; i<dim; i++)
      W(i,i) = 1.;

    return (iter < maxSearchSteps);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssLambda::searchIntegerBlocked(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d,
                                      Double sigmaMaxResolve, UInt searchBlockSize, UInt maxSearchSteps, IncompleteAction incompleteAction, Bool timing,
                                      Vector &isNotFixed, Double &sigma, Matrix &solutionSteps)
{
  try
  {
    const UInt dim = W.rows();
    Vector xInt(dim);
    isNotFixed     = Vector(dim, 1.);
    solutionSteps  = Matrix();
    sigma          = 0.;

    // find values to be solved
    UInt minIndex = dim;
    while((minIndex>0) && (d(minIndex-1,0) < sigmaMaxResolve*sigmaMaxResolve))
      minIndex--;
    if(minIndex >= dim)
      return xInt;
    UInt idxSolved = dim;

    // float solution
    std::vector<Vector> solutions(1, d);
    auto storeSolution = [&](Vector xBar, UInt idxSolved)
    {
      for(UInt i=dim; i-->idxSolved;)
        axpy(-xBar(i), W.slice(0,i,i,1), xBar.row(0,i));
      for(UInt i=0; i<idxSolved; i++)
        xBar(i) = std::fabs(std::remainder(xBar(i), 1));
      for(UInt i=idxSolved; i<dim; i++)
        xBar(i) = -std::fabs(xBar(i));
      solutions.push_back(xBar);
    };
    storeSolution(xFloat, idxSolved);

    UInt defaultBlockSize = std::min(dim-minIndex, searchBlockSize);
    UInt blockSize     = defaultBlockSize;
    UInt blockStart    = dim-blockSize;
    UInt blockStartOld = dim;
    UInt iter=0;
    Log::Timer timer((dim-minIndex)/(defaultBlockSize/2), 1, timing);
    for(;;)
    {
      timer.loopStep(iter++);

      // compute xBar for found solution so far
      Vector xBar = xFloat-xInt;
      for(UInt i=dim; i-->blockStart+blockSize;)
        axpy(-xBar(i), W.slice(0,i,i,1), xBar.row(0,i));

      Vector dxInt;
      Double ePe;
      Bool completed = searchInteger(xBar.row(blockStart, blockSize),
                                     W.slice(blockStart, blockStart, blockSize, blockSize),
                                     d.row(blockStart, blockSize), maxSearchSteps, dxInt, ePe);
      axpy(1., dxInt, xInt.row(blockStart, blockSize));
      copy(Vector(blockSize, 0), isNotFixed.row(blockStart, blockSize));
      idxSolved = blockStart;

      if(!completed)
      {
        if(incompleteAction == IncompleteAction::STOP)
        {
          if(timing) logWarning<<"searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): cannot find a solution -> stop searching"<<Log::endl;
          idxSolved = blockStart+blockSize;
          break;
        }
        else if(incompleteAction == IncompleteAction::SHRINKBLOCKSIZE)
        {
          if(timing) logWarning<<"searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): cannot find a solution -> restart with smaller block size"<<Log::endl;
          defaultBlockSize = (defaultBlockSize+1)/2;
          blockStart = std::min(blockStartOld, blockStart+(blockSize+1)/2);
          blockSize -= (blockSize+1)/2;
          iter--;
          continue;
        }
        else if(incompleteAction == IncompleteAction::EXCEPTION)
        {
          throw(Exception("Ambiguity resolution failed."));
        }
        // incompleteAction == IncompleteAction::IGNORE)
      }

      storeSolution(xFloat-xInt, idxSolved);

      // check consistence with solution of old block, if not restart with larger block
      if(quadsum(dxInt.row(blockStartOld-blockStart, blockStart+blockSize-blockStartOld)))
      {
        if(timing) logWarning<<"searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): not consistent -> restart"<<Log::endl;
        blockStartOld = blockStart+blockSize;
        blockSize    += std::min(defaultBlockSize/2, dim-blockStart-blockSize);
        iter--;
        continue;
      }

      blockStartOld = blockStart;
      blockSize     = defaultBlockSize;
      if(blockStart == minIndex) // at beginning?
        break;
      blockStart = std::max(blockStart, minIndex+defaultBlockSize/2) - defaultBlockSize/2;
    } // for(blocks)
    timer.loopEnd();

    xInt.row(0, idxSolved).fill(0);
    isNotFixed.row(0, idxSolved).fill(1);

    // compute norm
    Vector xBar = xFloat - xInt;
    triangularSolve(1., W, xBar);
    sigma = 0;
    for(UInt i=idxSolved; i<dim; i++)
      sigma += xBar(i)*xBar(i)/d(i,0);
    sigma = std::sqrt(sigma/(dim-idxSolved));

    // copy solutions into one matrix
    solutionSteps = Matrix(dim, solutions.size());
    for(UInt i=0; i<solutions.size(); i++)
      copy(solutions.at(i), solutionSteps.column(i));

    return xInt;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GnssLambda::Transformation::Transformation(UInt dim) : columnForward(dim), columnBackward(dim)
{
  try
  {
    // unity matrix in both directions
    for(UInt i=0; i<dim; i++)
    {
      columnForward[i].insert (std::make_pair(i, 1.0));
      columnBackward[i].insert(std::make_pair(i, 1.0));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// row(k) -= alpha * row(i)
inline void GnssLambda::Transformation::reduce(Double alpha, UInt i, UInt k)
{
  try
  {
    if(alpha == 0.)
      return;

    // --- lambda ---
    auto reduce = [](Double alpha, std::map<UInt, Double> &column1, std::map<UInt, Double> &column2)
    {
      for(auto &col : column1)
      {
        auto element = column2.find(col.first);
        if(element != column2.end())
        {
          element->second -= alpha * col.second;
          if(std::fabs(element->second) < 1e-9)
            column2.erase(element);
        }
        else
          column2.insert(std::make_pair(col.first, -alpha * col.second));
      }
    };
    // -----------

    reduce(+alpha, columnForward[i],  columnForward[k]);
    reduce(-alpha, columnBackward[k], columnBackward[i]);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssLambda::Transformation::swap(UInt i, UInt k)
{
  std::swap(columnForward.at(i),  columnForward.at(k));
  std::swap(columnBackward.at(i), columnBackward.at(k));
}

/***********************************************/

Matrix GnssLambda::Transformation::transform(const_MatrixSliceRef x) const
{
  try
  {
    if(x.rows() != columnForward.size())
      throw(Exception("Dimension error"));

    Matrix y(x.rows(), x.columns());
    for(UInt i=0; i<columnForward.size(); i++)
      for(const auto &col : columnForward[i])
        axpy(col.second, x.row(col.first), y.row(i));
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GnssLambda::Transformation::transformBack(const_MatrixSliceRef x) const
{
  try
  {
    Matrix y(x.rows(), x.columns());
    for(UInt i=0; i<columnBackward.size(); i++)
      for(const auto &col : columnBackward[i])
        axpy(col.second, x.row(i), y.row(col.first));
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GnssLambda::Transformation::distributeBack(const_MatrixSliceRef x) const
{
  try
  {
    Matrix y(x.rows(), x.columns());
    for(UInt i=0; i<columnBackward.size(); i++)
      for(UInt k=0; k<x.columns(); k++)
        if(x(i, k))
          for(const auto &col : columnBackward[i])
            y(col.first, k) = 1.;
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
