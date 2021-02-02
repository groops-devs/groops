/***********************************************/
/**
* @file normalsShortTimeStaticLongTime.cpp
*
* @brief Normal equations with short and long time gravity variations.
*
* @author Torsten Mayer-Guerr
* @date 2020-11-29
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "classes/observation/observation.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "misc/varianceComponentEstimation.h"
#include "normalsShortTimeStaticLongTime.h"

/***********************************************/

void NormalsShortTimeStaticLongTime::init(ObservationPtr observation, const std::vector<Time> &timesInterval,
                                          UInt defaultBlockSize, Parallel::CommunicatorPtr comm, Bool sortOtherParametersBeforeGravityParameters,
                                          UInt countShortTimeParameters, ParameterSelectorPtr parameterShortTime,
                                          ParametrizationTemporalPtr temporal, ParameterSelectorPtr parameterTemporal)
{
  try
  {
    // parameter names
    // ---------------
    std::vector<ParameterName> paraNameAll;
    observation->parameterName(paraNameAll);

    // examine intervals
    // -----------------
    std::vector<UInt> minInterval(paraNameAll.size(), timesInterval.size());
    std::vector<UInt> maxInterval(paraNameAll.size(), 0);
    std::vector<std::vector<UInt>> indexDesign(timesInterval.size()-1, std::vector<UInt>(paraNameAll.size(), NULLINDEX));

    for(UInt idInterval=0; idInterval+1<timesInterval.size(); idInterval++)
    {
      std::vector<ParameterName> paraNameInterval;
      observation->setInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1));
      observation->parameterName(paraNameInterval);

      auto iter = paraNameAll.begin();
      for(UInt k=0; k<paraNameInterval.size(); k++)
      {
        iter = std::find(iter, paraNameAll.end(), paraNameInterval.at(k));
        const UInt i = std::distance(paraNameAll.begin(), iter);
        indexDesign.at(idInterval).at(i) = k;
        minInterval.at(i) = std::min(minInterval.at(i), idInterval);
        maxInterval.at(i) = std::max(maxInterval.at(i), idInterval);
      }
    }
    observation->setInterval(timesInterval.front(), timesInterval.back());

    // --- lambda -------------------
    auto testSelectedParameters = [&](const std::string &text, std::vector<UInt>::const_iterator begin, std::vector<UInt>::const_iterator end)
    {
      if(begin == end)
        throw(Exception(text+": no parameters selected"));
      if(!std::is_sorted(begin, end) || ((*(end-1)-*begin+1) != static_cast<UInt>(std::distance(begin, end))))
        throw(Exception(text+": selected parameters must be in consecutive order"));
      if(!std::all_of(begin, end, [&](UInt i) {return (minInterval.at(i) == minInterval.at(*begin)) && (maxInterval.at(i) == maxInterval.at(*begin));}))
        throw(Exception(text+": all selected parameters must be in the same time interval"));
    };
    // ------------------------------

    // add short time gravity parameters
    // ---------------------------------
    if(countShortTimeParameters && parameterShortTime)
    {
      std::vector<UInt> index = parameterShortTime->indexVector(paraNameAll);
      if(index.size() % countShortTimeParameters)
        throw(Exception("Number of selected short time parameters disagree with dimension of AR model"));

      // is only one static block? -> copy static parameters for each interval
      if(index.size() == countShortTimeParameters)
      {
        testSelectedParameters("short time variations", index.begin(), index.end());
        paraNameAll.reserve(paraNameAll.size() + index.size() * (timesInterval.size()-1));
        minInterval.reserve(paraNameAll.size() + index.size() * (timesInterval.size()-1));
        maxInterval.reserve(paraNameAll.size() + index.size() * (timesInterval.size()-1));

        for(UInt idInterval=0; idInterval+1<timesInterval.size(); idInterval++)
        {
          blockIndexShortTime.push_back(paraNameAll.size());
          indexDesign.at(idInterval).reserve(paraNameAll.size() + index.size() * (timesInterval.size()-1));
          const ParameterName paraNameInterval("", "", "", timesInterval.at(idInterval), timesInterval.at(idInterval+1));
          for(UInt i=0; i<index.size(); i++)
          {
            paraNameAll.push_back(paraNameAll.at(index.at(i)));
            paraNameAll.back().combine(paraNameInterval);
            for(UInt idInterval2=0; idInterval2+1<timesInterval.size(); idInterval2++)
              indexDesign.at(idInterval2).push_back((idInterval == idInterval2) ? indexDesign.at(idInterval).at(index.at(i)) : NULLINDEX);
            minInterval.push_back(idInterval);
            maxInterval.push_back(idInterval);
          }
        }
      }
      else
        for(UInt i=0; i<index.size(); i+=countShortTimeParameters)
        {
          testSelectedParameters("short time variations", index.begin()+i, index.begin()+(i+countShortTimeParameters));
          blockIndexShortTime.push_back(index.at(i));
        }
    }

    // add long time gravity parameters
    // --------------------------------
    UInt countTemporalParameters = 0;
    if(temporal && parameterTemporal)
    {
      std::vector<UInt> index = parameterTemporal->indexVector(paraNameAll);
      testSelectedParameters("long time variations", index.begin(), index.end());
      countTemporalParameters = index.size();

      blockIndexTemporal = {index.at(0)}; // static
      for(UInt k=0; k<temporal->parameterCount(); k++)
        blockIndexTemporal.push_back(paraNameAll.size()+k*index.size());

      minInterval.resize(paraNameAll.size()+index.size()*temporal->parameterCount(), 0);
      maxInterval.resize(paraNameAll.size()+index.size()*temporal->parameterCount(), timesInterval.size()-2);
      for(UInt idInterval=0; idInterval+1<timesInterval.size(); idInterval++)
        indexDesign.at(idInterval).resize(paraNameAll.size()+index.size()*temporal->parameterCount(), NULLINDEX);

      std::vector<ParameterName> parameterNameBase(index.size());
      for(UInt i=0; i<index.size(); i++)
        parameterNameBase.at(i) = paraNameAll.at(index.at(i));
      temporal->parameterName(parameterNameBase, paraNameAll);

      factorTemporal.resize(1+temporal->parameterCount(), std::vector<Double>(timesInterval.size()-1)); // static + temporal
      for(UInt idInterval=0; idInterval<timesInterval.size()-1; idInterval++)
      {
        const Vector f = temporal->factors(0.5*(timesInterval.at(idInterval)+timesInterval.at(idInterval+1)));
        factorTemporal.at(0).at(idInterval) = 1.; // static
        for(UInt k=0; k<f.rows(); k++)
          factorTemporal.at(1+k).at(idInterval) = f(k);
      }
    }

    // sort parameters: interval parameters first
    // ------------------------------------------
    std::vector<UInt> sortIndex(paraNameAll.size());
    std::iota(sortIndex.begin(), sortIndex.end(), 0);
    blockIndexStatic = std::distance(sortIndex.begin(), std::find_if(sortIndex.begin(), sortIndex.end(), [&](UInt i)
                                                                     {return (minInterval.at(i) == 0) && (maxInterval.at(i) == timesInterval.size()-2);}));
    std::stable_sort(sortIndex.begin(), sortIndex.end(), [&] (UInt i, UInt k)
    {
      if(maxInterval.at(i) != maxInterval.at(k)) return (maxInterval.at(i) < maxInterval.at(k));
      if(minInterval.at(i) != minInterval.at(k)) return (minInterval.at(i) < minInterval.at(k));
      if(sortOtherParametersBeforeGravityParameters)
      {
        if(((observation->gravityParameterCount() <= i) && (i < observation->parameterCount())) && // is state parameter?
          !((observation->gravityParameterCount() <= k) && (k < observation->parameterCount())))   // is not state parameter?
        return TRUE;
      }
      return FALSE;
    });

    parameterNames.clear(); parameterNames.reserve(paraNameAll.size());
    for(UInt i=0; i<sortIndex.size(); i++)
      parameterNames.push_back(paraNameAll.at(sortIndex.at(i)));

    // divide into blocks
    // ------------------
    std::vector<UInt> blockIndex(1, 0);
    while(blockIndex.back() < paraNameAll.size())
    {
      // ---------------------
      auto mustSplit = [&](UInt i, UInt k)
      {
        for(UInt idInterval=0; idInterval+1<timesInterval.size(); idInterval++)
          if((indexDesign.at(idInterval).at(i) != NULLINDEX) || (indexDesign.at(idInterval).at(k) != NULLINDEX))
            if((indexDesign.at(idInterval).at(i) == NULLINDEX) || (indexDesign.at(idInterval).at(k) == NULLINDEX) || (indexDesign.at(idInterval).at(i) != indexDesign.at(idInterval).at(k)-1))
              return TRUE;
        if(std::find_if(blockIndexShortTime.begin(), blockIndexShortTime.end(), [&](UInt idx) {return (idx == k) || (idx+countShortTimeParameters-1 == i);}) != blockIndexShortTime.end())
          return TRUE;
        if(std::find_if(blockIndexTemporal.begin(), blockIndexTemporal.end(), [&](UInt idx) {return (idx == k) || (idx+countTemporalParameters-1 == i);}) != blockIndexTemporal.end())
          return TRUE;
        return FALSE;
      };
      // ---------------------

      blockIndex.push_back(blockIndex.back()+1);
      while((blockIndex.back() < paraNameAll.size()) && !mustSplit(sortIndex.at(blockIndex.back()-1), sortIndex.at(blockIndex.back())))
        blockIndex.back()++;
    }

    // divide blocks according to defaultBlockSize (if not a high frequency block)
    if(defaultBlockSize > 0)
      for(UInt i=0; i<blockIndex.size()-1; i++)
        if((blockIndex.at(i+1)-blockIndex.at(i) > defaultBlockSize*3/2) &&
           (std::find(blockIndexShortTime.begin(), blockIndexShortTime.end(), sortIndex.at(blockIndex.at(i))) == blockIndexShortTime.end()))
          blockIndex.insert(blockIndex.begin()+(i+1), blockIndex.at(i)+defaultBlockSize);

    // Init normal equations
    // ---------------------
    initEmpty(blockIndex, comm);
    n        = Matrix(parameterCount(), observation->rightSideCount());
    lPl      = Vector(observation->rightSideCount());
    obsCount = 0;
    if(countTemporalParameters)
    {
      normalsTemporal.resize(timesInterval.size()-1);
      for(UInt idInterval=0; idInterval<timesInterval.size()-1; idInterval++)
        normalsTemporal.at(idInterval).initEmpty(blockIndex, comm);
      nTemporal.clear();
      nTemporal.resize(timesInterval.size()-1, Matrix(parameterCount(), observation->rightSideCount()));
    }

    // compute correct block indices
    // -----------------------------
    blockIndexStatic = index2block(std::distance(sortIndex.begin(), std::find(sortIndex.begin(), sortIndex.end(), blockIndexStatic)));
    for(UInt &idx : blockIndexShortTime)
      idx = index2block(std::distance(sortIndex.begin(), std::find(sortIndex.begin(), sortIndex.end(), idx)));
    for(UInt &idx : blockIndexTemporal)
      idx = index2block(std::distance(sortIndex.begin(), std::find(sortIndex.begin(), sortIndex.end(), idx)));

    blockCountTemporal = 0;
    if(blockIndexTemporal.size())
      blockCountTemporal = index2block(blockIndex.at(blockIndexTemporal.at(0))+countTemporalParameters) - blockIndexTemporal.at(0);

    // index for each block in design matrix
    // -------------------------------------
    indexN.clear(); indexN.resize(timesInterval.size()-1);
    indexA.clear(); indexA.resize(timesInterval.size()-1);
    for(UInt idBlock=0; idBlock<blockIndex.size()-1; idBlock++)
    {
      const UInt idx = blockIndex.at(idBlock);
      for(UInt idInterval=0; idInterval+1<timesInterval.size(); idInterval++)
        if(indexDesign.at(idInterval).at(sortIndex.at(idx)) != NULLINDEX)
        {
          indexN.at(idInterval).push_back(idBlock);
          indexA.at(idInterval).push_back(indexDesign.at(idInterval).at(sortIndex.at(idx)));
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsShortTimeStaticLongTime::setBlocks(const std::vector<UInt> &arcsInterval)
{
  try
  {
    for(UInt idInterval=0; idInterval+1<arcsInterval.size(); idInterval++)
      if(arcsInterval.at(idInterval+1)-arcsInterval.at(idInterval) > 0)
        for(UInt i=0; i<indexN.at(idInterval).size(); i++)
          for(UInt k=i; k<indexN.at(idInterval).size(); k++)
          {
            const Bool isTemporal1 = blockCountTemporal && (blockIndexTemporal.at(0) <= indexN.at(idInterval).at(i)) && (indexN.at(idInterval).at(i) < blockIndexTemporal.at(0)+blockCountTemporal);
            const Bool isTemporal2 = blockCountTemporal && (blockIndexTemporal.at(0) <= indexN.at(idInterval).at(k)) && (indexN.at(idInterval).at(k) < blockIndexTemporal.at(0)+blockCountTemporal);
            if(isTemporal1 || isTemporal2)
            {
              normalsTemporal.at(idInterval).setBlock(indexN.at(idInterval).at(i), indexN.at(idInterval).at(k));
              normalsTemporal.at(idInterval).N(indexN.at(idInterval).at(i), indexN.at(idInterval).at(k)) = Matrix();
            }
            else
            {
              setBlock(indexN.at(idInterval).at(i), indexN.at(idInterval).at(k));
              N(indexN.at(idInterval).at(i), indexN.at(idInterval).at(k)) = Matrix();
            }
          }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsShortTimeStaticLongTime::setNull()
{
  try
  {
    MatrixDistributed::setNull();
    n.setNull();
    lPl.setNull();
    obsCount = 0;
    for(UInt idInterval=0; idInterval+1<normalsTemporal.size(); idInterval++)
      normalsTemporal.at(idInterval).setNull();
    for(UInt idInterval=0; idInterval+1<nTemporal.size(); idInterval++)
      nTemporal.at(idInterval).setNull();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsShortTimeStaticLongTime::accumulate(UInt idInterval, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    // if equations are orthogonaly transformed
    // additional residuals are appended to l
    Matrix l2;
    if(l.rows() > A.rows())
    {
      l2 = l.row(A.rows(), l.rows()-A.rows());
      l  = l.row(0, A.rows());
    }

    // eliminate arc dependent parameters
    // ----------------------------------
    Vector tau;
    if(B.size())
    {
      tau = QR_decomposition(B);
      QTransMult(B, tau, l); // transform observations: l:= Q'l
      QTransMult(B, tau, A); // transform design matrix A:=Q'A
    }
    // use only nullspace of design matrix B
    MatrixSlice A_bar( A.row(B.columns(), A.rows()-B.columns()) );
    MatrixSlice l_bar( l.row(B.columns(), l.rows()-B.columns()) );

    // build normals
    // -------------
    obsCount += l_bar.rows();
    for(UInt i=0; i<l_bar.columns(); i++)
      lPl(i) += quadsum(l_bar.column(i)) + quadsum(l2.column(i));
    for(UInt i=0; i<indexA.at(idInterval).size(); i++)
    {
      const UInt idxN1 = indexN.at(idInterval).at(i);
      const UInt idxA1 = indexA.at(idInterval).at(i);
      const Bool isTemporal1 = blockCountTemporal && (blockIndexTemporal.at(0) <= idxN1) && (idxN1 < blockIndexTemporal.at(0)+blockCountTemporal);

      // right hand sides
      Matrix &n2 = (isTemporal1) ? nTemporal.at(idInterval) : n;
      matMult(1., A_bar.column(idxA1, blockSize(idxN1)).trans(), l_bar, n2.row(blockIndex(idxN1), blockSize(idxN1)));

      // normal matrix diagonal block
      Matrix &N2 = (isTemporal1) ? normalsTemporal.at(idInterval).N(idxN1, idxN1) : N(idxN1, idxN1);
      if(N2.size() == 0)
        N2 = Matrix(blockSize(idxN1), Matrix::SYMMETRIC);
      rankKUpdate(1.0, A_bar.column(idxA1, blockSize(idxN1)), N2);

      // normal matrix, other blocks
      for(UInt k=i+1; k<indexA.at(idInterval).size(); k++)
      {
        const UInt idxN2 = indexN.at(idInterval).at(k);
        const UInt idxA2 = indexA.at(idInterval).at(k);
        const Bool isTemporal2 = blockCountTemporal && (blockIndexTemporal.at(0) <= idxN2) && (idxN2 < blockIndexTemporal.at(0)+blockCountTemporal);
        Matrix &N2 = (isTemporal1 || isTemporal2) ? normalsTemporal.at(idInterval).N(idxN1, idxN2) : N(idxN1, idxN2);
        if(N2.size() == 0)
          N2 = Matrix(blockSize(idxN1), blockSize(idxN2));
        matMult(1.0, A_bar.column(idxA1, blockSize(idxN1)).trans(), A_bar.column(idxA2, blockSize(idxN2)), N2);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsShortTimeStaticLongTime::reduceSum(Bool timing)
{
  try
  {
    Parallel::reduceSum(n,        0, communicator());
    Parallel::reduceSum(obsCount, 0, communicator());
    Parallel::reduceSum(lPl,      0, communicator());
    MatrixDistributed::reduceSum(timing);

    if(blockCountTemporal == 0)
      return;
    if(timing) logStatus<<"setup long time normal equations"<<Log::endl;

    // right hand side
    for(UInt idInterval=0; idInterval<nTemporal.size(); idInterval++)
      Parallel::reduceSum(nTemporal.at(idInterval), 0, communicator());

    if(Parallel::isMaster(communicator()))
      for(UInt blockRow=0; blockRow<blockCountTemporal; blockRow++)
        for(UInt k=0; k<blockIndexTemporal.size(); k++)
           for(UInt idInterval=0; idInterval<nTemporal.size(); idInterval++)
             axpy(factorTemporal.at(k).at(idInterval),
                  nTemporal.at(idInterval).row(blockIndex(blockRow+blockIndexTemporal.at(0)), blockSize(blockRow+blockIndexTemporal.at(0))),
                  n.row(blockIndex(blockRow-blockIndexTemporal.at(0)+blockIndexTemporal.at(k)), blockSize(blockRow+blockIndexTemporal.at(0))));

    for(UInt idInterval=0; idInterval<normalsTemporal.size(); idInterval++)
      normalsTemporal.at(idInterval).reduceSum(FALSE);

    // --- lambda -------------------
    auto send = [&](Matrix &N, UInt rankSource, UInt rankDest)
    {
      if(rankSource == rankDest)
        return;
      if(Parallel::myRank(communicator()) == rankSource)
        Parallel::send(N, rankDest, communicator());
      if(Parallel::myRank(communicator()) == rankDest)
        Parallel::receive(N, rankSource, communicator());
    };
    // ------------------------------

    for(UInt blockRow=0; blockRow<blockIndexTemporal.at(1); blockRow++)
      for(UInt blockCol=blockRow; blockCol<blockIndexTemporal.at(1); blockCol++)
      {
        // --- lambda -------------------
        auto reduceTemporal = [&](const std::vector<Double> &factor, Bool trans, UInt blocki, UInt blockk)
        {
          setBlock(blocki, blockk);
          Matrix N2;
          if(isMyRank(blockRow, blockCol))
          {
            N2 = Matrix(blockSize(blockRow), blockSize(blockCol));
            for(UInt idInterval=0; idInterval<normalsTemporal.size(); idInterval++)
              if(normalsTemporal.at(idInterval).isMyRank(blockRow, blockCol))
                axpy(factor.at(idInterval), normalsTemporal.at(idInterval).N(blockRow, blockCol), N2);
          }
          send(N2, rank(blockRow, blockCol), rank(blocki, blockk));
          if(isMyRank(blocki, blockk))
          {
            if(blockRow == blockCol)
              N2.setType(Matrix::SYMMETRIC);
            if((blockRow == blockCol) && (blocki != blockk))
            {
              fillSymmetric(N2);
              N2.setType(Matrix::GENERAL);
            }
            N(blocki, blockk) = (trans) ? N2.trans() : N2;
          }
        };
        // ------------------------------

        // upper blocks
        // ------------
        if((blockRow < blockIndexTemporal.at(0)) &&
           (blockIndexTemporal.at(0) <= blockCol) && (blockCol < blockIndexTemporal.at(0)+blockCountTemporal))
        {
          for(UInt k=0; k<blockIndexTemporal.size(); k++)
            reduceTemporal(factorTemporal.at(k), FALSE, blockRow, blockCol-blockIndexTemporal.at(0)+blockIndexTemporal.at(k));
        }

        // temporal blocks
        // ---------------
        if((blockIndexTemporal.at(0) <= blockRow) && (blockRow < blockIndexTemporal.at(0)+blockCountTemporal) &&
           (blockIndexTemporal.at(0) <= blockCol) && (blockCol < blockIndexTemporal.at(0)+blockCountTemporal))
        {
          for(UInt k=0; k<blockIndexTemporal.size(); k++)
            for(UInt l=k; l<blockIndexTemporal.size(); l++)
            {
              std::vector<Double> factor(factorTemporal.at(0).size());
              std::transform(factorTemporal.at(k).begin(), factorTemporal.at(k).end(), factorTemporal.at(l).begin(), factor.begin(), std::multiplies<Double>());
              reduceTemporal(factor, FALSE,
                             blockRow-blockIndexTemporal.at(0)+blockIndexTemporal.at(k),
                             blockCol-blockIndexTemporal.at(0)+blockIndexTemporal.at(l));
            }
        }

        // right side blocks
        // -----------------
        if((blockIndexTemporal.at(0) <= blockRow) && (blockRow < blockIndexTemporal.at(0)+blockCountTemporal) &&
           (blockIndexTemporal.at(0)+blockCountTemporal <= blockCol))
        {
          reduceTemporal(factorTemporal.at(0), FALSE, blockRow, blockCol);
          for(UInt k=1; k<blockIndexTemporal.size(); k++)
            reduceTemporal(factorTemporal.at(k), TRUE, blockCol, blockRow-blockIndexTemporal.at(0)+blockIndexTemporal.at(k));
        }

        // free memory
        for(UInt idInterval=0; idInterval<normalsTemporal.size(); idInterval++)
          if(normalsTemporal.at(idInterval).isBlockUsed(blockRow, blockCol))
            normalsTemporal.at(idInterval).N(blockRow, blockCol) = Matrix();
      }

    // fill symmetric
    // --------------
    for(UInt k=0; k<blockIndexTemporal.size(); k++)
      for(UInt l=k+1; l<blockIndexTemporal.size(); l++)
         for(UInt blockRow=0; blockRow<blockCountTemporal; blockRow++)
           for(UInt blockCol=blockRow+1; blockCol<blockCountTemporal; blockCol++)
           {
             const UInt rowSource = blockIndexTemporal.at(k)+blockRow;
             const UInt colSource = blockIndexTemporal.at(l)+blockCol;
             const UInt rowDest   = blockIndexTemporal.at(k)+blockCol;
             const UInt colDest   = blockIndexTemporal.at(l)+blockRow;

             Matrix N2 = N(rowSource, colSource).trans();
             setBlock(rowDest, colDest);
             send(N2, rank(rowSource, colSource), rank(rowDest, colDest));
             if(isMyRank(rowDest, colDest))
               N(rowDest, colDest) = N2;
           }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void NormalsShortTimeStaticLongTime::addShortTimeNormals(Double sigma2, const std::vector<std::vector<std::vector<Matrix>>> &normalsShortTime)
{
  try
  {
    if(!blockIndexShortTime.size())
      return;

    if(Parallel::isMaster(communicator()))
      obsCount += normalsShortTime.at(0).at(0).at(0).rows() * blockIndexShortTime.size();
    for(UInt id=0; id<blockIndexShortTime.size(); id++)
    {
      const UInt idx = std::min(id, normalsShortTime.size()-1);
      for(UInt i=0; i<normalsShortTime.at(idx).size(); i++)
        for(UInt k=i; k<normalsShortTime.at(idx).at(i).size(); k++)
        {
          setBlock(blockIndexShortTime.at(id+i-idx), blockIndexShortTime.at(id+k-idx));
          if(isMyRank(blockIndexShortTime.at(id+i-idx), blockIndexShortTime.at(id+k-idx)))
            axpy(1./sigma2, normalsShortTime.at(idx).at(i).at(k), N(blockIndexShortTime.at(id+i-idx), blockIndexShortTime.at(id+k-idx)));
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsShortTimeStaticLongTime::regularizeUnusedParameters(UInt countBlock)
{
  try
  {
    for(UInt i=0; i<countBlock; i++)
    {
      setBlock(i, i);
      if(isMyRank(i,i))
      {
        Matrix &N2 = N(i,i);
        for(UInt k=0; k<N2.rows(); k++)
          if(N2(k,k) == 0)
          {
            N2(k,k) = 1.;
            logWarning<<" parameters is not used: "<<parameterNames.at(blockIndex(i)+k).str()<<Log::endl;
          }
      }
    }
    Parallel::barrier(communicator());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double NormalsShortTimeStaticLongTime::solve(Matrix &x, Matrix &Wz)
{
  try
  {
    regularizeUnusedParameters(blockCount());
    x = MatrixDistributed::solve(n, TRUE/*timing*/);
    Parallel::broadCast(x, 0, communicator());

    // N contains now the cholesky decomposition
    Wz = Vce::monteCarlo(parameterCount(), 100); // monte carlo vector for VCE
    triangularSolve(Wz);
    Parallel::broadCast(Wz, 0, communicator());

    return std::sqrt((lPl(0)-inner(x.column(0), n.column(0)))/(obsCount-parameterCount()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NormalsShortTimeStaticLongTime::parameterStandardDeviation()
{
  try
  {
    cholesky2SparseInverse();

    Vector diagonal = Vector(dimension());
    for(UInt i=0; i<blockCount(); i++)
      if(isMyRank(i,i))
      {
        Matrix &N2 = N(i,i);
        for(UInt z=0; z<N2.rows(); z++)
          diagonal(blockIndex(i)+z) = std::sqrt(N2(z,z));
      }

    Parallel::reduceSum(diagonal, 0, communicator());
    return diagonal;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double NormalsShortTimeStaticLongTime::estimateShortTimeNormalsVariance(Double sigma2, const std::vector<std::vector<std::vector<Matrix>>> &normalsShortTime,
                                                                        const_MatrixSliceRef x, const_MatrixSliceRef Wz) const
{
  try
  {
    const UInt count = blockSize(blockIndexShortTime.at(0));
    Matrix Nx(x.rows(), x.columns());
    Matrix NWz(Wz.rows(), Wz.columns());
    for(UInt id=0; id<blockIndexShortTime.size(); id++)
    {
      const UInt idx = std::min(id, normalsShortTime.size()-1);
      for(UInt i=0; i<normalsShortTime.at(idx).size(); i++)
        for(UInt k=i; k<normalsShortTime.at(idx).at(i).size(); k++)
          if(isMyRank(blockIndexShortTime.at(id+i-idx), blockIndexShortTime.at(id+k-idx)))
          {
            const UInt idxi = blockIndex(blockIndexShortTime.at(id+i-idx));
            const UInt idxk = blockIndex(blockIndexShortTime.at(id+k-idx));
            matMult(1., normalsShortTime.at(idx).at(i).at(k),  x.row(idxk, count),  Nx.row(idxi, count));
            matMult(1., normalsShortTime.at(idx).at(i).at(k), Wz.row(idxk, count), NWz.row(idxi, count));
            if(k > i) // extend symmetric
            {
              matMult(1., normalsShortTime.at(idx).at(i).at(k).trans(),  x.row(idxi, count),  Nx.row(idxk, count));
              matMult(1., normalsShortTime.at(idx).at(i).at(k).trans(), Wz.row(idxi, count), NWz.row(idxk, count));
            }
          }
    } // for(id)
    Parallel::reduceSum(Nx,  0, communicator());
    Parallel::reduceSum(NWz, 0, communicator());

    const Double ePe = inner(x.column(0), Nx.column(0));
    const Double r   = count*blockIndexShortTime.size() - 1./sigma2*inner(Wz, NWz);
    Double s2  = ePe/r;
    Parallel::broadCast(s2, 0, communicator());
    return s2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsShortTimeStaticLongTime::designMatMult(UInt idInterval, Double factor, const_MatrixSliceRef A, const_MatrixSliceRef x, MatrixSliceRef Ax)
{
  try
  {
    for(UInt i=0; i<indexA.at(idInterval).size(); i++)
    {
      const UInt idxN = indexN.at(idInterval).at(i);
      const UInt idxA = indexA.at(idInterval).at(i);
      matMult(factor, A.column(idxA, blockSize(idxN)), x.row(blockIndex(idxN), blockSize(idxN)), Ax);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
