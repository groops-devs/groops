/***********************************************/
/**
* @file matrixDistributed.h
*
* @brief Positve definte matrix in distributed memory.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2011-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXDISTRIBUTED__
#define __GROOPS_MATRIXDISTRIBUTED__

#include "base/importStd.h"
#include "base/matrix.h"
#include "parallel/parallel.h"

/***** CLASS ***********************************/

/** @brief Representation of a positve definte matrix in distributed memory
* All algorithms operate only on the upper triangle of the matrix.
* @ingroup parallelGroup */
class MatrixDistributed
{
  Parallel::CommunicatorPtr             comm;
  std::function<UInt(UInt, UInt, UInt)> calcRank;
  std::vector<UInt>                    _blockIndex;

  std::vector<std::vector<std::pair<UInt, UInt>>> _row;      // each column, used row -> idx to _N and _rank
  std::vector<std::vector<std::pair<UInt, UInt>>> _column;   // each row, used column -> idx to _N and _rank
  std::vector<Matrix>               _N;        // unorderd list of used blocks
  std::vector<UInt>                 _rank;     // unorderd list of rank of used blocks

  Bool isMyRank(UInt idx) const {return (Parallel::myRank(comm) == _rank[idx]);}
  UInt index(UInt row, UInt col) const;        // NULLINDEX if not set
  // call: loopBlockColumn({rowStart, rowPastEnd}, col, [&](UInt row, UInt idx) {_N[idx] == N(row,col)});
  void loopBlockColumn(const std::array<UInt,2> &rows, UInt col, std::function<void(UInt, UInt)> block) const;
  // call: loopBlockRow(row, {colStart, colPastEnd}, [&](UInt col, UInt idx) {_N[idx] == N(row,col)});
  void loopBlockRow   (UInt row, const std::array<UInt,2> &cols, std::function<void(UInt, UInt)> block) const;

  // parallel
  // call: usedRanks = usedRanksInColumn({rowStart, rowPastEnd}, col);
  std::vector<Bool> usedRanksInColumn(const std::array<UInt,2> &rows, UInt col) const;
  // call: usedRanks = usedRanksInRow(row, {colStart, colPastEnd});
  std::vector<Bool> usedRanksInRow(UInt row, const std::array<UInt,2> &cols) const;
  void broadCast(Matrix &x, UInt idx, const std::vector<Bool> &usedRank);
  void reduceSum(Matrix &x, UInt idx, const std::vector<Bool> &usedRank, Bool free=TRUE);


  /* @brief Solve an triangular system of equations \f$ \mathbf{W}\mathbf{y} = \mathbf{x}\f$
  * The input must be initialized at all nodes.
  * The sum other all nodes defines the input (reduceSum is called internally)
  * Output is valid at master only. */
  void triangularSolve(std::vector<Matrix> &x) {triangularSolve(x, 0, blockCount());}
  void triangularSolve(std::vector<Matrix> &x, UInt startBlock, UInt countBlock);

  /* @brief Solve a part of a triangular system of equations \f$ \mathbf{W}^T\mathbf{y} = \mathbf{x}\f$
  * This corresponds to the partly @a cholesky function.
  * The input must be initialized at all nodes.
  * The sum other all nodes defines the input (reduceSum is called internally)
  * Output is valid at master only. */
  void triangularTransSolve(std::vector<Matrix> &x) {triangularTransSolve(x, 0, blockCount(), TRUE);}
  void triangularTransSolve(std::vector<Matrix> &x, UInt startBlock, UInt countBlock, Bool collect);

public:
  /// Default constructor.
  MatrixDistributed();

  /** @brief Initialize matrix structure (number of blocks, block sizes, communicator, and rank calculation).
  * The @a blockIndex is a vector of ascending indices starting with zero, where each adjacent pair of elements define the position and shape of the subblocks.
  * Hence, all elements \f$ n_{ij} \f$  which satisfy \f$ blockIndex[k] \leq  i < blockIndex[k+1] \land blockIndex[l] \leq j < blockIndex[l+1] \f$ are found in block (@a k, @a l).
  * @param blockIndex: boundary indices of the sub-blocks.
  * @param comm: Parallel communicator of the matrix (default: MPI_COMM_WORLD).
  * @param calcRank: function handler to determine the process rank of block(i,k) (default: block cyclic distribution). */
  explicit MatrixDistributed(const std::vector<UInt> &blockIndex, Parallel::CommunicatorPtr comm, const std::function<UInt(UInt, UInt, UInt)> &calcRank=nullptr);

  /** @copydoc MatrixDistributed(const std::vector<UInt> &, Parallel::CommunicatorPtr, std::function<UInt(UInt, UInt, UInt)>)
  * This method allocates the upper block triangle of a symmetric matrix with zero matrices. */
  void init(const std::vector<UInt> &blockIndex, Parallel::CommunicatorPtr comm, const std::function<UInt(UInt, UInt, UInt)> &calcRank=nullptr);

  /** @copydoc MatrixDistributed(const std::vector<UInt> &, Parallel::CommunicatorPtr, std::function<UInt(UInt, UInt, UInt)>)
  * No memory is allocated. To assign blocks to processes and allocate memory, @a setBlock(UInt i, UInt k, UInt rank) has to be called. */
  void initEmpty(const std::vector<UInt> &blockIndex, Parallel::CommunicatorPtr comm, const std::function<UInt(UInt, UInt, UInt)> &calcRank=nullptr);

  // =========================================

  /** @brief Set a new handler to calculate the rank of new blocks.
  * If @p calcRank is nullptr a default block cyclic distribution is assumed.
  * This function must be called by all processes within the communcator in the matrix! */
  void setCalculateRank(const std::function<UInt(UInt, UInt, UInt)> &calcRank);

  /** @brief returns the handler to calculate the rank of new blocks. */
  const std::function<UInt(UInt, UInt, UInt)> &getCalculateRank() const {return calcRank;}

  /** @brief The default calculateRank handler.
  * Block cyclic distribution.
  * This function must be called by all processes within the communcator in the matrix! */
  static UInt calculateRankBlockCyclic(UInt i, UInt k, UInt commSize);

  /** @brief Assigns block (@a i, @a k) to process @a rank and allocates memory (only at the parent process).
  * If no rank is specified, the calculateRank handler is called for new blocks.
  * This function must be called by all processes within the communcator in the matrix! */
  UInt setBlock(UInt i, UInt k, UInt rank=NULLINDEX);

  /** @brief Removes the block rows/columns.
  * This function erases the block rows/columns containing the diagonal blocks @a startIndex to @a startIndex+count-1
  * effectively reducing the size by blockSize[startIndex] + ... + blockSize[startIndex+count-1] elements.
  * This function must be called by all processes within the communcator in the matrix!
  * @param startIndex index of block row/column to remove
  * @param count number of rows/columns to remove  */
  void eraseBlocks(UInt startIndex, UInt count=1);

  // =========================================

  Parallel::CommunicatorPtr communicator() const {return comm;}
  std::vector<UInt>         blockIndex()   const {return _blockIndex;}

  UInt parameterCount()            const {return _blockIndex.at(blockCount());}          //!< Number of rows/columns (dimension) of distributed matrix
  UInt dimension()                 const {return parameterCount();}                      //!< Number of rows/columns (dimension) of distributed matrix
  UInt blockIndex(UInt i)          const {return _blockIndex.at(i);}                     //!< Start index of block @a i
  UInt blockSize(UInt i)           const {return _blockIndex.at(i+1)-_blockIndex.at(i);} //!< Size of block @a i
  UInt blockCount()                const {return _blockIndex.size()-1;}                  //!< Number of block rows/columns
  UInt rank(UInt i, UInt k)        const;                                                //!< Rank of process holding block (@a i, @a k)
  Bool isMyRank(UInt i, UInt k)    const {return (Parallel::myRank(comm) == rank(i,k));} //!< Returns TRUE if the calling process holds block (@a i, @a k) within the underlying communicator of the matrix
  Bool isBlockUsed(UInt i, UInt k) const {return (rank(i,k) != NULLINDEX);}              //!< Returns TRUE if block (@a i, @a k) is assigned to a process
  UInt index2block(UInt i)         const;                                                //!< Returns the index of the block which holds element @a i

  Matrix       &N(UInt i, UInt k);       //!< Returns a writable reference to block (@a i, @a k). Throws an exception if the block is not assigned to a process.
  const Matrix &N(UInt i, UInt k) const; //!< Returns a read only reference to block (@a i, @a k). Throws an exception if the block is not assigned to a process.

  /// Fill all matrix blocks with zero.
  void setNull();

  /// Reduce block (@a i, @a k) on its parent process. After the operation, the memory on all other processes is freed.
  void reduceSum(UInt i, UInt k);

  /// Reduce all assigned blocks of the matrix on their parent processes. After the operation, the memory on all other processes is freed.
  void reduceSum(Bool timing=TRUE);

  // =========================================

  /** @brief (In-place) Cholesky decomposition of the distributed matrix. */
  void cholesky(Bool timing=TRUE) {cholesky(timing, 0, blockCount(), TRUE);}

  /** @brief Performs a part of the Cholesky decomposition.
  * The Cholesky decomposition must be already performed for the blocks before @p startBlock.
  * The blocks after @p startBlock + @p countBlock contain at ouput the normal matrix where all parameters before
  * are eliminated. If not @p collect, @a reduceSum must be called afterwards for these blocks. */
  void cholesky(Bool timing, UInt startBlock, UInt countBlock, Bool collect);

  /** @brief Solve the system of equations \f$ \mathbf{N}\mathbf{x} = \mathbf{n}\f$
  * Performs @a cholesky, @a triangularTransSolve, and @a triangularSolve.
  * The input must be valid at master only. Output is valid at master only. */
  Matrix solve(const_MatrixSliceRef n, Bool timing=TRUE);

  /** @brief Solve a triangular system of equations \f$ \mathbf{W}\mathbf{y} = \mathbf{x}\f$
  * \f$ \mathbf{W} \f$ is assumed to be an upper triangular matrix.
  * The input must be valid at master only. Output is valid at master only. */
  void triangularSolve(MatrixSliceRef x);

  /** @brief Solve a triangular system of equations \f$ \mathbf{W}^T\mathbf{y} = \mathbf{x}\f$
  * \f$ \mathbf{W} \f$ is assumed to be an upper triangular matrix.
  * The input must be valid at master only. Output is valid at master only. */
  void triangularTransSolve(MatrixSliceRef x) {triangularTransSolve(x, 0, blockCount());}
  void triangularTransSolve(MatrixSliceRef x, UInt startBlock, UInt countBlock);

  /** @brief Compute the inverse of the (upper triangular) Cholesky  factor \f$ \mathbf{W}\f$ */
  void choleskyInverse(Bool timing=TRUE) {choleskyInverse(timing, 0, blockCount());}
  void choleskyInverse(Bool timing, UInt startBlock, UInt countBlock);

  /** @brief Compute product \f$ \mathbf{U}\mathbf{U}^T \f$, where \f$ \mathbf{U} \f$ is either a Cholesky factor or its inverse. */
  void choleskyProduct(Bool timing=TRUE);

  /** @brief Compute the sparse inverse of the matrix from the Cholesky  factor \f$ \mathbf{W}\f$.
  * Computes the upper triangle of the inverse of the symmetric matrix. The structure of the Cholesky factor is preserved.  */
  void cholesky2SparseInverse(Bool timing=TRUE);

  // =========================================

  /** @brief Reorder parameters in a symmetric matrix.
  * The parameter @a index contains the indices of the elements in the reordered matrix.
  * Indices with NULLINDEX are inserted as zero elements into the new matrix.
  * @param index: indices of the original matrix elements in the reordered matrix
  * @param blockIndex: boundary indices of the sub-blocks.
  * @param calcRank: function handler to determine the process rank of block(i,k) (default: block cyclic distribution). */
  void reorder(const std::vector<UInt> &index, const std::vector<UInt> &blockIndex, const std::function<UInt(UInt, UInt, UInt)> &calcRank=nullptr);

  // =========================================

  /** @brief Compute boundary indices for distributed blocks from parameter count and block size.
  * If blockSize is zero, the matrix consists of a single block. */
  static std::vector<UInt> computeBlockIndex(UInt parameterCount, UInt blockSize=2048);

  friend class GnssProcessingStep;
  friend class GnssParametrizationAmbiguities;
};

//***********************************************/

#endif
