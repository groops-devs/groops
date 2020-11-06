/***********************************************/
/**
* @file parallel.h
*
* @brief Wrapper for Message Passing Interface (MPI).
* All functions are empty statements in case
* of the single processor version (parallelSingle.cpp).
*
* @author Torsten Mayer-Guerr
* @date 2004-11-13
*
*/
/***********************************************/

#ifndef __GROOPS_PARALLEL__
#define __GROOPS_PARALLEL__

#include "base/import.h"
#include "inputOutput/archiveBinary.h"
#include "inputOutput/logging.h"

/***********************************************/

class GnssType;

/***********************************************/

/** @brief Wrapper for Message Passing Interface (MPI).
* All functions are empty statements in case
* of the single processor version (parallelSingle.cpp).
* @ingroup parallelGroup */
namespace Parallel
{
  class Communicator;
  typedef std::shared_ptr<Communicator> CommunicatorPtr;

  /** @brief Must be called firstly in main. */
  void init(int argc, char *argv[]);

  /** @brief Must be called at the end of main. */
  void finalize();

  /** @brief Best attempt to abort tasks in comm. */
  void abort();

  // =========================================================

  /** @brief The global communicator.
  * Equivalent to MPI_COMM_WORLD. */
  CommunicatorPtr globalCommunicator();

  /** @brief The default communicator.
  * It is used if comm is a nullptr. */
  CommunicatorPtr defaultCommunicator();

  /** @brief Set a new default communicator.
  * @return the old communicator */
  CommunicatorPtr setDefaultCommunicator(CommunicatorPtr comm);

  /** @brief Creates new communicators.
  * a new group is created for each different @a color.
  * the ranks in the groups are sorted by the @a key. */
  CommunicatorPtr splitCommunicator(UInt color, UInt key, CommunicatorPtr comm=nullptr);

  /** @brief Creates new communicators.
  * Must be called by every process in @a comm. */
  CommunicatorPtr createCommunicator(std::vector<UInt> ranks, CommunicatorPtr comm=nullptr);

  /** @brief The communicator that refers to the own process only. */
  CommunicatorPtr selfCommunicator();

  // =========================================================

  /** @brief Number of processes. */
  UInt size(CommunicatorPtr comm=nullptr);

  /** @brief Process index. */
  UInt myRank(CommunicatorPtr comm=nullptr);

  /** @brief Is ths the master process (rank==0)? */
  inline Bool isMaster(CommunicatorPtr comm=nullptr) {return (myRank(comm) == 0);}

  /** @brief Blocks until all process have reached this routine. */
  void barrier(CommunicatorPtr comm=nullptr);

  // =========================================================

  /** @brief Send raw data @a x to process with rank @a process. */
  void send(const Byte *x, UInt size, UInt process, CommunicatorPtr comm=nullptr);

  /** @brief receive raw data @a x from prozess with rank @a process.
  * If @a process = NULLINDEX then receive from an arbitrary process. */
  void receive(Byte *x, UInt size, UInt process, CommunicatorPtr comm=nullptr);

  /** @brief Distribute raw data @a x at @a process to all other processes. */
  void broadCast(Byte *x, UInt size, UInt process, CommunicatorPtr comm=nullptr);

  // =========================================================

  /** @brief Send @a x to process with rank @a process. */
  ///@{
  template<typename T> void send(const T &x, UInt process, CommunicatorPtr comm=nullptr);
  template<> void send(const UInt     &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Double   &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Bool     &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Angle    &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Time     &x, UInt process, CommunicatorPtr comm);
  template<> void send(const GnssType &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Vector3d &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Vector   &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Matrix   &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief receive @a x from prozess with rank @a process.
  * If @a process = NULLINDEX then receive from an arbitrary process. */
  ///@{
  template<typename T> void receive(T &x, UInt process, CommunicatorPtr comm=nullptr);
  template<> void receive(UInt     &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Double   &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Bool     &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Angle    &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Time     &x, UInt process, CommunicatorPtr comm);
  template<> void receive(GnssType &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Vector3d &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Vector   &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Matrix   &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief Distribute @a x at @a process to all other processes. */
  ///@{
  template<typename T> void broadCast(T &x, UInt process=0, CommunicatorPtr comm=nullptr);
  template<> void broadCast(UInt     &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Double   &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Bool     &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Angle    &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Time     &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(GnssType &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Vector3d &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Vector   &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Matrix   &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief Sum up @a x at all processes (also rank 0) and send the result to @a process. */
  ///@{
  void reduceSum(UInt    &x, UInt process=0, CommunicatorPtr comm=nullptr);
  void reduceSum(Double  &x, UInt process=0, CommunicatorPtr comm=nullptr);
  void reduceSum(Bool    &x, UInt process=0, CommunicatorPtr comm=nullptr);
  void reduceSum(Matrix  &x, UInt process=0, CommunicatorPtr comm=nullptr);
  ///@}

  /** @brief Find min/max of @a x at all processes (also rank 0) and send the result to @a process. */
  ///@{
  void reduceMin(UInt   &x, UInt process=0, CommunicatorPtr comm=nullptr);
  void reduceMin(Double &x, UInt process=0, CommunicatorPtr comm=nullptr);
  void reduceMax(UInt   &x, UInt process=0, CommunicatorPtr comm=nullptr);
  void reduceMax(Double &x, UInt process=0, CommunicatorPtr comm=nullptr);
  ///@}

  // =========================================================

  /** @brief Parallelized loop.
  * Calls @a func(i) for every @a i in [0,count).
  * The different calls are distributed other the processes (without master).
  * @return The process number for @a i is returned (valid at master). */
  template<typename T> std::vector<UInt> forEach(UInt count, T func, CommunicatorPtr comm=nullptr, Bool timing=TRUE);

  /** @brief Parallelized loop.
  * Calls @a vec[i]=func(i) for every @a i in [0,vec.size()).
  * The different calls are distributed other the processes (without master).
  * The result in @a vec is only valid at master.
  * @return The process number for @a i is returned (valid at master). */
  template<typename A, typename T> std::vector<UInt> forEach(std::vector<A> &vec, T func, CommunicatorPtr comm=nullptr, Bool timing=TRUE);

  /** @brief Parallelized loop.
  * Calls @a func(i) for every @a i in [0,count).
  * The different calls are distributed other the processes (without master).
  * @return The process number for @a i is returned (valid at master). */
  template<typename T> std::vector<UInt> forEachInterval(UInt count, const std::vector<UInt> &interval, T func, CommunicatorPtr comm=nullptr, Bool timing=TRUE);

  /** @brief Parallelized loop.
  * Calls @a vec[i]=func(i) for every @a i in [0,vec.size()).
  * The different calls are distributed other the processes (without master).
  * The result in @a vec is only valid at master.
  * @return The process number for @a i is returned (valid at master). */
  template<typename A, typename T> std::vector<UInt> forEachInterval(std::vector<A> &vec, const std::vector<UInt> &interval, T func, CommunicatorPtr comm=nullptr, Bool timing=TRUE);

  /** @brief Parallelized loop.
  * Calls @a func(i) for every @a i in [0,count).
  * The different calls are distributed using @a processNo (without master).
  * The result in @a vec is only valid at master. */
  template<typename T> void forEachProcess(UInt count, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm=nullptr, Bool timing=TRUE);

  /** @brief Parallelized loop.
  * Calls @a vec[i]=func(i) for every @a i in [0,vec.size()).
  * The different calls are distributed using @a processNo (without master).
  * The result in @a vec is only valid at master. */
  template<typename A, typename T> void forEachProcess(std::vector<A> &vec, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm=nullptr, Bool timing=TRUE);
} // end namespace Parallel

/***********************************************/
/***** INLINES ***********************************/
/***********************************************/

template<typename T>
inline void Parallel::send(const T &x, UInt process, CommunicatorPtr comm)
{
  if(size(comm)<=1)
    return;
  std::stringstream stream; //(std::ios::binary);
  OutArchiveBinary oa(stream, "", MAX_UINT);
  oa<<nameValue("xxx", x);
  std::string str = stream.str();
  UInt size = str.size();
  send(size, process, comm);
  send(str.data(), size, process, comm);
}

/***********************************************/

template<typename T>
inline void Parallel::receive(T &x, UInt process, CommunicatorPtr comm)
{
  if(size(comm)<=1)
    return;
  UInt size;
  receive(size, process, comm);
  Byte *str = new Byte[size+1];
  receive(str, size, process, comm);
  std::stringstream stream(std::string(str, size)); //, std::ios::binary);
  InArchiveBinary ia(stream);
  ia>>nameValue("xxx", x);
  delete[] str;
}

/***********************************************/

template<typename T>
inline void Parallel::broadCast(T &x, UInt process, CommunicatorPtr comm)
{
  if(size(comm)<=1)
    return;
  if(Parallel::myRank(comm) == process)
  {
    std::stringstream stream; //(std::ios::binary);
    OutArchiveBinary oa(stream, "", MAX_UINT);
    oa<<nameValue("xxx", x);
    std::string str = stream.str();
    UInt size = str.size();
    broadCast(size, process, comm);
    broadCast(const_cast<Byte*>(str.data()), size, process, comm);
  }
  else
  {
    UInt size;
    broadCast(size, process, comm);
    Byte *str = new Byte[size+1];
    broadCast(str, size, process, comm);
    std::stringstream stream(std::string(str, size)); //, std::ios::binary);
    InArchiveBinary ia(stream);
    ia>>nameValue("xxx", x);
    delete[] str;
  }
}

/***********************************************/
/***********************************************/

template<typename T>
inline std::vector<UInt> Parallel::forEach(UInt count, T func, CommunicatorPtr comm, Bool timing)
{
  return forEachInterval(count, {0, count}, func, comm, timing);
}

/***********************************************/

template<typename A, typename T>
inline std::vector<UInt> Parallel::forEach(std::vector<A> &vec, T func, CommunicatorPtr comm, Bool timing)
{
  return forEachInterval(vec, {0, vec.size()}, func, comm, timing);
}

/***********************************************/

template<typename T>
inline std::vector<UInt> Parallel::forEachInterval(UInt count, const std::vector<UInt> &interval, T func, CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::vector<UInt> processNo(count, 0);

    // single process version
    // ----------------------
    if(size(comm)<3)
    {
      // single process version
      if(isMaster(comm))
      {
        if(timing) logTimerStart;
        for(UInt i=0; i<count; i++)
        {
          if(timing) logTimerLoop(i,count);
          func(i);
        }
        if(timing) logTimerLoopEnd(count);
      }
      return processNo;
    }

    if(count!=interval.back())
      throw(Exception("interval size and count differ"));

    // parallel version
    // ----------------
    if(isMaster(comm))
    {
      std::vector<UInt> countInInterval(interval.size()-1, 0);
      std::vector<UInt> processedInterval(size(comm), NULLINDEX);

      // master distributes the loop numbers
      UInt process, index;
      if(timing) logTimerStart;
      for(UInt i=0; i<count; i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process

        // can we compute func in the same interval?
        UInt idInterval = processedInterval.at(process);
        if((idInterval==NULLINDEX) || (countInInterval.at(idInterval) >= interval.at(idInterval+1)-interval.at(idInterval)))
        {
          // search new interval to compute
          UInt maxLeft = 0;
          for(UInt k=0; k<countInInterval.size(); k++)
          {
            UInt left = interval.at(k+1)-interval.at(k)-countInInterval.at(k);
            // interval not used?
            if((countInInterval.at(k) == 0) && (left>0))
            {
              idInterval = k;
              break;
            }
            if(left>maxLeft)
            {
              maxLeft = left;
              idInterval = k;
            }
          }
          processedInterval.at(process) = idInterval;
        }

        UInt id = interval.at(idInterval) + countInInterval.at(idInterval);
        countInInterval.at(idInterval)++;
        send(id, process, comm);       // send new loop number to be computed at process
        processNo.at(id) = process;
        if(timing) logTimerLoop(i,count);
      }
      // send to all processes the end signal (NULLINDEX)
      for(UInt i=1; i<size(comm); i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process
        send(NULLINDEX, process, comm);    // end signal
      }
      if(timing) logTimerLoopEnd(count);
    }
    else // clients
    {
      send(myRank(comm), 0, comm);
      send(NULLINDEX, 0, comm); // no results computed yet
      for(;;)
      {
        UInt i;
        receive(i,0, comm);
        if(i==NULLINDEX)
          break;
        func(i);
        send(myRank(comm), 0, comm);
        send(i, 0, comm);
      }
    }

    broadCast(processNo, 0, comm);
    return processNo;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename A, typename T>
inline std::vector<UInt> Parallel::forEachInterval(std::vector<A> &vec, const std::vector<UInt> &interval, T func, CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::vector<UInt> processNo(vec.size(), 0);

    // single process version
    // ----------------------
    if(size(comm)<3)
    {
      // single process version
      if(isMaster(comm))
      {
        if(timing) logTimerStart;
        for(UInt i=0; i<vec.size(); i++)
        {
          if(timing) logTimerLoop(i,vec.size());
          vec[i] = func(i);
        }
        if(timing) logTimerLoopEnd(vec.size());
      }
      return processNo;
    }

    if(vec.size()!=interval.back())
      throw(Exception("interval size and vec.size() differ"));

    // parallel version
    // ----------------
    if(isMaster(comm))
    {
      std::vector<UInt> countInInterval(interval.size()-1, 0);
      std::vector<UInt> processedInterval(size(comm), NULLINDEX);

      // master distributes the loop numbers
      UInt process, index;
      if(timing) logTimerStart;
      for(UInt i=0; i<vec.size(); i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process
        if(index!=NULLINDEX)
        {
          receive(vec[index], process, comm); // receive result
        }

        // can we compute func in the same interval?
        UInt idInterval = processedInterval.at(process);
        if((idInterval==NULLINDEX) || (countInInterval.at(idInterval) >= interval.at(idInterval+1)-interval.at(idInterval)))
        {
          // search new interval to compute
          UInt maxLeft = 0;
          for(UInt k=0; k<countInInterval.size(); k++)
          {
            UInt left = interval.at(k+1)-interval.at(k)-countInInterval.at(k);
            // interval not used?
            if((countInInterval.at(k) == 0) && (left>0))
            {
              idInterval = k;
              break;
            }
            if(left>maxLeft)
            {
              maxLeft = left;
              idInterval = k;
            }
          }
          processedInterval.at(process) = idInterval;
        }

        UInt id = interval.at(idInterval) + countInInterval.at(idInterval);
        countInInterval.at(idInterval)++;
        send(id, process, comm);           // send new loop number to be computed at process
        processNo.at(id) = process;
        if(timing) logTimerLoop(i,vec.size());
      }
      // send to all processes the end signal (NULLINDEX)
      for(UInt i=1; i<size(comm); i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process
        if(index!=NULLINDEX)
          receive(vec[index], process, comm); // receive result
        send(NULLINDEX, process, comm);    // end signal
      }
      if(timing) logTimerLoopEnd(vec.size());
    }
    else // clients
    {
      send(myRank(comm), 0, comm);
      send(NULLINDEX, 0, comm); // no results computed yet
      for(;;)
      {
        UInt i;
        receive(i,0, comm);
        if(i==NULLINDEX)
          break;
        vec[i] = func(i);
        send(myRank(comm), 0, comm);
        send(i, 0, comm);
        send(vec[i], 0, comm);
      }
    }

    broadCast(processNo, 0, comm);
    return processNo;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename T>
inline void Parallel::forEachProcess(UInt count, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm, Bool timing)
{
  try
  {
    UInt idx;
    if(timing) logTimerStart;
    for(UInt i=0; i<count; i++)
    {
      if(timing) logTimerLoop(i, count);
      if(myRank(comm) == processNo.at(i))
      {
        func(i);
        if(!isMaster(comm)) send(i, 0, comm);
      } // if(arcs.at(i))
      else if(isMaster(comm))
      {
        receive(idx, NULLINDEX, comm);
      }
    } // for(i)
    if(timing) logTimerLoopEnd(count);
    barrier(comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename A, typename T>
inline void Parallel::forEachProcess(std::vector<A> &vec, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm, Bool timing)
{
  try
  {
    UInt idx = 0;
    if(timing) logTimerStart;
    for(UInt i=0; i<vec.size(); i++)
    {
      if(timing) logTimerLoop(i, vec.size());
      if(myRank(comm) == processNo.at(i))
      {
        vec[i] = func(i);
        if(!isMaster(comm)) send(i, 0, comm);
        if(!isMaster(comm)) send(vec[i], 0, comm);
      } // if(arcs.at(i))
      else if(isMaster(comm))
      {
        receive(idx, NULLINDEX, comm);
        receive(vec[idx], processNo.at(idx), comm);
      }
    } // for(i)
    if(timing) logTimerLoopEnd(vec.size());
    barrier(comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
