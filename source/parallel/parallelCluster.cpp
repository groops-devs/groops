/***********************************************/
/**
* @file parallelCluster.cpp
*
* @brief Wrapper for Message Passing Interface (MPI).
*
* @author Torsten Mayer-Guerr
* @date 2004-11-13
*
*/
/***********************************************/

#include "base/import.h"
#include "base/gnssType.h"
#include "parallel/parallel.h"

// #undef SEEK_SET
// #undef SEEK_END
// #undef SEEK_CUR
#include <mpi.h>

/***** FUNCTIONS *******************************/

namespace Parallel
{

class Communicator
{
public:
  MPI_Comm comm;

  Communicator(const MPI_Comm &comm_=MPI_COMM_WORLD) : comm(comm_) {}
 ~Communicator()
  {
    if((comm != MPI_COMM_WORLD) && (comm != MPI_COMM_SELF))
    {
      MPI_Group group;
      MPI_Comm_group(comm, &group);
      MPI_Comm_free(&comm);
      MPI_Group_free(&group);
    }
  }
};

static CommunicatorPtr commDefault;

/***********************************************/
/***********************************************/

void init(int argc, char *argv[])
{
  try
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    commDefault = globalCommunicator();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void finalize()
{
  MPI_Finalize();
}

/***********************************************/

void abort()
{
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}


/***********************************************/
/***********************************************/

inline MPI_Comm getComm(CommunicatorPtr comm)
{
  return (comm != nullptr) ? comm->comm : commDefault->comm;
}

/***********************************************/

inline void check(int errorcode)
{
  if(errorcode != MPI_SUCCESS)
  {
    int  resultlen;
    char errorStr[MPI_MAX_ERROR_STRING];
    MPI_Error_string(errorcode, errorStr, &resultlen);
    throw(Exception("MPI Error: "+std::string(errorStr,resultlen)));
  }
}

/***********************************************/
/***********************************************/

CommunicatorPtr globalCommunicator()
{
  return std::make_shared<Communicator>(MPI_COMM_WORLD);
}

/***********************************************/

CommunicatorPtr defaultCommunicator()
{
  return commDefault;
}

/***********************************************/

CommunicatorPtr setDefaultCommunicator(CommunicatorPtr comm)
{
  try
  {
    if(comm == nullptr)
      throw(Exception("null pointer"));
    std::swap(commDefault, comm);
    return comm;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

CommunicatorPtr splitCommunicator(UInt color, UInt key, CommunicatorPtr comm)
{
  try
  {
    CommunicatorPtr commNew = std::make_shared<Communicator>(MPI_COMM_WORLD);
    check(MPI_Comm_split(getComm(comm), ((color==NULLINDEX) ? (MPI_UNDEFINED) : (color)), key, &commNew->comm));
    if(color==NULLINDEX)
      return nullptr;
    return commNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

CommunicatorPtr selfCommunicator()
{
  try
  {
    CommunicatorPtr commNew = std::make_shared<Communicator>();
    commNew->comm = MPI_COMM_SELF;
    return commNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

CommunicatorPtr createCommunicator(std::vector<UInt> ranks, CommunicatorPtr comm)
{
  try
  {
    CommunicatorPtr commNew = std::make_shared<Communicator>();
    std::vector<int> mpiRanks(ranks.size());
    for(UInt i=0; i<ranks.size(); i++)
      mpiRanks.at(i) = static_cast<int>(ranks.at(i));

    MPI_Group group, newgroup;
    check(MPI_Comm_group(getComm(comm), &group));
    check(MPI_Group_incl(group, mpiRanks.size(), mpiRanks.data(), &newgroup));
#if MPI_VERSION >= 3
    check(MPI_Comm_create_group(getComm(comm), newgroup, 99, &commNew->comm));
#else
    check(MPI_Comm_create(getComm(comm), newgroup, &commNew->comm));
#endif
    check(MPI_Group_free(&group));
    check(MPI_Group_free(&newgroup));
    if(commNew->comm == MPI_COMM_NULL)
        return nullptr;
    return commNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

UInt myRank(CommunicatorPtr comm)
{
  int rank;
  check(MPI_Comm_rank(getComm(comm), &rank));
  return static_cast<UInt>(rank);
}

/***********************************************/

UInt size(CommunicatorPtr comm)
{
  int size;
  check(MPI_Comm_size(getComm(comm), &size));
  return static_cast<UInt>(size);
}

/***********************************************/

void barrier(CommunicatorPtr comm)
{
  try
  {
    check(MPI_Barrier(getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void send(const Byte *x, UInt size, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Send(x, size, MPI_CHAR, process, 11, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void receive(Byte *x, UInt size, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Status status;
    check(MPI_Recv(x, size, MPI_CHAR, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 11, getComm(comm), &status));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void broadCast(Byte *x, UInt size, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Bcast(x, size, MPI_CHAR, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void send(const UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    check(MPI_Send(&y, 1, MPI_UNSIGNED, process, 1, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Send(&x, 1, MPI_DOUBLE, process, 2, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    check(MPI_Send(&y, 1, MPI_UNSIGNED, process, 2, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Angle &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Send(&x, 1, MPI_DOUBLE, process, 2, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Time &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(static_cast<Double>(x.mjdInt()), process, comm);
    send(x.mjdMod(), process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const GnssType &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x.type);
    check(MPI_Send(&y, 1, MPI_UNSIGNED, process, 1, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Vector3d &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Send(x.vector().field(), 3, MPI_DOUBLE, process, 3, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Vector &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(x.rows(), process, comm);
    if(x.size()!=0)
      check(MPI_Send(x.field(), x.size(), MPI_DOUBLE, process, 4, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt type  = static_cast<UInt>(x.getType());
    UInt upper = x.isUpper();
    send(x.rows(),    process, comm);
    send(x.columns(), process, comm);
    send(type,        process, comm);
    send(upper,       process, comm);
    if(x.size()!=0)
      check(MPI_Send(x.field(), x.size(), MPI_DOUBLE, process, 4, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void receive(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y;
    MPI_Status status;
    check(MPI_Recv(&y, 1, MPI_UNSIGNED, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 1, getComm(comm), &status));
    x = static_cast<UInt>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Status status;
    check(MPI_Recv(&x, 1, MPI_DOUBLE, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 2, getComm(comm), &status));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y;
    MPI_Status status;
    check(MPI_Recv(&y, 1, MPI_UNSIGNED, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 1, getComm(comm), &status));
    x = static_cast<Bool>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Angle &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Status status;
    check(MPI_Recv(&x, 1, MPI_DOUBLE, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 2, getComm(comm), &status));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Time &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    Double mjdInt, mjdMod;
    receive(mjdInt, process, comm);
    receive(mjdMod, process, comm);
    x = Time(static_cast<Int>(mjdInt), mjdMod);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Vector3d &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    Double v[3];
    MPI_Status status;
    check(MPI_Recv(v, 3, MPI_DOUBLE, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 3, getComm(comm), &status));
    x.x() = v[0];
    x.y() = v[1];
    x.z() = v[2];
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Vector &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Status status;
    UInt rows;
    receive(rows, process, comm);
    x = Vector(rows);
    if(x.size()!=0)
      check(MPI_Recv(x.field(), x.size(), MPI_DOUBLE, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 4, getComm(comm), &status));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Status status;
    UInt rows, columns;
    UInt  type, upper;
    receive(rows,    process, comm);
    receive(columns, process, comm);
    receive(type,    process, comm);
    receive(upper,   process, comm);
    x = Matrix(rows, columns);
    x.setType(static_cast<Matrix::Type>(type), (upper) ? Matrix::UPPER : Matrix::LOWER);
    if(x.size()!=0)
      check(MPI_Recv(x.field(), x.size(), MPI_DOUBLE, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 4, getComm(comm), &status));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(GnssType &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y;
    MPI_Status status;
    check(MPI_Recv(&y, 1, MPI_UNSIGNED, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 1, getComm(comm), &status));
    x = GnssType(static_cast<UInt>(y));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void broadCast(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    check(MPI_Bcast(&y, 1, MPI_UNSIGNED, process, getComm(comm)));
    x = static_cast<UInt>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Bcast(&x, 1, MPI_DOUBLE, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    check(MPI_Bcast(&y, 1, MPI_UNSIGNED, process, getComm(comm)));
    x = static_cast<Bool>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Angle &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    check(MPI_Bcast(&x, 1, MPI_DOUBLE, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Time &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    Double mjdInt = x.mjdInt();
    Double mjdMod = x.mjdMod();
    broadCast(mjdInt, process, comm);
    broadCast(mjdMod, process, comm);
    x = Time(static_cast<Int>(mjdInt), mjdMod);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(GnssType &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x.type);
    check(MPI_Bcast(&y, 1, MPI_UNSIGNED, process, getComm(comm)));
    x = GnssType(static_cast<UInt>(y));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Vector &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt rows = x.rows();
    broadCast(rows, process, comm);
    if(x.rows() != rows)
      x = Vector(rows);
    if(x.size()!=0)
      check(MPI_Bcast(x.field(), x.size(), MPI_DOUBLE, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt rows = x.rows();
    UInt cols = x.columns();
    UInt type  = static_cast<UInt>(x.getType());
    UInt upper = x.isUpper();

    broadCast(rows,  process, comm);
    broadCast(cols,  process, comm);
    broadCast(type,  process, comm);
    broadCast(upper, process, comm);

    if((x.rows() != rows) || (x.columns() != cols))
      x = Matrix(rows, cols);
    x.setType(static_cast<Matrix::Type>(type), (upper) ? Matrix::UPPER : Matrix::LOWER);

    constexpr UInt BLOCKSIZE = 200*1024*1024/sizeof(Double); // 200 Mb
    UInt index = 0;
    while(index<x.size())
    {
      const UInt size = std::min(x.size()-index, BLOCKSIZE);
      check(MPI_Bcast(x.field()+index, size, MPI_DOUBLE, process, getComm(comm)));
      index += size;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void reduceSum(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int >(x);
      unsigned int y = 0;
      check(MPI_Reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_SUM, process, getComm(comm)));
      x = static_cast<UInt>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      check(MPI_Reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_SUM, process, getComm(comm)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      Double tmp = x;
      x = 0;
      check(MPI_Reduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, process, getComm(comm)));
    }
    else
      check(MPI_Reduce(&x, nullptr, 1, MPI_DOUBLE, MPI_SUM, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int >(x);
      unsigned int y = 0;
      check(MPI_Reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_SUM, process, getComm(comm)));
      x = static_cast<Bool>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      check(MPI_Reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_SUM, process, getComm(comm)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    constexpr UInt BLOCKSIZE = 50*1024*1024/sizeof(Double); // 50 Mb

    UInt index = 0;
    while(index<x.size())
    {
      const UInt size = std::min(x.size()-index, BLOCKSIZE);
      if(myRank(comm) == process)
      {
        Vector tmp(size);
        check(MPI_Reduce(x.field()+index, tmp.field(), size, MPI_DOUBLE, MPI_SUM, process, getComm(comm)));
        std::copy_n(tmp.field(), size, x.field()+index);
      }
      else
        check(MPI_Reduce(x.field()+index, nullptr, size, MPI_DOUBLE, MPI_SUM, process, getComm(comm)));

      index += size;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void reduceMin(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int>(x);
      unsigned int y   = tmp;
      check(MPI_Reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_MIN, process, getComm(comm)));
      x = static_cast<UInt>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      check(MPI_Reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_MIN, process, getComm(comm)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceMin(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      Double tmp = x;
      check(MPI_Reduce(&tmp, &x, 1, MPI_DOUBLE, MPI_MIN, process, getComm(comm)));
    }
    else
      check(MPI_Reduce(&x, nullptr, 1, MPI_DOUBLE, MPI_MIN, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceMax(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int>(x);
      unsigned int y   = tmp;
      check(MPI_Reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_MAX, process, getComm(comm)));
      x = static_cast<UInt>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      check(MPI_Reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_MAX, process, getComm(comm)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceMax(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      Double tmp = x;
      check(MPI_Reduce(&tmp, &x, 1, MPI_DOUBLE, MPI_MAX, process, getComm(comm)));
    }
    else
      check(MPI_Reduce(&x, nullptr, 1, MPI_DOUBLE, MPI_MAX, process, getComm(comm)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

} // namespace Parallel
