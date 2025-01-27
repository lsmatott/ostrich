/******************************************************************************
File      : mpi_stub.h
Author    : L. Shawn Matott
Copyright : 2002, L. Shawn Matott

Stubs out mpi functions, so that one can compile the Ostrich code without using
the MPI libraries.

Version History
11-18-02    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
******************************************************************************/
#ifndef USE_MPI_STUB
#include <mpi.h>
#else
#ifndef MPI_INCLUDED
#define MPI_INCLUDED

/* Keep C++ compilers from getting confused */
#ifdef __cplusplus
extern "C"
{
#endif
 
 typedef int MPI_Comm;
 typedef int MPI_Datatype;
 typedef int MPI_Op;

 typedef struct MPI_STATUS 
 {
	 int MPI_SOURCE;
	 int MPI_TAG;
 }MPI_Status;

 #define MPI_SUM 0
 #define MPI_DOUBLE 0 
 #define MPI_INTEGER 0
 #define MPI_INT 0
 #define MPI_CHAR 0
 #define MPI_ANY_SOURCE 0
 #define MPI_ANY_TAG 0
 #define MPI_COMM_WORLD 91

 int MPI_Init(int * argc, char *** argv);
 int MPI_Abort(MPI_Comm comm, int errorcode);

 int MPI_Comm_size(MPI_Comm comm, int * size);
 int MPI_Comm_rank(MPI_Comm comm, int * rank);

 int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                void *recvbuf, int *recvcnts, int *displs, 
                MPI_Datatype recvtype, int root, MPI_Comm comm);

 int MPI_Allgatherv (void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcounts, int *displs, 
                     MPI_Datatype recvtype, MPI_Comm comm );

int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm);

 int MPI_Barrier (MPI_Comm comm);

 int MPI_Bcast(void * buf, int count, MPI_Datatype datatype, int root, 
	            MPI_Comm comm);

 int MPI_Reduce(void * sendbuf, void * recvbuf, int count, 
	             MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

 int MPI_Allreduce(void * sendbuf, void * recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

 int MPI_Recv(void * buf, int count, MPI_Datatype datatype, int source, 
	           int tag, MPI_Comm comm, MPI_Status * status);

 int MPI_Send(void * buf, int count, MPI_Datatype datatype, int dest, int tag, 
	           MPI_Comm comm);

 int MPI_Finalize(void);
#ifdef __cplusplus
}
#endif

#endif /* MPI_INCLUDED */
#endif /* USE_MPI_STUB */
