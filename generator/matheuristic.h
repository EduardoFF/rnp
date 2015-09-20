#ifndef _MATHEURISTIC_H
#define _MATHEURISTIC_H

#include "socket_lib.h"
#include "fcntl.h"
#include <sys/msg.h>
#include "rnpsolution.h"
typedef enum 
{
  TIME_INFO,
  TEST,
  SOLUTION_PACKET,
} MsgType;


typedef struct
{
  MsgType  type;
  UInt32   size;
  Real     timestamp;
  UInt32   seqn;

} HeaderPacket;

typedef struct
{
  int time;
} TimeInfoPacket;


typedef struct
{
  double testvalue;
} TestPacket;


#define MSGQUEUEMAX 8000
#define MAX_PAYLOAD_SIZE 32768
typedef struct
{
  int len;
  char data[MAX_PAYLOAD_SIZE];
} SolutionPacket;

typedef struct  {
  long msgtype;
  char    mtext[MSGQUEUEMAX];
} message_buf;
typedef enum
{
  COMM_PIPES,
  COMM_SOCKETS,
  COMM_MSGQUEUE,
  COMM_UDP
} CommType;
class Matheuristic
{
  CommType m_comm;


  Socket m_clientSocket;
  Socket m_serverSocket;

  int m_CSPipe[2]; // Client -> Server
  int m_SCPipe[2]; // Server -> Client
  int m_readPipe;
  int m_sendPipe;

  /// For MSG_QUEUE
  int m_SCqueue, m_CSqueue;

  FILE *m_readStream, *m_writeStream;
  bool m_amServer;
  int m_lastSync;
  int m_nextSync;
  int m_seqn;
  int m_lastRcvSeqn;
  double m_lastPushObjval;
  double m_lastPullObjval;
  double m_lastRcvTest, m_lastSentTest;
  CpuTime m_cpuTime;
  ofstream m_output;
  bool m_outputRequested;

  int sendSolution(LpSolutionPtr );
  public:
  Matheuristic();
  void setServer(bool b);
  void setOutput(string fname);
  int pushSolution(LpSolutionPtr);
  LpSolutionPtr pullSolution(bool *);
  int testPush(double);
  int testPull(double*);
  double lastPull(){ return m_lastPullObjval;}
  double lastPush(){ return m_lastPushObjval;}


};



#endif
