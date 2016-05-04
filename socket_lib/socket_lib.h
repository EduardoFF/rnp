/***********************************************************************************/
/***  Utility functions for socket communications                    ***/
/*    Author: Gianni Di Caro, IDSIA                                    */
/*                                                                     */
/***********************************************************************************/

#ifndef SOCKET_LIB_H 
#define SOCKET_LIB_H 

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <netdb.h>
#include <math.h>
#include <fcntl.h>
#include <errno.h>

using namespace std;

#define _SOCKETLIB_MAX_MSG_SIZE 65536

const int max_data = _SOCKETLIB_MAX_MSG_SIZE;


typedef float Real;
typedef signed short SInt16;
typedef unsigned short UInt16;

typedef signed int SInt32;
typedef unsigned int UInt32;

typedef signed long long SInt64;
typedef unsigned long long UInt64;


typedef struct
{
  int sock;
  int addr_len;
  struct sockaddr_in server;
  struct sockaddr_in client;
  struct sockaddr_un server_un;
  struct sockaddr_un client_un;
  struct hostent *hp;

} Socket;



/***********************************************************************************/
/***  Common functions for socket communications                    ***/
/*                                                                    */
/***********************************************************************************/

/// OpenBlockingUDPSocket simply means that we do not set the O_NONBLOCK flag
/// a non-blocking operation can still be achieved using the MSG_DONTWAIT flag
/// when receiving data, using the ReceivedData*NonBlocking functions below
int OpenBlockingUDPSocketServerSide(Socket &new_socket, const int port);

int OpenNonBlockingUDPSocketServerSide(Socket &new_socket, const int port);

int OpenBlockingUDPSocketClientSide(Socket &new_socket, const int port, const char *server_host);

int OpenUnixSocketServerSide(Socket &new_socket, const char *server_path);

int OpenUnixSocketClientSide(Socket &new_socket, const char *server_path);


int OpenUnixNonBlockingSocketPair(Socket &socket_a, Socket &socket_b);
int OpenUnixSocketPair(Socket &, Socket &);

int SendDataToUnixSocket(Socket &socket, const char *msg, const int msg_size);

//int ClientSendDataToUnixSocket(Socket &socket, const char *msg, const int msg_size);


int ServerSendDataToSocket(Socket &socket, const char *msg, const int msg_size);

int ClientSendDataToSocket(Socket &socket, const char *msg, const int msg_size);

int ReceiveDataFromUnixSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ServerReceiveDataFromSocketNonBlocking(Socket &socket, char **msg, int &msg_size, const int max_data);
int ServerReceiveDataFromSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ServerReceiveDataFromSocketBlocking(Socket &socket, char **msg, int &msg_size, const int max_data);
int ServerReceiveDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ServerPeekDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);
int ServerPeekDataFromSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ServerReceiveFixedDataLenFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int data_len);

int ClientReceiveDataFromSocketNonBlocking(Socket &socket, char **msg, int &msg_size, const int max_data);
int ClientReceiveDataFromSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ClientReceiveDataFromSocketBlocking(Socket &socket, char **msg, int &msg_size, const int max_data);
int ClientReceiveDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ClientPeekDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data);

int ClientReceiveFixedDataLenFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int data_len);

void error(const char *msg);

#endif
