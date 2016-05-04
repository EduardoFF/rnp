/***********************************************************************************/
/***  Utility functions for socket communications                    ***/
/*    Author: Gianni Di Caro, IDSIA                                    */
/*                                                                     */
/***********************************************************************************/

#include "socket_lib.h"

int OpenBlockingUDPSocketServerSide(Socket &new_socket, const int port)
{
  if( (new_socket.sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    error("Opening socket");
  
  new_socket.addr_len = sizeof(struct sockaddr_in);
  bzero(&new_socket.server, new_socket.addr_len);
  
  new_socket.server.sin_family      = AF_INET;
  new_socket.server.sin_addr.s_addr = INADDR_ANY;
  new_socket.server.sin_port        = htons(port);
  int tr = 1;
  if (setsockopt(new_socket.sock,SOL_SOCKET,SO_REUSEADDR,&tr,sizeof(int)) == -1) {
    error("setsockopt");
    exit(1);
  }

  if ( bind(new_socket.sock,(struct sockaddr *)&new_socket.server, new_socket.addr_len) < 0 ) 
    error("Binding");
 
  return new_socket.sock;
}


int 
OpenNonBlockingUDPSocketServerSide(Socket &new_socket, const int port)
{
  if( (new_socket.sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    error("Opening socket");
  
  new_socket.addr_len = sizeof(struct sockaddr_in);
  bzero(&new_socket.server, new_socket.addr_len);
  
  new_socket.server.sin_family      = AF_INET;
  new_socket.server.sin_addr.s_addr = INADDR_ANY;
  new_socket.server.sin_port        = htons(port);
  int tr = 1;
  if (setsockopt(new_socket.sock,SOL_SOCKET,SO_REUSEADDR,&tr,sizeof(int)) == -1) {
    error("setsockopt");
    exit(1);
  }


  int flags = fcntl(new_socket.sock, F_GETFL, 0);
  if ( fcntl(new_socket.sock, F_SETFL, flags | O_NONBLOCK) )
  {
    error("fcntl");
  } else {
    printf("Sockets set NONBLOCK\n");
  }
  if ( bind(new_socket.sock,(struct sockaddr *)&new_socket.server, new_socket.addr_len) < 0 ) 
    error("Binding");
 
  return new_socket.sock;
}


int OpenUnixSocketPair(Socket &socket_a, Socket &socket_b)
{
	int sv[2];
	if (socketpair(AF_UNIX, SOCK_STREAM, 0, sv) == -1) {
		error("socketpair");
	}
	socket_a.sock = sv[0];
	socket_b.sock = sv[1];
}


int OpenUnixNonBlockingSocketPair(Socket &socket_a, Socket &socket_b)
{
	int sv[2];
	if (socketpair(AF_UNIX, SOCK_STREAM , 0, sv) == -1) {
		error("socketpair");
	}
	int flags = fcntl(sv[0], F_GETFL, 0);
	if ( fcntl(sv[0], F_SETFL, flags | O_NONBLOCK) )
	{
	  error("fcntl");
	} else {
	  printf("Sockets set NONBLOCK\n");
	}


	flags = fcntl(sv[1], F_GETFL, 0);
	if( fcntl(sv[1], F_SETFL, flags | O_NONBLOCK) )
	{
	  error("fcntl");
	} else 
	{
	  printf("Sockets set NONBLOCK\n");
	}

	///Mar 20, 2014 
	/// Get rid of the annoying "Resource temporary unavailable" error
	/// due to buffer overflow
	/// NOTE: It seems that this might not work, because the socket takes a
	/// default size equals (almost) to the maximum
	/// In ARES: is equal to the maximum /proc/sys/net/core/wmem_default
	//int sendbuff;
	//socklen_t optlen = sizeof(sendbuff);
	//int res = getsockopt(sv[0], SOL_SOCKET, SO_SNDBUF, &sendbuff, &optlen);

	//if(res == -1)
	  //printf("Error getsockopt one");
	//else
	  //printf("send buffer size = %d\n", sendbuff);

	///// Set buffer size
	//sendbuff = 131071;

	//printf("sets the send buffer to %d\n", sendbuff);
	//res = setsockopt(sockfd, SOL_SOCKET, SO_SNDBUF, &sendbuff, sizeof(sendbuff));

	//if(res == -1)
	  //printf("Error setsockopt");

	socket_a.sock = sv[0];
	socket_b.sock = sv[1];
}


int OpenUnixSocketServerSide(Socket &new_socket, const char *server_path)
{
       /********************************************************************/
      /* The socket() function returns a socket descriptor, which represents   */
      /* an endpoint.  The statement also identifies that the UNIX        */
      /* address family with the stream transport (SOCK_STREAM) will be   */
      /* used for this socket.                                            */
      /********************************************************************/
      if((new_socket.sock = socket(AF_UNIX, SOCK_DGRAM, 0)) < 0)
      {
         error("socket() failed");
 
      }

      new_socket.addr_len = sizeof(struct sockaddr_un);
      bzero(&new_socket.client_un, new_socket.addr_len);

      /********************************************************************/
      /* After the socket descriptor is created, a bind() function gets a */
      /* unique name for the socket.                                      */
      /********************************************************************/

      new_socket.addr_len =  strlen(server_path) + sizeof(new_socket.client_un.sun_family);
      new_socket.client_un.sun_family = AF_UNIX;
      strcpy(new_socket.client_un.sun_path, server_path);


      unlink(new_socket.client_un.sun_path);   /* in case it already exists */

  

      if(bind(new_socket.sock, (struct sockaddr *)&new_socket.client_un, new_socket.addr_len) < 0)
      {
         error("bind() failed");
      }
      return new_socket.sock;
}

int OpenUnixSocketClientSide(Socket &new_socket, const char *server_path)
{

	/********************************************************************/
	/* The socket() function returns a socket descriptor, which represents   */
	/* an endpoint.  The statement also identifies that the UNIX  */
	/* address family with the stream transport (SOCK_STREAM) will be   */
	/* used for this socket.                                            */
	/********************************************************************/
	if((new_socket.sock = socket(AF_UNIX, SOCK_DGRAM, 0)) < 0)
	{
		error("socket() failed");
	}
	new_socket.addr_len = sizeof(struct sockaddr_un);
	bzero(&new_socket.server_un, new_socket.addr_len);

	/********************************************************************/
	/* If an argument was passed in, use this as the server, otherwise  */
	/* use the #define that is located at the top of this program.      */
	/********************************************************************/
	new_socket.server_un.sun_family = AF_UNIX;
	strcpy(new_socket.server_un.sun_path, server_path);
	return new_socket.sock;
}




int OpenBlockingUDPSocketClientSide(Socket &new_socket, const int port, const char *server_host)
{
  if( (new_socket.sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    error("Opening socket");
  
  new_socket.server.sin_family = AF_INET;

  if ( (new_socket.hp = gethostbyname(server_host)) == 0 )
    error("Unknown host");
  
  bcopy((char *)new_socket.hp->h_addr, (char *)&new_socket.server.sin_addr, new_socket.hp->h_length);
  new_socket.server.sin_port = htons(port);

  new_socket.addr_len = sizeof(struct sockaddr_in);

  return new_socket.sock;
}


int ServerSendDataToSocket(Socket &socket, const char *msg, const int msg_size)
{
  //fprintf(stderr, "Server sending %d bytes to socket %d\n", msg_size, socket.sock);
  int r = sendto(socket.sock, msg, msg_size, 0, (struct sockaddr *)&socket.client, socket.addr_len);
  if( r != msg_size )
  {
    printf("sendto returned %d\n", r);
    error("Server send");
  }
  return r;
}

int SendDataToUnixSocket(Socket &socket, const char *msg, const int msg_size)
{
  //fprintf(stderr, "Server sending %d bytes to socket %d\n", msg_size, socket.sock);
  printf("sending %d bytes\n", msg_size);
  int bytes_sent = 0;
  int err = 0;
  while( bytes_sent < msg_size )
  {
    printf("... %d bytes\n", bytes_sent);
    fflush(stdout);
    int r = send(socket.sock, msg + bytes_sent, msg_size - bytes_sent, 0);
    if( r < 0 )
    {
      fprintf(stderr, 
	      "send failed with error %d errno %d and error %s\n",
	      r,
	      errno,
	      strerror(errno));
      err = errno;
      break;
    }
    bytes_sent += r;
  }
  printf("Done with socket send\n");
  fflush(stdout);
  return bytes_sent;
}


int ClientSendDataToSocket(Socket &socket, const char *msg, const int msg_size)
{
  //fprintf(stderr, "Client sending %d bytes to socket %d\n", msg_size, socket.sock);
  int r = sendto(socket.sock, msg, msg_size, 0, (struct sockaddr *)&socket.server, socket.addr_len);
  if( r != msg_size )
    error("Client send");
  return r;
}


int ServerReceiveDataFromSocketNonBlocking(Socket &socket, char **msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Server On receive max %d bytes from socket %d\n", max_data, socket.sock);

  msg_size = recvfrom(socket.sock, *msg, max_data, MSG_DONTWAIT, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size <= 0)
    msg_size = 0;
  else
    {
      //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
    }

  return msg_size;
}

int ServerReceiveDataFromSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Server On receive max %d bytes from socket %d\n", max_data, socket.sock);

  msg_size = recvfrom(socket.sock, msg, max_data, MSG_DONTWAIT, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size <= 0)
    msg_size = 0;
  else
    {
      //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
    }

  return msg_size;
}
int ReceiveDataFromUnixSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Server On receive max %d bytes from socket %d\n", max_data, socket.sock);

  msg_size = recv(socket.sock, msg, max_data, MSG_DONTWAIT);

  if (msg_size <= 0)
    msg_size = 0;
  else
    {
      //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
    }

  return msg_size;
}


int ServerPeekDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Server peek max %d bytes from socket %d\n", max_data, socket.sock);

  msg_size = recvfrom(socket.sock, msg, max_data, MSG_PEEK, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size < 0)
    error("Server Recvfrom peek");

  //fprintf(stderr, "[Peek] Received %d bytes from socket %d\n", msg_size, socket.sock);

  return msg_size;
}


int ServerPeekDataFromSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Server peek max %d bytes from socket %d\n", max_data, socket.sock);

  msg_size = recvfrom(socket.sock, msg, max_data, MSG_PEEK | MSG_DONTWAIT, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size <= 0)
    msg_size = 0;
  else
    {
      //fprintf(stderr, "[Peek] Received %d bytes from socket %d\n", msg_size, socket.sock);
    }

  return msg_size;
}



int ServerReceiveFixedDataLenFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int data_len)
{
  //fprintf(stderr, "Server On receive %d bytes from socket %d\n", data_len, socket.sock);

  msg_size = recvfrom(socket.sock, msg, data_len, MSG_WAITALL, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size < 0)
    error("Server Recvfrom fixed");

  //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);

  return msg_size;
}


int ServerReceiveDataFromSocketBlocking(Socket &socket, char **msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Server On receive max %d bytes from socket %d\n", max_data, socket.sock);

  msg_size = recvfrom(socket.sock, *msg, max_data, 0, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size < 0) 
    error("Server Recvfrom");

  //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
  return msg_size;
}


int ServerReceiveDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{

  //fprintf(stderr, "On receive max %d bytes from socket %d\n", max_data, socket.sock);


  msg_size = recvfrom(socket.sock, msg, max_data, 0, (struct sockaddr *)&socket.client, (socklen_t *)&socket.addr_len);

  if (msg_size < 0) 
    error("Server Recvfrom");

  //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
  return msg_size;
}


int ClientReceiveDataFromSocketNonBlocking(Socket &socket, char **msg, int &msg_size, const int max_data)
{

  //fprintf(stderr, "Client On receive max %d bytes from socket %d\n", max_data, socket.sock);

  struct sockaddr from;

  msg_size = recvfrom(socket.sock, *msg, max_data, MSG_DONTWAIT, (struct sockaddr *)&from, (socklen_t *)&socket.addr_len);

  if (msg_size <= 0)
    msg_size = 0;
  else
    {
      //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
    }

  return msg_size;
}


int ClientReceiveDataFromSocketNonBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{

  //fprintf(stderr, "Client On receive max %d bytes from socket %d\n", max_data, socket.sock);

  struct sockaddr from;

  msg_size = recvfrom(socket.sock, msg, max_data, MSG_DONTWAIT, (struct sockaddr *)&from, (socklen_t *)&socket.addr_len);

  if (msg_size <= 0)
    msg_size = 0;
  else
    {
      //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
    }

  return msg_size;
}


int ClientReceiveFixedDataLenFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int data_len)
{
  //fprintf(stderr, "Client On receive %d bytes from socket %d\n", data_len, socket.sock);

  struct sockaddr from;

  msg_size = recvfrom(socket.sock, msg, data_len, MSG_WAITALL, (struct sockaddr *)&from, (socklen_t *)&socket.addr_len);

  if (msg_size < 0)
    error("Client Recvfrom fixed");

  //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);

  return msg_size;
}


int ClientReceiveDataFromSocketBlocking(Socket &socket, char **msg, int &msg_size, const int max_data)
{
  //fprintf(stderr, "Client On receive max %d bytes from socket %d\n", max_data, socket.sock);

  struct sockaddr from;

  msg_size = recvfrom(socket.sock, *msg, max_data, 0, (struct sockaddr *)&from, (socklen_t *)&socket.addr_len);

  if (msg_size < 0) 
    error("Client Recvfrom");

  //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
  return msg_size;
}


int ClientReceiveDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{

  //fprintf(stderr, "Client On receive max %d bytes from socket %d\n", max_data, socket.sock);

  struct sockaddr from;

  msg_size = recvfrom(socket.sock, msg, max_data, 0, (struct sockaddr *)&from, (socklen_t *)&socket.addr_len);

  if (msg_size < 0) 
    error("Client Recvfrom");

  //fprintf(stderr, "Received %d bytes from socket %d\n", msg_size, socket.sock);
  return msg_size;
}


int ClientPeekDataFromSocketBlocking(Socket &socket, char *msg, int &msg_size, const int max_data)
{

  //fprintf(stderr, "Client peek max %d bytes from socket %d\n", max_data, socket.sock);

  struct sockaddr from;

  msg_size = recvfrom(socket.sock, msg, max_data, MSG_PEEK, (struct sockaddr *)&from, (socklen_t *)&socket.addr_len);

  if (msg_size < 0) 
    error("Client Recvfrom");

  //fprintf(stderr, "[Peek] Received %d bytes from socket %d\n", msg_size, socket.sock);
  return msg_size;
}

void error(const char *msg)
{
    perror(msg);
    exit(-1);
}

