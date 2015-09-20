#include "main.h"


#include "socket_lib.h"

char rcv_buffer[_SOCKETLIB_MAX_MSG_SIZE];
char send_buffer[_SOCKETLIB_MAX_MSG_SIZE];

const int MAX_SOCKETS = 64;
int main ( int argc, char *argv[] )
{
  char client_msg[_SOCKETLIB_MAX_MSG_SIZE];
  char client_active[_SOCKETLIB_MAX_MSG_SIZE];
  int client_msg_size;
  bool ok;

  int n_sockets = 3;
  int base_port = 12300;
  Socket socket[MAX_SOCKETS];
  for(int i=0; i< n_sockets; i++)
  {
    OpenNonBlockingUDPSocketServerSide(socket[i], base_port + i);
    client_active[i] = false;
  }
  
  
  while( true ) 
  {
    for(int i=0; i< n_sockets; i++)
    {
      if( ServerReceiveDataFromSocketNonBlocking(socket[i], client_msg, client_msg_size, _SOCKETLIB_MAX_MSG_SIZE) )
      {
	printf("Received msg in socket %d with len %d\n",i, client_msg_size);
	client_active[i] = true;
	for(int j=0; j< n_sockets; j++)
	{
	  if( j==i || !client_active[j] )
	    continue;
	  printf("Trying to send to %d\n", j);
	  ok = ( ServerSendDataToSocket(socket[j], client_msg, client_msg_size) == client_msg_size);
	  if(ok)
	  {
	    printf("Packet sent successfully to %d!\n", j);
	  }
	  else
	  {
	    printf("ERROR Sending packet to %d\n", j);
	  }
	}
      }
    }
  }
  return 0;
}









