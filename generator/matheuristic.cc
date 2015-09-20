#include "main.h"
#include "generator.h"
#include "solver.h"
#include "graph.h"
#include "network.h"
#include "metric.h"
#include "matheuristic.h"

#include <ext/stdio_filebuf.h>
#include <errno.h>
//#inclaude "stdio_filebuf.h"
#include "socket_lib.h"
extern param prob;

char rcv_buffer[_SOCKETLIB_MAX_MSG_SIZE];
char send_buffer[_SOCKETLIB_MAX_MSG_SIZE];

char msgq_buffer[MSGQUEUEMAX];
string read_from_pipe (int file_descriptor)
{
  // We create a C++ stream from a file descriptor
  // stdio_filebuf is not synced with stdio.
  // From GCC 3.4.0 on exists in addition stdio_sync_filebuf
  // You can also create the filebuf from a FILE* with
  // FILE* f = fdopen(file_descriptor, mode);
  __gnu_cxx::stdio_filebuf<char> filebuf(file_descriptor,
					 std::ios_base::in);
  istream stream(&filebuf);
  // You can also do:
  // ostringstream stream;
  // stream << &filebuf;
  // return stream.str();

  string line;
  if (stream.good())
    getline(stream, line);

  return line;
}
void write_to_pipe (int file_descriptor, const string& line)
{
  __gnu_cxx::stdio_filebuf<char> filebuf(file_descriptor,
					 std::ios_base::out);
  ostream stream(&filebuf);

  if (stream.good())
    stream << line;
}

void createMsgHeader(HeaderPacket &p_header,MsgType p_msgType, int p_size, int p_seqn)
{
  p_header.type = p_msgType;
  p_header.size = p_size;
  p_header.seqn = p_seqn;
}

  int 
Matheuristic::testPull(double *val)
{
  bool gotSomething = false;

  if( m_comm == COMM_SOCKETS 
      || m_comm == COMM_UDP)
  {
    int msg_size;
    char rcv_msg_packet[_SOCKETLIB_MAX_MSG_SIZE];

    if(m_amServer && m_comm == COMM_SOCKETS)
    {
      gotSomething = (ReceiveDataFromUnixSocketNonBlocking(m_serverSocket,rcv_msg_packet, msg_size, _SOCKETLIB_MAX_MSG_SIZE));
    }
    else
    {
      if( m_comm == COMM_SOCKETS )
      {
      	gotSomething = (ReceiveDataFromUnixSocketNonBlocking(m_clientSocket,rcv_msg_packet, msg_size, _SOCKETLIB_MAX_MSG_SIZE));
      } 
      else
      {
	gotSomething = (ClientReceiveDataFromSocketBlocking(m_clientSocket, rcv_msg_packet, msg_size, _SOCKETLIB_MAX_MSG_SIZE));
      }
    }

    if(gotSomething)
    {
      // Got something from server
      printf("Got something from server\n");
      HeaderPacket *msg = (HeaderPacket *)rcv_msg_packet;
      //int seqn = msg->seqn;
      //printf("Seqno is %d\n", seqn);
      //if( seqn <= m_lastRcvSeqn)
      //{
	//printf("Last received seqno is %d > %d - return 0\n",
	       //m_lastRcvSeqn, seqn);
	//return 0;
      //}
      //m_lastRcvSeqn = seqn;
      if(msg->type == TEST)
      {

	TestPacket *msgInfo = 
	  (TestPacket *)(rcv_msg_packet + sizeof(HeaderPacket));
	printf("Test value received is %f\n",
	       msgInfo->testvalue);
	*val =  msgInfo->testvalue;
	return 1;
      }
    }
    return 0;
  } 
  else if(m_comm == COMM_PIPES)
  {
    //		string gots = read_from_pipe(m_readPipe);
    //		cout << gots << endl;

    char *line = NULL;
    size_t len = 0;
    ssize_t b_read;

    FILE *fp = m_readStream;
    if (fp == NULL)
      exit(EXIT_FAILURE);
    //		while(read(m_readPipe, &buffer, 1)>0)
    //		{

    //		b_read = read(m_readPipe, &buffer, 1);
    //		if( b_read < 0)
    //			printf("Error reading!\n");
    //		else if(b_read == 0)
    //			printf("No DATA\n");
    //		else
    //			printf("got char %c\n",buffer[0]);
    //		}

    while(fgets (rcv_buffer, sizeof (rcv_buffer), fp) != NULL)
    {
      printf("Retrieved line of length %zu :\n", strlen(rcv_buffer));
      printf("%s", rcv_buffer);
    }
    //		while ((read = getline(&line, &len, fp)) != -1) {
    //			printf("Retrieved line of length %zu :\n", read);
    //			printf("%s", line);
    //		}
    if( line)
      free(line);



  }


}


#define SERIALIZE_BINARY 1
#define SOCKET_JITTER 250
  int 
Matheuristic::sendSolution(LpSolutionPtr lpsol)
{

  int wait_ms = 200;
  int max_trials = 5; /// 1 second wait
  int buf_len = 65536;
  int msg_size;
  string toSend;
  std::vector<char> buff_vec;

  printf("Attempting to push LpSolution with obj val %f\n", lpsol->objval);
  if( SERIALIZE_BINARY)
  {
    save_lpsolution(*lpsol, buff_vec, buf_len);
    printf("Serialized BINARY size %d\n", buff_vec.size());
    msg_size = buff_vec.size();
    memcpy(send_buffer, buff_vec.data(), msg_size);

  } else
  {
    save_lpsolution(*lpsol, toSend);
    msg_size = toSend.size()+1;
    memcpy(send_buffer, toSend.c_str(), msg_size-1);
    send_buffer[ msg_size] = '\0';
    printf("Serialized TEXT size %d\n", toSend.size());
  }

  if( msg_size >= _SOCKETLIB_MAX_MSG_SIZE)
  {
    fprintf(stderr, "Solution size too big. %d (%d)\n", msg_size, _SOCKETLIB_MAX_MSG_SIZE);
    if(m_outputRequested)
    {
      m_output << "PUSH_ERROR " << m_cpuTime.cpu_time_elapsed() << " " 
	<< lpsol->objval << endl;
    }
    return 1;
  }

  if( m_comm == COMM_PIPES)
  {
    printf("Matheuristic: SENDING  through PIPE. msg_size %d\n", msg_size);
    int bytes_written = 0;
    int trials = 0;
    while( bytes_written < msg_size)
    {
      int remaining = msg_size - bytes_written;
      int written = write(m_sendPipe, (void *) (send_buffer + bytes_written), remaining);


      if( written < 0)
      {
	//fprintf(stderr, "PIPE write errno %s\n", strerror(errno));
	if ( errno == EAGAIN)
	{
	  trials++;
	  if(trials < max_trials)
	  {
	    usleep(wait_ms*1000);
	    //printf("Got EAGAIN\n");
	    continue;
	  }else
	    fprintf(stderr, "PIPE max trials\n");

	}

	printf("ERROR Sending to PIPE\n");

	fprintf(stderr,"ERROR Sending to PIPE\n");


	if(m_outputRequested)
	{
	  m_output << "PUSH_ERROR " << m_cpuTime.cpu_time_elapsed() << " " 
	    << lpsol->objval << endl;
	}
	return 1;
      }

      bytes_written += written;
      printf("PIPE: Written %d (%d) - trials %d\n", written, bytes_written, trials);
    }
  } 
  else if( m_comm == COMM_MSGQUEUE)
  {
    /// MSG QUEUE NOT WORKING PROPERLY

    message_buf my_msg;
    int buf_length;
    int msg_size = toSend.size()+1;
    buf_length = MSGQUEUEMAX;
    memcpy(my_msg.mtext, toSend.c_str(), msg_size-1);
    my_msg.msgtype = 2;
    my_msg.mtext[ msg_size] = '\0';
    int msqid;
    if( m_amServer)
    {
      msqid = m_SCqueue;
    } else 
    {
      msqid = m_CSqueue;
    }

    if (msgsnd(msqid, &my_msg, buf_length, IPC_NOWAIT) < 0) {
      //printf ("%d, %d, %s, %d\n", msqid, sbuf.mtype, sbuf.mtext, buf_length);
      perror("msgsnd");
      exit(1);
    }
    else 
      printf("Message Sent\n");
  } else if( m_comm == COMM_SOCKETS ||
	     m_comm == COMM_UDP)
  {

    char send_msg_packet[_SOCKETLIB_MAX_MSG_SIZE];
    bool ok;

    HeaderPacket mymsg;
    SolutionPacket content;
    int content_size = sizeof(SolutionPacket);

    createMsgHeader(mymsg,SOLUTION_PACKET,content_size, m_seqn);
    memcpy(content.data, send_buffer, msg_size);
    content.len = msg_size;


    int total_size = sizeof(HeaderPacket) + content_size;
    memcpy( send_msg_packet, (void *)&mymsg, sizeof(HeaderPacket));
    memcpy( send_msg_packet + sizeof(HeaderPacket),(void *)&content, content_size);

    printf("Matheuristic: SENDING though SOCKETS (total_size %d) \n", total_size);
    if(m_amServer && m_comm == COMM_SOCKETS)
    {
      ok = ( SendDataToUnixSocket(m_serverSocket, send_msg_packet, total_size) == total_size);
    } 
    else 
    {
      if( m_comm == COMM_SOCKETS )
      {
      	ok = ( SendDataToUnixSocket(m_clientSocket, send_msg_packet, total_size) == total_size);
      }
      else if(m_comm == COMM_UDP)
      {
      	ok = ( ClientSendDataToSocket(m_clientSocket, send_msg_packet, total_size) == total_size);
      }
    }
    if(ok)
    {
      printf("SOLUTION_PACKET sent successful (%d)!\n", m_seqn);
      m_seqn++;
      /// Flow control, sleep for 250ms
      usleep(SOCKET_JITTER*1000);
      //			m_lastSentTest = val;
    }else
      printf("Error sending TEST packet\n");
  }


  if(m_outputRequested)
  {
    m_output << "PUSH_OK " << m_cpuTime.cpu_time_elapsed() << " " 
      << lpsol->objval << endl;
  }
  return 0;
}


  int
Matheuristic::pushSolution(LpSolutionPtr lpsol)
{
  m_lastPushObjval = lpsol->objval;
  //m_sendingQueue.push_back(lpsol);
  //while(
  return sendSolution(lpsol);

}


  LpSolutionPtr 
Matheuristic::pullSolution(bool *valid)
{
  bool gotData = false;
  int data_len;
  if(m_comm == COMM_PIPES)
  {

    FILE *fp = m_readStream;
    if (fp == NULL)
    {
      printf("Matheuristic: There was an ERROR\n");
      exit(EXIT_FAILURE);
    }
    if(fgets (rcv_buffer, sizeof (rcv_buffer), fp) != NULL)
    {
      printf("Retrieved line of length %zu :\n", strlen(rcv_buffer));
      gotData = true;
      data_len = strlen(rcv_buffer);	
    }
  } else if( m_comm == COMM_MSGQUEUE)
  {
    message_buf my_msg;
    int msqid;
    if( m_amServer)
    {
      msqid = m_CSqueue;
    } else 
    {
      msqid = m_SCqueue;
    }

    if (msgrcv(msqid, (void *)&my_msg, sizeof(message_buf), 2, IPC_NOWAIT) < 0){
      //printf ("%d, %d, %s, %d\n", msqid, sbuf.mtype, sbuf.mtext, buf_length);

      if ( errno == ENOMSG)
      {
	printf("No message avaliable\n");
      } else
      {
	perror("msgrcv");
	exit(1);
      }
    }
    else 
    {
      printf("Message Received len %d\n", strlen(my_msg.mtext));
      data_len = strlen(my_msg.mtext);
      memcpy( rcv_buffer, my_msg.mtext, data_len); 

    }
  } 
  else if( m_comm == COMM_SOCKETS ||
	   m_comm == COMM_UDP)
  {
    int msg_size;
    bool gotSomething = false;
    char rcv_msg_packet[_SOCKETLIB_MAX_MSG_SIZE];
    HeaderPacket header;
    SolutionPacket msgInfo;
    int max_msg_size = sizeof(HeaderPacket) + sizeof(SolutionPacket);
    if(m_amServer && m_comm == COMM_SOCKETS)
    {
      gotSomething = 
	(ReceiveDataFromUnixSocketNonBlocking(m_serverSocket,rcv_msg_packet, 
					      msg_size, max_msg_size ));
    }
    else
    {
      if( m_comm == COMM_SOCKETS )
      {
      	gotSomething = (ReceiveDataFromUnixSocketNonBlocking(m_clientSocket,rcv_msg_packet, 
							   msg_size, max_msg_size));
      }
      else if( m_comm == COMM_UDP)
      {
	gotSomething = (ClientReceiveDataFromSocketNonBlocking(m_clientSocket, rcv_msg_packet, msg_size, max_msg_size));
      }
    }

    if(gotSomething){
      // Got something from server
      printf("SOCKETS: Got something from server msg_size (%d)\n",
	     msg_size);
      memcpy(&header, (HeaderPacket *)rcv_msg_packet, sizeof(HeaderPacket));
      int seqn = header.seqn;
      printf("sequence number %d\n", 
	     seqn);

      m_lastRcvSeqn = seqn;
      if(header.type == SOLUTION_PACKET)
      {
	printf("Got SOLUTION PACKET\n");
	memcpy(&msgInfo, (SolutionPacket *)(rcv_msg_packet + sizeof(HeaderPacket)),
	       sizeof(SolutionPacket));
	printf("Solution length %d\n", msgInfo.len);

	memcpy(rcv_buffer, msgInfo.data, msgInfo.len);
	data_len = msgInfo.len;
	gotData = true;
      } else
      {
	printf("Got STRANGE PACKET\n");
      }
    } else {
      printf("SOCKETS: No message available\n");
    }
  }
  /// Here, we assume data is placed in rcv_buffer
  /// and variable data_len has the number of valid bytes
  // Return empty ptr
  *valid = false;
  LpSolution *lp_sol = new LpSolution();
  if( gotData)
  {
    printf("GOT DATA (data_len %d)\n", data_len);
    if( SERIALIZE_BINARY)
    {
      ///Currently, only works with SOCKETS
      restore_lpsolution(*lp_sol, rcv_buffer, data_len);
    } else
    {
      string lp_str(rcv_buffer);
      restore_lpsolution(*lp_sol, lp_str);
    }
    printf("Matheuristic: got sol with val %f\n",lp_sol->objval);	
    if(lp_sol->solved)
    {
      LpSolutionPtr lp_ptr(lp_sol);
      *valid = true;
      if(m_outputRequested)
      {
	m_output << "PULL_OK " << m_cpuTime.cpu_time_elapsed() << " " 
	  << lp_ptr->objval <<  " " << lp_ptr->solved << endl;
      }
      m_lastPullObjval = lp_sol->objval;
      return lp_ptr;
    }
    else
    {
      if(m_outputRequested)
      {
	m_output << "PULL_ERROR " << m_cpuTime.cpu_time_elapsed() << " 0" << endl;
      }
    }
  }
  LpSolutionPtr lp_ptr(lp_sol);

  return lp_ptr;
}

  int 
Matheuristic::testPush(double val)
{
  if( m_comm == COMM_SOCKETS  
      || m_comm == COMM_UDP )
  {
    char send_msg_packet[_SOCKETLIB_MAX_MSG_SIZE];
    bool ok;

    HeaderPacket mymsg;
    TestPacket content;
    int content_size = sizeof(TestPacket);

    createMsgHeader(mymsg,TEST,content_size, m_seqn);

    content.testvalue = val;

    int total_size = sizeof(HeaderPacket) + content_size;
    memcpy( send_msg_packet, (void *)&mymsg, sizeof(HeaderPacket));
    memcpy( send_msg_packet + sizeof(HeaderPacket),(void *)&content, content_size);

    printf("Matheurustic: Trying to send testvalue %f\n", val);
    if(m_amServer && m_comm == COMM_SOCKETS)
    {
      ok = ( SendDataToUnixSocket(m_serverSocket, send_msg_packet, total_size) == total_size);
    } 
    else 
    {
      if( m_comm == COMM_SOCKETS )
      {
      	ok = ( SendDataToUnixSocket(m_clientSocket, send_msg_packet, total_size) == total_size);
      }
      else
      {
      	ok = ( ClientSendDataToSocket(m_clientSocket, send_msg_packet, total_size) == total_size);
      }
    }
    if(ok)
    {
      printf("TEST packet sent successful!\n");
      m_seqn++;
      m_lastSentTest = val;
    }else
      printf("Error sending TEST packet\n");
  }
  if( m_comm == COMM_PIPES)
  {
    //		stringstream ss;
    //		ss << "Hello " << val << endl;

    sprintf(send_buffer, "Hello %f", val);
    int msg_size = strlen(send_buffer)+1;
    printf("Matheuristic: TEST through PIPE (%s) len %d\n", send_buffer, msg_size);
    if( write(m_sendPipe, send_buffer, msg_size) != msg_size)
    {
      printf("ERROR Sending to PIPE\n");
    }
    //		write_to_pipe(m_sendPipe, ss.str());
  }
}


void 
Matheuristic::setServer(bool b )
{ 
  m_amServer = b;
  if( m_comm == COMM_PIPES)
  {
    if(m_amServer)
    {
      m_sendPipe = m_SCPipe[1];
      m_readPipe = m_CSPipe[0];
      close(m_SCPipe[0]);
      close(m_CSPipe[1]);
    } else
    {

      m_sendPipe = m_CSPipe[1];
      m_readPipe = m_SCPipe[0];
      close(m_SCPipe[1]);
      close(m_CSPipe[0]);
    }
    m_readStream = fdopen (m_readPipe, "r");
    if(setvbuf(m_readStream, NULL, _IONBF, 65536))
      printf("Error setting inout buffer\n");
    //m_writeStream = fdopen(m_sendPipe, "w");

  } else if( m_comm == COMM_MSGQUEUE)
  {

  }
}

  void
Matheuristic::setOutput(string fname)
{
  m_outputRequested = true;
  m_output.open(fname.c_str(), ios_base::out | ios_base::trunc );
}


Matheuristic::Matheuristic()
{
  m_seqn = 0;
  m_lastRcvSeqn= -1;

  m_lastPullObjval = 1000000;
  m_lastPushObjval = 1000000;
  m_cpuTime;
  m_cpuTime.start();

  m_outputRequested = false;
  if( prob.matheuristic_params.comm == "pipes")
  {
    printf("Matheuristic: Opening Pipes\n");
    m_comm = COMM_PIPES;
    if(pipe(m_CSPipe))
    {
      printf("Error creating pipe\n");
      exit(1);
    } else
    {
      if(fcntl(m_CSPipe[0],F_SETFL,O_NONBLOCK ) < 0)
      {
	fprintf(stderr, "error setting Pipe\n");
	exit(1);
      }
      if(fcntl(m_CSPipe[1],F_SETFL,O_NONBLOCK ) < 0)
      {
	fprintf(stderr, "error setting Pipe\n");
	exit(1);
      }
    }

    if(pipe(m_SCPipe))
    {
      printf("Error creating Pipe\n");
      exit(1);
    } else
    {
      if(fcntl(m_SCPipe[0],F_SETFL,O_NONBLOCK ) < 0)
      {
	fprintf(stderr, "error setting Pipe\n");
	exit(1);
      }
      if(fcntl(m_SCPipe[1],F_SETFL,O_NONBLOCK ) < 0)
      {
	fprintf(stderr, "error setting Pipe\n");
	exit(1);
      }
    }


    /*
      if(pipe2(m_CSPipe,O_NONBLOCK)!=0)
      {
      printf("Error creating pipe\n");
      }

      if(pipe2(m_SCPipe,O_NONBLOCK)!=0)
      {
      printf("Error creating pipe\n");
      }
      */

  } else if(prob.matheuristic_params.comm == "sockets")
  {
    printf("Matheuristic: Opening sockets\n");
    m_comm = COMM_SOCKETS;
    OpenUnixSocketPair(m_clientSocket, m_serverSocket);
  }
  else if( prob.matheuristic_params.comm == "udp" )
  {
    m_comm = COMM_UDP;
    printf("Matheuristic: Connecting to %s:%d\n", 
	   prob.matheuristic_params.socket_address.c_str(),
	   prob.matheuristic_params.socket_port);
    OpenBlockingUDPSocketClientSide(m_clientSocket, 
				    prob.matheuristic_params.socket_port,
				    prob.matheuristic_params.socket_address.c_str());
  }
  else if( prob.matheuristic_params.comm == "msgqueue")
  {
    m_comm = COMM_MSGQUEUE;

    /// Create server -> client msg queue
    string msg_file = prob.output_path + prob.instance_id + ".msgqueue.1";
    printf("Creating msg queue file in %s\n", msg_file.c_str());
    key_t key = ftok(msg_file.c_str(), 'b');
    int msgflg = IPC_CREAT | 0666;
    if ((m_SCqueue = msgget(key, msgflg )) < 0) {
      perror("msgget");
      exit(1);
    } else
    {
      fprintf(stderr,"msgget: msgget succeeded: m_SCqueue = %d\n", m_SCqueue);
    }



    /// Create client -> server msg queue
    msg_file = prob.output_path + prob.instance_id + ".msgqueue.2";

    printf("Creating msg queue file in %s\n", msg_file.c_str());
    key = ftok(msg_file.c_str(), 'b');
    if ((m_CSqueue = msgget(key, msgflg )) < 0) {
      perror("msgget");
      exit(1);
    } else
    {
      fprintf(stderr,"msgget: msgget succeeded: m_CSqueue = %d\n", m_CSqueue);
    }



  }
  //	m_amServer = server;
  //	string address = "/tmp/myaddr";
  //if(server)
  //	{
  //		OpenUnixSocketServerSide(m_socket, address.c_str());
  //	}
  //	else
  //	{
  //		OpenUnixSocketClientSide(m_socket, address.c_str());
  //	}

}

#if 0

int Matheuristic::Synchronize(int time)
{


  const int sync_gap = 2;
  char send_msg_packet[_SOCKETLIB_MAX_MSG_SIZE];
  char rcv_msg_packet[_SOCKETLIB_MAX_MSG_SIZE];
  int ret_val=-1;


  int msg_size;
  //Server waits for client
  if(m_amServer)
  {
    HeaderPacket mymsg;
    TimeInfoPacket content;
    int content_size = sizeof(TimeInfoPacket);

    createMsgHeader(mymsg,TIME_INFO,content_size);

    content.time = time;

    int total_size = sizeof(HeaderPacket) + content_size;
    memcpy( send_msg_packet, (void *)&mymsg, sizeof(HeaderPacket));
    memcpy( send_msg_packet + sizeof(HeaderPacket),(void *)&content, content_size);

    printf("Server trying to send time %d\n", time);
    if( SendDataToUnixSocket(m_serverSocket, send_msg_packet, total_size) == total_size)
    {
      printf("Server sent OK\n");
    }
    else{
      printf("Server sent BAD\n");
    }

    // Now wait for client
    if(ReceiveDataFromUnixSocketNonBlocking(m_serverSocket,rcv_msg_packet, msg_size, MAX_MSG_SIZE))
    {
      printf("Got something from client\n");

      HeaderPacket *msg = (HeaderPacket *)rcv_msg_packet;
      if(msg->type == TIME_INFO)
      {

	TimeInfoPacket *msgInfo = 
	  (TimeInfoPacket *)(rcv_msg_packet + sizeof(HeaderPacket));
	printf("client time is %d my time is %d\n",
	       msgInfo->time, time);

	if(abs(msgInfo->time - time) < sync_gap)
	{
	  ret_val = 0;
	}
	else
	{
	  if( time > msgInfo->time ) // I'm ahead!
	  {
	    ret_val = 1;
	  }
	  else // I'm behind, I should hurry
	  {
	    ret_val = 0;
	  }
	}

      }

    }

  }
  else
  {


    if(ReceiveDataFromUnixSocketNonBlocking(m_clientSocket,rcv_msg_packet, msg_size, MAX_MSG_SIZE))
    {
      // Got something from server
      printf("Got something from server\n");
      HeaderPacket *msg = (HeaderPacket *)rcv_msg_packet;
      if(msg->type == TIME_INFO)
      {

	TimeInfoPacket *msgInfo = 
	  (TimeInfoPacket *)(rcv_msg_packet + sizeof(HeaderPacket));
	printf("server time is %d my time is %d\n",
	       msgInfo->time, time);

	if(abs(msgInfo->time - time) < sync_gap)
	{
	  ret_val = 0;
	}
	else
	{
	  if( time > msgInfo->time ) // I'm ahead!
	  {
	    ret_val = 1;
	  }
	  else // I'm behind, I should hurry
	  {
	    ret_val = 0;
	  }
	}

      }

    }


    HeaderPacket mymsg;
    TimeInfoPacket content;
    int content_size = sizeof(TimeInfoPacket);

    createMsgHeader(mymsg,TIME_INFO,content_size);

    content.time = time;

    int total_size = sizeof(HeaderPacket) + content_size;
    memcpy( send_msg_packet, (void *)&mymsg, sizeof(HeaderPacket));

    memcpy( send_msg_packet + sizeof(HeaderPacket),(void *)&content, content_size);
    printf("Client trying to send %d\n", time);
    if( SendDataToUnixSocket(m_clientSocket, send_msg_packet, total_size) == total_size)
    {
      printf("Client sent OK\n");
    }
    else
    {
      printf("Client sent BAD\n");
    }
  }

  m_lastSync = time;
  m_nextSync = time + 5;
  return ret_val;
}
#endif





