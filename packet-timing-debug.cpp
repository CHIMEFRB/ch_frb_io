#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <string>
#include <unordered_map>
#include "ch_frb_io_internals.hpp"
#include "ch_frb_io.hpp"
#include "chlog.hpp"

using namespace std;
using namespace ch_frb_io;

int main() {
  int sockfd;
  struct timeval tv_ini = xgettimeofday();

  sockfd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
  if (sockfd < 0)
    throw runtime_error(string("ch_frb_io: socket() failed: ") + strerror(errno));

  int flags = fcntl(sockfd, F_GETFD);
  flags |= FD_CLOEXEC;
  if (fcntl(sockfd, F_SETFD, flags) < 0)
    throw runtime_error(string("ch_frb_io: couldn't set close-on-exec flag on socket file descriptor") + strerror(errno));

  int socket_bufsize = 128 * 1024 * 1024;
  // bufsize
  int err = setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, (void *)&socket_bufsize, sizeof(socket_bufsize));
  if (err < 0)
    throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVBUF) failed: ") + strerror(errno));

  // timeout
  int socket_timeout_usec = 10000;                // 0.01 sec
  const struct timeval tv_timeout = { 0, socket_timeout_usec };
  err = setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, &tv_timeout, sizeof(tv_timeout));
  if (err < 0)
    throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVTIMEO) failed: ") + strerror(errno));

  //int udp_port = 6677;
  int udp_port = 1313;

  // We used to just hard-code the ipaddr....
  // string ipaddr = "10.6.201.11";  // cf1n1
  // err = inet_pton(AF_INET, ipaddr.c_str(), &server_address.sin_addr);
  // if (err <= 0)
  //   throw runtime_error(ipaddr + ": inet_pton() failed (note that no DNS lookup is done, the argument must be a numerical IP address)");

  struct sockaddr_in server_address;
  memset(&server_address, 0, sizeof(server_address));

  // instead, now do this overly-elaborate dance...
  char hostname[1024];
  hostname[1023] = '\0';
  if (gethostname(hostname, 1023) < 0)
    throw runtime_error(string("ch_frb_io: gethostname failed: ") + strerror(errno));
  chlog("hostname: " << string(hostname));
  struct addrinfo hints, *info, *p;
  memset(&hints, 0, sizeof hints);
  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_flags = AI_CANONNAME;
  if (getaddrinfo(hostname, "http", &hints, &info) != 0)
      throw runtime_error(string("getaddrinfo failed: ") + strerror(errno));
  for (p = info; p != NULL; p = p->ai_next) {
      //chlog("hostname: " << p->ai_canonname);
      if (p->ai_family == AF_INET) {
          memcpy(&server_address, p->ai_addr, sizeof(server_address));
          break;
      }
  }
  freeaddrinfo(info);
  string ipaddr = string(inet_ntoa(server_address.sin_addr));

  server_address.sin_family = AF_INET;
  server_address.sin_port = htons(udp_port);

  err = ::bind(sockfd, (struct sockaddr *) &server_address, sizeof(server_address));
  if (err < 0)
    throw runtime_error(string("ch_frb_io: bind() failed (" + ipaddr + ":" + to_string(udp_port) + "): " + strerror(errno)));

  chlog("Bound port: " << ipaddr << ":" << udp_port);

  uint64_t first_fpga = 0;
  // Wait 5 seconds...
  uint64_t wait_fpga = 384 * 1024 * 5;

  double t0 = 0;

  unordered_map<string, double> first_from;

  for (;;) {
    char packet_data[10000];
    // Record the sender IP & port here
    sockaddr_in sender_addr;
    int slen = sizeof(sender_addr);
    int packet_nbytes = ::recvfrom(sockfd, packet_data, sizeof(packet_data), 0,
				   (struct sockaddr *)&sender_addr, (socklen_t *)&slen);
    //chlog("Received " << packet_nbytes);
    struct timeval curr_tv;
    curr_tv = xgettimeofday();
    double curr_timestamp = usec_between(tv_ini, curr_tv);

    // Check for error or timeout in read()
    if (packet_nbytes < 0) {
      if ((errno == EAGAIN) || (errno == ETIMEDOUT))
	continue;  // normal timeout
      if (errno == EINTR)
	continue; // this can happen when running in gdb
      throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
    }

    intensity_packet packet;
    if (!packet.decode(reinterpret_cast<uint8_t*>(packet_data), packet_nbytes)) {
      chlog("Bad packet!");
      continue;
    }

    //chlog("Packet from " << ip_to_string(sender_addr) << " at " << curr_timestamp/1e6 << " sec"
    //<< ": FPGA " << packet.fpga_count);

    if (!first_fpga) {
      first_fpga = packet.fpga_count;
      chlog("Received packet with FPGA count " << first_fpga << "; waiting " << wait_fpga << " counts");
      continue;
    }
    if (packet.fpga_count < first_fpga + wait_fpga) {
      // waiting...
      continue;
    }
    if (packet.fpga_count == first_fpga + wait_fpga) {
      string sender = ip_to_string(sender_addr);
      if (t0 == 0.)
	t0 = curr_timestamp/1e6;
      auto p = first_from.find(sender);
      if (p == first_from.end()) {
	chlog("First packet from " << sender << " at " << curr_timestamp/1e6-t0);
	first_from[sender] = curr_timestamp/1e6 - t0;
      }
      continue;
    }
    if (packet.fpga_count > first_fpga + 2 * wait_fpga) {
      // report & reset!
      //for (auto p: first_from) {
      //chlog("First packet from " << p.first << " at " << p.second);
      //}
      first_fpga = 0;
      first_from.clear();
      t0 = 0.;
      //continue;
      break;
    }
    

  }


  return 0;
}

