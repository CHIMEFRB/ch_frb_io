#include <string>
#include <thread>
#include <sstream>
#include <iostream>
#include <mutex>
#include <unistd.h>
#include <zmq.hpp>
#include "chlog.hpp"

using namespace std;
using namespace ch_frb_io;

std::mutex cout_mutex;

static string msg_string(zmq::message_t &msg) {
    return string(static_cast<const char*>(msg.data()), msg.size());
}

static void* server_main(zmq::context_t* ctx, string port, int sid) {
    zmq::socket_t sock(*ctx, ZMQ_SUB);
    sock.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    sock.bind(port);
    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        cout << "server_main " << sid << ": bound " << port << endl;
    }
    for (;;) {
        zmq::message_t msg;
        try {
            if (!sock.recv(&msg)) {
                cout << "log server: failed to receive message" << endl;
                break;
            }
        } catch (const zmq::error_t& e) {
            std::lock_guard<std::mutex> lock(cout_mutex);
            cout << "log server: error receiving message: " << e.what() << endl;
            break;
        }

        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            cout << "Server " << sid << ": " << msg_string(msg) << endl;
        }
    }
    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        cout << "server_main " << sid << " done." << endl;
    }
    return NULL;
}

void log_client(int num) {
    chime_log_set_thread_name(stringprintf("Client-%i", num));
    for (int i=0; i<10; i++) {
        chlog("Hello I am client " << num << ", message " << i);
    }
}

int main() {
    zmq::context_t ctx;

    // Try opening the socket with a given context.
    chime_log_set_name("myputer");
    chime_log_set_thread_name("Main");
    chime_log_open_socket(&ctx);

    chime_log_local(false);

    chime_log_server ser(cout, NULL, "127.0.0.1");
    cout << "Server addr: " << ser.get_address() << endl;
    ser.start();
    chime_log_add_server(ser.get_address());

    //chime_log_server ser2;
    //cout << "Server addr: " << ser2.get_address() << endl;

    /*
    string port = "tcp://127.0.0.1:6666";
    thread serverthread(std::bind(server_main, &ctx, port, 1));
    serverthread.detach();
    chime_log_add_server(port);

    string port2 = "tcp://127.0.0.1:6667";
    thread serverthread2(std::bind(server_main, &ctx, port2, 2));
    serverthread2.detach();
    chime_log_add_server(port2);
     */

    // What if we connect to a server that doesn't exist?
    string port3 = "tcp://127.0.0.1:6668";
    chime_log_add_server(port3);

    usleep(100000);

    chlog("Hello world");
    chlog("Hello " << 1 << ", " << 2 << ", 3");
    chdebug("Debug" << 42+43);

    chlogf("Gotta love %s style; %03i", "printf", 7);

    usleep(100000);

    /*
    chime_log_close_socket();

    // Now re-open with a new context
    chime_log_open_socket(NULL);

    chime_log_add_server(port);
    chime_log_add_server(port2);
    // (wait for connect)
    usleep(100000);
     */

    thread logger1(std::bind(log_client, 1));
    thread logger2(std::bind(log_client, 2));

    logger1.join();
    logger2.join();

    usleep(100000);

    chime_log_close_socket();

    chime_log_local(true);

    chassert(4 < 7);

    // Don't assert(false) when you want your unit tests to succeed!
    //chassert(7 < 4);

    cout << "main() finished" << endl;
    return 0;
}
