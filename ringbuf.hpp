#ifndef _CH_FRB_IO_RINGBUF_HPP
#define _CH_FRB_IO_RINGBUF_HPP

/**

 Templated "ring buffer".

 This class is meant to hold a data stream.  As a way of controlling
 the amount of memory used, it is allowed to have a maximum number of
 objects alive at any time.  When clients retrieve objects from the
 ring buffer, the objects are kept alive (via shared_ptrs), and
 continue to count toward the ring buffer's quota of live objects.
 Only when the last shared_ptr is dropped does the ring buffer get to
 add a new object.

 The use case is that we need to support callbacks (RPCs) that
 retrieve a lot of data and write it to disk, or send it over the
 network, which could take a long time.  While this is happening, we
 don't want to use more and more memory in the ring buffer; we want to
 control the total memory footprint.  We do this by assuming that the
 objects in the ring buffer are all the same size, so then we just
 need to limit the number of objects that are "live".  The RPC client
 will keep objects live via shared pointers.  When the ring buffer is
 asked to accept more data (via push()) and it is already at its
 limit, then it will begin dropping the oldest data from the buffer.
 If no client is holding a pointer to these objects, then they will be
 freed and space will become available for the new data.  If the whole
 buffer is dropped and space still does not become available, then
 space cannot be allocated for the new object and the push() fails.

 This is implemented via shared_ptr Deleter functionality: the Deleter
 is just a function that is called immediately before an object is
 about to be deleted.  The ring buffer decrements its count of the
 number of live objects when this happens; it then knows it is allowed
 to allocate one more new object.

 The CHIME FRB L1 code subclasses this class; see l1-ringbuf.hpp.

 */


#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <deque>
#include <memory>
#include <vector>

// forward decl
template<class T>
class RingbufDeleter;

template<class T>
class Ringbuf {
    friend class RingbufDeleter<T>;

public:
    /*
     Creates a new Ringbuf with the given maximum quota size.
     */
    Ringbuf(int maxsize);

    virtual ~Ringbuf();

    // Noncopyable
    Ringbuf(const Ringbuf<T> &) = delete;
    Ringbuf<T> &operator=(const Ringbuf<T> &) = delete;

    /*
     Try to push a new frame onto the ring buffer.  If the quota of
     frames has been reached, frames will be popped off the front
     until the buffer is empty.  If space is still not available, will
     return an empty shared_ptr.

     Returns the enqueued shared_ptr if successful.
     */
    std::shared_ptr<T> push(T* t);

    /*
     Check whether there is room to push a new frame onto the ring
     buffer.  This method will try to free up space by dropping frames
     off the front of the queue.  Returns false if no space can be found.
     */
    bool can_push();

    /*
     Pops the oldest frame from this buffer and returns it.  Returns
     an empty shared_ptr if the buffer is empty.
     */
    std::shared_ptr<T> pop();

    /*
     Returns the number of frames queued.
     */
    int size();

    /*
     Returns the maximum number of frames that can be alive.
     */
    int maxsize();

    /*
     Retrieves frames from this ring buffer that satisfy the
     (optional) selection function.
     */
    std::vector<std::shared_ptr<T> > snapshot(bool (*testFunc)(const std::shared_ptr<T>) = NULL);

    /*
     Retrieves frames from this ring buffer (into the given vector)
     that satisfy the (optional) selection function.
     */
    void snapshot(std::vector<std::shared_ptr<T> > &v,
                  bool (*testFunc)(const std::shared_ptr<T>) = NULL);

    void snapshot(std::vector<std::shared_ptr<T> > &v,
                  std::function<bool(const std::shared_ptr<T>)>);

    void visit(std::function<void(std::shared_ptr<T>)> func);

protected:
    // Small helper class that calls deleted() when one of our frames
    // is deleted.
    RingbufDeleter<T> _deleter;

    // The actual queue of frames.
    std::deque<std::shared_ptr<T> > _q;

    // (and the mutex for it)
    pthread_mutex_t _q_lock;

    // Number of live frames.
    size_t _live;

    // Maximum allowable number of live frames.
    size_t _maxsize;

    // (and the mutex for it)
    pthread_mutex_t _live_lock;

    // We're about to drop this frame.
    virtual void dropping(std::shared_ptr<T> t);

    // Called by the RingbufDeleter when a shared_ptr is deleted
    void deleted(T* t);

};

// Helper class that is the Deleter for our shared_ptr<frames>.  Calls
// Ringbuf.deleted() to track number of live frames.
template<class T>
class RingbufDeleter {
public:
    RingbufDeleter(Ringbuf<T>* rb) : _ringbuf(rb) {};
    void operator()(T* t) {
        _ringbuf->deleted(t);
    }

protected:
    Ringbuf<T>* _ringbuf;
};

// Implementation... must appear here because Ringbuf is templated!
#include "ringbuf-impl.hpp"


#endif
