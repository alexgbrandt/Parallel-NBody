
#ifndef _ONEWAY_EVENT_THREAD_HPP_
#define _ONEWAY_EVENT_THREAD_HPP_

#if defined(SERIAL) && SERIAL
#define ONEWAY_EVENT_THREAD_PARALLEL 0
#else
#define ONEWAY_EVENT_THREAD_PARALLEL 1
#endif

#if ONEWAY_EVENT_THREAD_PARALLEL
#include "Synchronized.hpp"
#include "AsyncObjectStream.hpp"
#include <thread>
#include <condition_variable>
#endif

#include <vector>
#include <map>
#include <iostream>


/**
 * The OneWayEventThread class represents a long-running, event-loop based
 * thread for unified access to some sub-module. It supports
 * many requesters and exaclty one responder--the event thread.
 *
 * The design is such that it automatically sets up a
 * processing thread on creation.
 *
 * This class is abstract. One should subclass it and specialize
 * the processTask method to actually perform the work required
 * to turn the Request object into a Response object.
 * If no response is required, one simply returns false.
 * See the documentation of the processTask method.
 *
 * The class can easily be used with the intent that
 * the processing thread will be a single point of
 * access to some sub-system, simply make a singleton
 * of the sub-class
 *
 * The class is templated by the object that it should
 * receive as a request and return as a response.
 *
 */
template<class Request>
class OneWayEventThread {

protected:

#if ONEWAY_EVENT_THREAD_PARALLEL
	AsyncObjectStream<Request> requestQueue;

	std::thread m_worker;

#endif

#if ONEWAY_EVENT_THREAD_PARALLEL
	virtual void eventLoop() {
		Request reqObj;
		while (requestQueue.getNextObject(reqObj)) {
			// std::cerr << std::this_thread::get_id() << " is about to process request" << std::endl;
			processRequest(reqObj);
		}

	}

#endif

protected:

	/**
	 * Clean up the inter-thread communication resources
	 * and any resources used by the thread.
	 *
	 * Dervived classes should call this super-class method
	 * after they have cleaned up resources specific to their
	 * dervied implementation.
	 */
	virtual void threadCleanup() {
#if ONEWAY_EVENT_THREAD_PARALLEL
		requestQueue.resultsFinished();
		if (m_worker.joinable()) {
			m_worker.join();
		}
#endif
	}

	virtual void processRequest(const Request& reqObj) = 0;

public:

	/**
	 * Create a new event thread. Spins up the
	 * worker thread and the collector thread.
	 */
	OneWayEventThread()
#if ONEWAY_EVENT_THREAD_PARALLEL
	:	requestQueue(),
		m_worker()
	{
		m_worker = std::thread(&OneWayEventThread::eventLoop, this);
	}
#else
	{}
#endif


	/**
	 * Destructor.
	 */
	virtual ~OneWayEventThread() {
		threadCleanup();
	}


	/* Event threads are not copyable */
	OneWayEventThread(const OneWayEventThread& other) = delete;

	OneWayEventThread& operator=(const OneWayEventThread&) = delete;


	/**
	 * Implements an asynchronous request send.
	 *
	 * @param reqObj: the Request to send.
	 */
	virtual void sendRequest(const Request& reqObj) {
#if ONEWAY_EVENT_THREAD_PARALLEL
		Request sendCopy = reqObj;
		requestQueue.addResult(sendCopy);
#else
		processRequest(reqObj);
#endif
	}


};

#endif
