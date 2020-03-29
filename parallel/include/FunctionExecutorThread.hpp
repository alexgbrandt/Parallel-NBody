
#ifndef _FUNCTION_EXECUTOR_THREAD_HPP_
#define _FUNCTION_EXECUTOR_THREAD_HPP_

#include "OneWayEventThread.hpp"
#include <exception>

class ExecutorThreadPool;

/**
 * The function executor class implements a long-running
 * thread to execute various functors one after each other.
 */
class FunctionExecutorThread : public OneWayEventThread<std::function<void()>> {

private:

	ExecutorThreadPool* owningPool;

public:

	/**
	 * Construct a stand-alone FunctionExecutorThread.
	 */
	FunctionExecutorThread() :
		OneWayEventThread<std::function<void()>>(),
		owningPool(NULL) {}

	/**
	 * Create a FunctionExecutorThread which belongs
	 * to a the input thread pool.
	 *
	 * @param pool, the owning pool.
	 */
	FunctionExecutorThread(ExecutorThreadPool* pool) :
		OneWayEventThread<std::function<void()>>(),
		owningPool(pool) {}


	/**
	 * Set the owning thread pool of this thread after construction.
	 *
	 * @param pool, the owning pool.
	 */
	void setPool(ExecutorThreadPool* pool) {
		owningPool = pool;
	}

	/**
	 * Destructor.
	 */
	virtual ~FunctionExecutorThread() {}

protected:

	/**
	 * Override for the OneWayEventThread process function.
	 *
	 * @param f, the functor to process.
	 */
	virtual void processRequest(const std::function<void()>& f);
};



#endif
