
#ifndef _EXECUTOR_THREAD_POOL_HPP_
#define _EXECUTOR_THREAD_POOL_HPP_


class FunctionExecutorThread;

#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <thread>

/**
 * A class implementing a thread pool and FunctionExecutorThread,
 * i.e. long-running threads executing functors asynchronously.
 */
class ExecutorThreadPool {

private:

	size_t nThreads;
	FunctionExecutorThread* poolThreads;

	std::deque<FunctionExecutorThread*> threadPool;
	std::deque<std::function<void()>> taskPool;

	size_t nRetiredThreads;
	size_t nPriorityThreads;
	size_t maxPriorityThreads;
	std::vector<FunctionExecutorThread*> priorityThreads;

	std::mutex m_mutex;
	std::condition_variable m_cv;

	void putbackThread(FunctionExecutorThread* t);

	void tryPullTask();

	friend FunctionExecutorThread;

	/**
	 * Create a pool of FunctionExecutorThreads of size n.
	 * n: size of thread pool.
	 */
	ExecutorThreadPool(int n = std::thread::hardware_concurrency());

	~ExecutorThreadPool();

public:

	/**
	 * Obtain a reference to the singleton ExecutorThreadPool object.
	 * @return a reference to the singleton ExecutorThreadPool.
	 */
	static ExecutorThreadPool& getThreadPool();

	/**
	 * The number of threads in the thread pool.
	 */
	static int maxThreads;

	/**
	 * Add a task to the executor.
	 * @param f, the function task to execute.
	 */
	void addTask(std::function<void()>& f);

	/**
	 * Add a task to the executor for a thread at a particular index.
	 * This is an advanced option to allow for data locality within a thread.
	 * If the thread of index idx is already busy, this method does nothing.
	 *
	 * @param f, the function task to execute.
	 * @param idx, the index of the thread to execute this task.
	 */
	void addTaskAtIdx(std::function<void()>& f, int idx);


	/**
	 * Add a task to th executor pool which should be executed before other tasks.
	 * If the thread pool is empty, a temporary new thread will be created to service this task.
	 * This new priority thread will replace a standard thread once the
	 * standard thread returns to the thread pool.
	 *
	 * @param f, the function task to execute.
	 */
	void addPriorityTask(std::function<void()>& f);


	/**
	 * A query to see if all threads are busy and the pool is empty.
	 *
	 * @return true iff the pool is empty.
	 */
	bool allThreadsBusy();


	/**
	 * A blocking function call to wait for all tasks and threads to finish.
	 */
	void waitForAllThreads();
};


#endif
