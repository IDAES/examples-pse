"""
Utility class to run things in parallel.
"""
import logging
from multiprocessing import Process, Queue
from queue import Empty
from random import randint
import threading
import time
from typing import List


def sleepy(id_, item, logq, **kwargs):
    """Example worker function that interprets its input as an integer and
       just sleeps for that amount of time (in seconds).
    """
    logq.put((logging.INFO, f"{id_}: Begin sleep {item}s"))
    time.sleep(item)
    logq.put((logging.INFO, f"{id_}: End sleep {item}s"))
    return "z" * item


class Bossy:
    def __init__(self, work=None, num_workers=2, worker_function=sleepy, output_log=None,
                 **kwargs):
        self.n = num_workers
        # create queues
        self.work_q = Queue()
        self.log_q = Queue()
        self.result_q = Queue()
        self.worker_fn = worker_function
        self.log = output_log
        self.is_done = False
        self._worker_kwargs = kwargs
        self.processes, self.logging_thread = None, None
        # add work, including sentinels to stop processes, and wait for it all to be processed
        self._add_work(work)
        sentinels = [None for i in range(self.n)]
        self._add_work(sentinels)

    def run(self) -> List:
        results = [None] * self.n
        if self.processes or self.logging_thread:
            return results
        # start workers and thread to tail the log messages
        self.processes = self._create_worker_processes(self._worker_kwargs)
        self.logging_thread = self._create_logging_thread()
        self._join_worker_processes()
        # stop logging thread
        self.is_done = True
        self._join_logging_thread()
        self.processes, self.logging_thread = None, None
        # return results as a vector
        for i in range(self.n):
            id_, result = self.result_q.get_nowait()
            index = id_ - 1
            results[index] = result
        return results

    def _add_work(self, items):
        for item in items:
            self.work_q.put(item)

    def _create_worker_processes(self, kwargs):
        processes = []
        for i in range(self.n):
            p = Process(target=self.worker, args=(i + 1, self.work_q, self.log_q, self.result_q, self.worker_fn,
                                                  kwargs))
            p.start()
            processes.append(p)
        return processes

    @staticmethod
    def worker(id_, q, log_q, result_q, func, kwargs):
        result = None
        try:
            while True:
                item = q.get()
                if item is None:  # sentinel
                    break
                result = func(id_, item, log_q, **kwargs)
        except:
            pass
        result_q.put((id_, result))

    def _create_logging_thread(self):
        t = threading.Thread(target=self._tail_messages, args=(), daemon=True)
        t.start()
        return t

    def _tail_messages(self):
        while True:
            try:
                level, msg = self.log_q.get(True, 2)
                self.log.log(level, msg)
            except Empty:
                if self.is_done:
                    return

    def _join_worker_processes(self):
        for p in self.processes:
            p.join()

    def _join_logging_thread(self):
        self.logging_thread.join()


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("usage: bossy.py <NUM>")
        sys.exit(1)
    try:
        n = int(sys.argv[1])
        if n < 1:
            raise ValueError()
    except ValueError:
        print(f"{sys.argv[1]} is not a positive integer")
        print("usage: qtest.py <NUM>")
        sys.exit(1)

    delays = [randint(3, 8) for i in range(n)]
    olog = logging.getLogger()
    h = logging.StreamHandler()
    olog.addHandler(h)
    olog.setLevel(logging.INFO)
    b = Bossy(work=delays, num_workers=n, output_log=olog)
    r = b.run()
    print(f"Results: {r}")

    sys.exit(0)
