import logging
from multiprocessing import Process, Queue
from queue import Empty
from random import randint
import threading
import time


def sleepy(id_, item, logq):
    """Example worker function that interprets its input as an integer and
       just sleeps for that amount of time (in seconds).
    """
    logq.put((logging.INFO, f"{id_}: Begin sleep {item}s"))
    time.sleep(item)
    logq.put((logging.INFO, f"{id_}: End sleep {item}s"))


class Bossy:
    def __init__(self, work=None, num_workers=2, worker_function=sleepy, output_log=None):
        self.n = num_workers
        # create queues
        self.work_q = Queue()
        self.log_q = Queue()
        self.worker_fn = worker_function
        self.log = output_log
        # start workers and thread to tail the log messages
        self.processes = self._create_worker_processes()
        self.logging_thread = self._create_logging_thread()
        # add work, including sentinels to stop processes, and wait for it all to be processed
        self.is_done = False
        self._add_work(work)
        sentinels = [None for i in range(num_workers)]
        self._add_work(sentinels)
        self._join_worker_processes()
        # stop logging thread
        self.is_done = True
        self._join_logging_thread()
        # done!

    def _add_work(self, items):
        for item in items:
            self.work_q.put(item)

    def _create_worker_processes(self):
        processes = []
        for i in range(self.n):
            p = Process(target=self.worker, args=(i + 1, self.work_q, self.log_q, self.worker_fn))
            p.start()
            processes.append(p)
        return processes

    @staticmethod
    def worker(id_, q, log_q, func):
        try:
            while True:
                item = q.get()
                if item is None:  # sentinel
                    break
                func(id_, item, log_q)
        except:
            pass

    def _create_logging_thread(self):
        t = threading.Thread(target=self._tail_messages, args=())
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
        print("usage: qtest.py <NUM>")
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
    Bossy(work=delays, num_workers=n, output_log=olog)
    sys.exit(0)
