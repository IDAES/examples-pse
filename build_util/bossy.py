"""
Utility class to run things in parallel.
"""
import logging
from multiprocessing import Process, Queue
from queue import Empty
from random import randint
import signal
import sys
import threading
import time
from typing import List


class WorkerInterrupted(Exception):
    def __init__(self, why):
        super().__init__(f"Worker interrupted: {why}")


def sleepy(id_, item, logq, **kwargs):
    """Example worker function that interprets its input as an integer and
       just sleeps for that amount of time (in seconds).
    """
    logq.put((logging.INFO, f"{id_}: Begin sleep {item}s"))
    time.sleep(item)
    logq.put((logging.INFO, f"{id_}: End sleep {item}s"))
    return "z" * item


class Bossy:
    def __init__(
        self,
        work=None,
        num_workers=2,
        worker_function=sleepy,
        output_log=None,
        **kwargs,
    ):
        self.n = num_workers
        # create queues
        self.work_q = Queue()
        self.log_q = Queue()
        self.done_q = Queue()
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
        self._handle_signals()

    def _handle_signals(self):
        if hasattr(signal, 'SIGINT'):
            signal.signal(signal.SIGINT, self.signal_handler)
        if hasattr(signal, 'SIGBREAK'):
            signal.signal(signal.SIGBREAK, self.signal_handler)

    def signal_handler(self, sig, frame):
        """User interrupted the program with a Control-C or by sending it a signal (UNIX).
        """
        self.log.warning("Interrupted. Killing child processes.")
        if self.processes:
            self.log.warning(f"Killing {len(self.processes)} processes")
            for p in self.processes:
                p.kill()

    def run(self) -> List:
        self.log.info("Run workers: begin")
        results = []
        if self.processes or self.logging_thread:
            return results
        # start workers and thread to tail the log messages
        self.processes = self._create_worker_processes(self._worker_kwargs)
        self.logging_thread = self._create_logging_thread()
        self._join_worker_processes()
        self.processes = None
        # stop logging thread
        self.is_done = True
        self._join_logging_thread()
        self.logging_thread = None
        # return results as a list of tuples (worker, result)
        self.log.debug(f"Collect {self.n} results: begin")
        while not self.result_q.empty():
            id_, result_list = self.result_q.get_nowait()
            self.log.debug(f"Worker [{id_}]: got result")
            for r in result_list:
                results.append((id_, r))
            self.log.debug(f"Worker [{id_}]: Recorded result")
        self.log.debug(f"Collect {self.n} results: end")
        self.log.info("Run workers: end")
        return results

    def _add_work(self, items):
        for item in items:
            self.work_q.put(item)

    def _create_worker_processes(self, kwargs):
        self.log.debug(f"Create worker processes. kwargs={kwargs}")
        processes = []
        g_log = self.log
        for i in range(self.n):
            p = Process(
                target=self.worker,
                args=(
                    i,
                    self.work_q,
                    self.log_q,
                    self.result_q,
                    self.done_q,
                    self.worker_fn,
                    {},
                ),
            )  # kwargs
            self.log.debug(f"Worker [{i + 1}]: Starting process {p}")
            p.start()
            processes.append(p)
        return processes

    @staticmethod
    def worker(id_, q, log_q, result_q, done_q, func, kwargs):
        pfx = f"[Worker {id_} Main-Loop]"

        def log_info(m, pfx=pfx):
            log_q.put((logging.INFO, f"{pfx}: {m}"))

        def log_debug(m, pfx=pfx):
            log_q.put((logging.DEBUG, f"{pfx}: {m}"))

        def log_error(m, pfx=pfx):
            log_q.put((logging.ERROR, f"{pfx}: {m}"))

        log_info("begin")
        result_list = []
        while True:
            item = q.get()
            log_debug("Got next item of work from queue")
            if item is None:  # sentinel
                log_info("No more work: Stop.")
                break
            try:
                log_debug(f"Run worker function: begin")
                result = func(id_, item, log_q, **kwargs)
                log_debug(f"Run worker function: end")
                result_list.append(result)
            except KeyboardInterrupt:
                log_debug(f"Run worker function: end (keyboard interrupt)")
                raise WorkerInterrupted("Keyboard interrupt")
            except Exception as err:
                log_error(f"Run worker function: end (exception): {err}")
        # Put results on the queue
        if result_list:
            result_q.put((id_, result_list))
        done_q.put(id_)

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
        """Unfortunately, simply joining processes doesn't seem to do the trick all the time.
        Instead, we use a special queue to keep track of which workers have finished, then,
        if they have not stopped on their own, forcibly terminate the associated processes
        in order to join() them. If even after that the join fails, we mark the process as
        'unjoinable' and give up (!)
        """
        num_joined, num_unjoinable, num_proc = 0, 0, len(self.processes)
        while num_joined + num_unjoinable < num_proc:
            self.log.debug(f"Waiting for {num_proc - num_joined - num_unjoinable} processes to finish")
            try:
                id_ = self.done_q.get(timeout=60)
            except Empty:
                # Interruptible wait allowing for control-c
                time.sleep(1)
                continue
            proc = self.processes[id_]
            if proc.is_alive():
                self.log.info(f"Terminating process: {proc}")
                proc.terminate()
                t0 = time.time()
                while proc.is_alive():
                    time.sleep(1)
                    if time.time() - t0 > 10:
                        break
            if proc.is_alive():
                self.log.error(f"Could not terminate process: {proc}")
                num_unjoinable += 1
            else:
                self.log.debug(f"Joining process: {proc}")
                proc.join()
                num_joined += 1
        if num_unjoinable > 0:
            self.log.error(f"{num_unjoinable} processes could not be joined")

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
