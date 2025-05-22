import time

# //////////////////////////////////////////////////////////////////////////////
class Timer:
    def __init__(self, display = ''):
        if display: print(display, end = ' ', flush = True)
        self.t = None

    def start(self):
        self.t = time.time()

    def end(self):
        if self.t is None:
            raise RuntimeError("Timer not started. Call start() before end().")
        elapsed = time.time() - self.t
        print(f"({int(elapsed // 60)}m {elapsed % 60:.2f}s)", flush = True)


# //////////////////////////////////////////////////////////////////////////////
