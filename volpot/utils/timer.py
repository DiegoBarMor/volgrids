import time

# //////////////////////////////////////////////////////////////////////////////
class Timer:
    def __init__(self, display = ''):
        if display: print(display, end = ' ', flush = True)
        self.start = time.time()

    def end(self):
        elapsed = time.time() - self.start
        print(f"({int(elapsed // 60)}m {elapsed % 60:.2f}s)", flush = True)


# //////////////////////////////////////////////////////////////////////////////
