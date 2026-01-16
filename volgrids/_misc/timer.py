import time

# //////////////////////////////////////////////////////////////////////////////
class Timer:
    def __init__(self, text = '', end = ' '):
        if text: print(text, end = end, flush = True)
        self.t = None

    # --------------------------------------------------------------------------
    def start(self) -> "Timer":
        self.t = time.time()
        return self

    # --------------------------------------------------------------------------
    def end(self, text = '', end = '\n', minus: float = 0.0) -> float:
        if self.t is None:
            raise RuntimeError("Timer not started. Call start() before end().")
        elapsed = (time.time() - self.t) - minus
        preffix = f"{text}: " if text else ''
        print(
            f"({preffix}{int(elapsed // 60)}m {elapsed % 60:.2f}s)",
            end = end, flush = True
        )
        return elapsed


# //////////////////////////////////////////////////////////////////////////////
