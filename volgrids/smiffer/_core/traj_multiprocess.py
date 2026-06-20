import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class TrajMultiprocess:
    def __init__(self, app: "sm.AppSmiffer"):
        ### Set by the parent before spawning the trajectory MP pool; inherited by workers via fork.
        self._worker_app: "sm.AppSmiffer" = app


    # --------------------------------------------------------------------------
    def run(self, n_frames: int) -> None:
        import multiprocessing as mp
        from collections import deque

        ### "fork" so children inherit module state (configs, _WORKER_APP, etc.) without pickling
        ctx = mp.get_context("fork")
        vg.MP_CMAP_LOCK = ctx.Lock()

        nproc = min(self._worker_app.nproc, n_frames)
        try:
            with ctx.Pool(nproc, initializer = self._worker_init) as pool:
                deque(pool.imap_unordered(self._worker_process_frame, range(n_frames)), maxlen = 0)
        finally:
            vg.MP_CMAP_LOCK = None


    # ------------------------------------------------------------------------------
    def _worker_init(self):
        """Re-create the MolSystem in each worker so file handles aren't shared across forks."""
        self._worker_app.ms = sm.MolSystem(sm.PATH_STRUCT, self._worker_app.path_traj)


    # ------------------------------------------------------------------------------
    def _worker_process_frame(self, frame_idx: int) -> None:
        return self._worker_app.process_frame(frame_idx)



# //////////////////////////////////////////////////////////////////////////////
