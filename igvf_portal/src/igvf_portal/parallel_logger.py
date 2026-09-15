import dataclasses
import logging
from contextlib import nullcontext
from multiprocessing.synchronize import Lock as ProcessLock
from threading import Lock as ThreadLock


@dataclasses.dataclass(slots=True)
class ParallelLogger:
    """Logger for parallel environments that keeps threads/processes from clobbering each other."""

    logger: logging.Logger
    lock: ThreadLock | ProcessLock | nullcontext

    @classmethod
    def new(
        cls,
        logger: logging.Logger,
        lock: ThreadLock | ProcessLock | nullcontext | None = None,
    ) -> ParallelLogger:
        """Create a ParallelLogger, using a no-op lock if none is supplied.

        Args:
            logger: Underlying logger that messages are forwarded to.
            lock: Lock shared by parallel workers. If None, a no-op lock is used.

        Returns:
            New ParallelLogger wrapping `logger`.
        """
        return ParallelLogger(logger=logger, lock=nullcontext() if lock is None else lock)

    def setLevel(self, level: int | str) -> None:  # noqa: N802 - mirrors logging.Logger API
        """Set the level of the underlying logger."""
        self.logger.setLevel(level)

    def getEffectiveLevel(self) -> int:  # noqa: N802 - mirrors logging.Logger API
        """Get the effective level of the underlying logger."""
        return self.logger.getEffectiveLevel()

    def debug(self, message: str) -> None:
        """Log `message` at DEBUG level while holding the lock."""
        with self.lock:
            self.logger.debug(message)

    def info(self, message: str) -> None:
        """Log `message` at INFO level while holding the lock."""
        with self.lock:
            self.logger.info(message)

    def warning(self, message: str) -> None:
        """Log `message` at WARNING level while holding the lock."""
        with self.lock:
            self.logger.warning(message)

    def error(self, message: str) -> None:
        """Log `message` at ERROR level while holding the lock."""
        with self.lock:
            self.logger.error(message)

    def critical(self, message: str) -> None:
        """Log `message` at CRITICAL level while holding the lock."""
        with self.lock:
            self.logger.critical(message)

    def fatal(self, message: str) -> None:
        """Log `message` at FATAL level while holding the lock."""
        with self.lock:
            self.logger.fatal(message)
