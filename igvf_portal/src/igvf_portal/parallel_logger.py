import dataclasses
import logging
from contextlib import nullcontext
from multiprocessing.synchronize import Lock as ProcessLock
from threading import Lock as ThreadLock


@dataclasses.dataclass(slots=True)
class ParallelLogger:
    logger: logging.Logger
    lock: ThreadLock | ProcessLock | nullcontext

    @classmethod
    def new(
        cls,
        logger: logging.Logger,
        lock: ThreadLock | ProcessLock | nullcontext | None = None,
    ) -> ParallelLogger:
        return ParallelLogger(
            logger=logger, lock=nullcontext() if lock is None else lock
        )

    def debug(self, message: str) -> None:
        with self.lock:
            self.logger.debug(message)

    def info(self, message: str) -> None:
        with self.lock:
            self.logger.info(message)

    def warning(self, message: str) -> None:
        with self.lock:
            self.logger.warning(message)

    def error(self, message: str) -> None:
        with self.lock:
            self.logger.error(message)

    def critical(self, message: str) -> None:
        with self.lock:
            self.logger.critical(message)

    def fatal(self, message: str) -> None:
        with self.lock:
            self.logger.fatal(message)
