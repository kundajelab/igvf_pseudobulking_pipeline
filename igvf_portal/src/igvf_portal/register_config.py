import dataclasses
import multiprocessing
from collections.abc import Sequence
from functools import cached_property
from multiprocessing.synchronize import Lock as ProcessLock
from threading import Lock as ThreadLock
from typing import (
    Final,
)

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.enums import Concurrency, IgvfMode

DRY_RUN: Final[bool] = False
NUM_TRIES: Final[int] = 3
DELAY: Final[float] = 5.0
BACKOFF: Final[float] = 2.0
OVERWRITE_ARRAY_VALUES: Final[bool] = False
REMOVE_PROPERTIES: Final[tuple[str, ...]] = ()
UPLOAD_FILE: Final[bool] = True
UPLOAD_DUPLICATE: Final[bool] = True
CONTINUE_ON_FAILED_CREDENTIALS: Final[bool] = True
EXPECT_PATCH: Final[bool] = False


@dataclasses.dataclass(slots=False, kw_only=True)
class RegisterConfig:
    """Settings for registering (posting or patching) records on the IGVF Portal."""

    igvf_mode: IgvfMode
    profile_id: str
    dry_run: bool = DRY_RUN
    num_tries: int = NUM_TRIES
    delay: float = DELAY
    backoff: float = BACKOFF
    overwrite_array_values: bool = OVERWRITE_ARRAY_VALUES
    remove_properties: Sequence[str] = REMOVE_PROPERTIES
    upload_file: bool = UPLOAD_FILE
    upload_duplicate: bool = UPLOAD_DUPLICATE
    continue_on_failed_credentials: bool = CONTINUE_ON_FAILED_CREDENTIALS
    expect_patch: bool = EXPECT_PATCH

    concurrency: Concurrency = Concurrency.NONE

    @cached_property
    def thread_lock(self) -> ThreadLock:
        """Lock shared by connections when running with thread concurrency."""
        return ThreadLock()

    @cached_property
    def process_lock(self) -> ProcessLock:
        """Lock shared by connections when running with process concurrency."""
        return multiprocessing.Lock()

    @property
    def rm_patch(self) -> bool:
        """Whether any properties should be removed when patching."""
        return len(self.remove_properties) > 0

    @property
    def new_connection(self) -> PConnection:
        """New submission connection, using the lock appropriate to the concurrency mode."""
        match self.concurrency:
            case Concurrency.NONE:
                lock = None
            case Concurrency.THREAD:
                lock = self.thread_lock
            case Concurrency.PROCESS:
                lock = self.process_lock
        return PConnection.new(
            igvf_mode=self.igvf_mode,
            submission=True,
            dry_run=self.dry_run,
            lock=lock,
            continue_on_failed_credentials=self.continue_on_failed_credentials,
            upload_retry=utils.RetryPolicy(
                num_tries=self.num_tries, delay=self.delay, backoff=self.backoff
            ),
        )
