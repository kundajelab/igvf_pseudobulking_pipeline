import dataclasses
import json
import logging
import multiprocessing
from collections.abc import Collection
from functools import cached_property
from multiprocessing.synchronize import Lock as ProcessLock
from threading import Lock as ThreadLock
from typing import Callable

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.enums import Concurrency, IgvfMode


@dataclasses.dataclass(slots=False, kw_only=True)
class RegisterConfig:
    igvf_mode: IgvfMode
    dry_run: bool
    profile_id: str
    num_tries: int
    delay: float
    backoff: float
    overwrite_array_values: bool
    remove_properties: list[str]
    upload_file: bool
    upload_duplicate: bool
    continue_on_failed_credentials: bool = True
    concurrency: Concurrency = Concurrency.NONE

    @cached_property
    def thread_lock(self) -> ThreadLock:
        return ThreadLock()

    @cached_property
    def process_lock(self) -> ProcessLock:
        return multiprocessing.Lock()

    @property
    def rm_patch(self) -> bool:
        return len(self.remove_properties) > 0

    @property
    def cleaned_profile_id(self) -> str:
        profile_id = self.profile_id.strip("/").split("/", 1)[0].lower()
        # Multi-word profile names are hypen-separated, i.e. genetic-modifications.
        profile_id = profile_id.replace("-", "_")
        return profile_id

    @property
    def new_connection(self) -> PConnection:
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
        )

    def retry[**P, R](
        self,
        no_retry_exceptions: Collection[type] = (json.decoder.JSONDecodeError,),
        logger: logging.Logger | None = None,
    ) -> Callable[[Callable[P, R]], Callable[P, R]]:
        return utils.retry(
            num_tries=self.num_tries,
            delay=self.delay,
            backoff=self.backoff,
            no_retry_exceptions=no_retry_exceptions,
            logger=logger,
        )
