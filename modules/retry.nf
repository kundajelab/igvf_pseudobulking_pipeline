// Shared out-of-memory retry policy for processes whose memory needs cannot be predicted up front.
//
// Such a process estimates a starting requirement as `baseMem` in its script block, then uses these
// two helpers for its memory and maxRetries directives:
//
//     include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'
//
//     process EXAMPLE {
//         memory { oomMemoryOf(task, baseMem) }
//         maxRetries { oomMaxRetriesOf(task, baseMem) }
//         ...
//         script:
//         baseMem = 2.GB + 1.GB * someInputSize
//
// NOTE: these have to be called from inside a directive closure rather than returning a closure of
// their own. A bare `oomMemoryOf()` in a process body is intercepted by nextflow's ProcessConfig,
// which treats any unknown method call as a directive name and fails with "Unknown process
// directive". Keeping the closure at the call site also means `task` and `baseMem` resolve against
// the task context in the usual way, with no delegate juggling.

// Whether Slurm recorded killing processes of the job in workDir for running out of memory, in the
// job's .command.log: "error: Detected N oom_kill events in StepId=...". Only the end of the log is
// read, which is where Slurm writes it.
def slurmReportedOutOfMemory(workDir) {
    if (workDir == null) {
        return false
    }
    def log_file = java.nio.file.Paths.get(workDir.toString(), ".command.log")
    try {
        if (!java.nio.file.Files.exists(log_file)) {
            return false
        }
        // NOTE: a ByteBuffer rather than a byte array, which nextflow's script syntax cannot create
        java.nio.channels.FileChannel.open(log_file).withCloseable { log ->
            def tail_size = Math.min(log.size(), 65536L)
            def tail = java.nio.ByteBuffer.allocate(tail_size as int)
            // a read of a regular file fills the buffer, up to the end of the file
            log.read(tail, log.size() - tail_size)
            return new String(tail.array(), 0, tail.position(), "UTF-8").contains("oom_kill")
        }
    } catch (IOException _exception) {
        return false
    }
}

// Whether the previous attempt looks like it ran out of memory: killed with SIGKILL (137) or a
// signal up to 140, or with SIGBUS (135), which a job can get when it touches memory the kernel
// cannot provide under memory pressure. Or exit status 1 when Slurm logged that it killed processes
// for running out of memory: if that also kills the memory monitor of nextflow's trace wrapper, the
// wrapper fails to stop it and exits 1 itself, hiding the task's 137.
def wasOutOfMemory(task) {
    task.exitStatus == 135 || task.exitStatus in 137..140 || (
        task.exitStatus == 1
        && slurmReportedOutOfMemory(task.previousTrace?.get("workdir") ?: task.workDir)
    )
}

// Memory to request for this attempt: the estimate on the first try, then double it each time the
// previous attempt looked like it ran out of memory. A previous attempt that died for some other
// reason (preemption, say) is retried at the same size it already had.
def oomMemoryOf(task, baseMem) {
    if (task.previousTrace) {
        wasOutOfMemory(task) ? 2 * task.previousTrace.memory : task.previousTrace.memory
    } else {
        baseMem
    }
}

// Retries to allow. Because oomMemoryOf doubles the request, the number of out-of-memory retries so
// far is the log base 2 of how far the current request has grown past baseMem. Once that reaches
// oom_max_retries and we are still running out of memory, stop; the estimate is wrong in a way more
// memory will not fix. Otherwise allow enough retries to also absorb preemptions, which only happen
// on a cluster queue.
def oomMaxRetriesOf(task, baseMem) {
    def wasOom = wasOutOfMemory(task)
    // An attempt that ran out of memory while already pinned to the executor's memory ceiling
    // cannot be helped by trying again: oomMemoryOf would double the request, resourceLimits would
    // clamp it straight back to the same ceiling, and the attempt would be identical. Give up now
    // rather than spend the rest of the budget re-running it. This also keeps previousOomCount
    // meaningful, since a clamped request stops growing and would otherwise stall the check below.
    def memoryCeiling = task.resourceLimits?.memory
    def previousMem = task.previousTrace?.memory
    if (wasOom && memoryCeiling && previousMem && previousMem >= memoryCeiling.toBytes()) {
        return 0
    }
    def previousOomCount = (
        task.previousTrace ?
        (Math.log(previousMem as Float / baseMem.toBytes() as Float) / Math.log(2.0)).round() :
        0
    ) as Integer
    task.executor == 'local' || (
        previousOomCount >= params.oom_max_retries && wasOom
    ) ?
    params.oom_max_retries :
    params.preemptible_max_retries + params.oom_max_retries
}
