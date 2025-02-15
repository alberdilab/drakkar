import sys
import json
from tqdm import tqdm

# Store progress
completed_jobs = 0
total_jobs = None
bar = None

def log_handler(msg):
    global completed_jobs, total_jobs, bar

    # Parse Snakemake log event
    log_event = json.loads(msg)

    # Initialize progress bar when total jobs are known
    if log_event.get("event") == "job_info" and total_jobs is None:
        total_jobs = log_event.get("num_jobs")
        if total_jobs is None:
            return  # Don't start until we know the total number of jobs
        bar = tqdm(total=total_jobs, position=0, leave=True, desc="Workflow Progress")

    # Update progress on job completion
    if log_event.get("event") == "job_finished":
        completed_jobs += 1
        bar.n = completed_jobs
        bar.refresh()

    # Handle workflow completion
    if log_event.get("event") == "workflow_error":
        bar.close()
        print("❌ Workflow failed.")
    elif log_event.get("event") == "workflow_end":
        bar.close()
        print("✅ Workflow completed!")
