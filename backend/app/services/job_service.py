import uuid
from typing import Dict, Any, Optional

class JobService:
    """
    Manages the lifecycle and data (status, context, results) for background jobs.
    """
    def __init__(self):
        # In-memory storage for job status and context.
        # Consider replacing with a persistent store (e.g., Redis, DB) for production.
        self._jobs_status: Dict[str, Dict[str, Any]] = {}
        self._job_contexts: Dict[str, Dict[str, Any]] = {}

    def create_job(self, initial_context: Optional[Dict[str, Any]] = None) -> str:
        """
        Creates a new job entry, initializes its status, and stores initial context.

        Args:
            initial_context: Optional dictionary containing initial data for the job context.

        Returns:
            The unique ID generated for the job.
        """
        job_id = str(uuid.uuid4())
        self._jobs_status[job_id] = {
            "status": "pending",
            "progress": 0.0,
            "message": "Job created and waiting to start.",
            "results": None,
        }
        self._job_contexts[job_id] = initial_context if initial_context else {}
        print(f"JobService: Created job {job_id}")
        return job_id

    def job_exists(self, job_id: str) -> bool:
        """Checks if a job ID exists in the status tracker."""
        return job_id in self._jobs_status

    def update_job_status(
        self,
        job_id: str,
        status: Optional[str] = None,
        progress: Optional[float] = None,
        message: Optional[str] = None,
        results: Optional[Dict[str, Any]] = None,
    ):
        """
        Updates the status, progress, message, and results for a given job.

        Args:
            job_id: The ID of the job to update.
            status: The new status string (e.g., 'running', 'completed', 'failed').
            progress: The completion progress (0.0 to 1.0).
            message: A descriptive message about the current status.
            results: A dictionary containing the final results of the job (if completed).
        """
        if not self.job_exists(job_id):
            print(f"JobService: Attempted to update non-existent job {job_id}")
            # Optionally raise an error or handle gracefully
            return

        if status is not None:
            self._jobs_status[job_id]["status"] = status
        if progress is not None:
            # Clamp progress between 0.0 and 1.0
            self._jobs_status[job_id]["progress"] = max(0.0, min(1.0, progress))
        if message is not None:
            self._jobs_status[job_id]["message"] = message
        if results is not None:
            self._jobs_status[job_id]["results"] = results # Store results directly in status dict

        print(f"JobService: Updated status for job {job_id}: {self._jobs_status[job_id]}")


    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """
        Retrieves the full status dictionary for a given job.

        Args:
            job_id: The ID of the job.

        Returns:
            The status dictionary, or None if the job doesn't exist.
        """
        return self._jobs_status.get(job_id)

    def store_job_context(self, job_id: str, context_data: Dict[str, Any]):
        """
        Stores or updates context data associated with a job.

        Args:
            job_id: The ID of the job.
            context_data: The dictionary containing context data to store.
        """
        if not self.job_exists(job_id):
             print(f"JobService: Attempted to store context for non-existent job {job_id}")
             # Optionally raise an error
             return
        # Merge new data with existing context, overwriting keys if they exist
        self._job_contexts[job_id].update(context_data)
        print(f"JobService: Stored context for job {job_id}")


    def get_job_context(self, job_id: str) -> Optional[Dict[str, Any]]:
        """
        Retrieves the context data associated with a job.

        Args:
            job_id: The ID of the job.

        Returns:
            The context dictionary, or None if the job doesn't exist.
        """
        if not self.job_exists(job_id):
            print(f"JobService: Attempted to get context for non-existent job {job_id}")
            return None
        return self._job_contexts.get(job_id)

# Singleton instance (or manage via FastAPI dependency injection)
# For simplicity here, we use a global instance. FastAPI DI is preferred.
job_service_instance = JobService()

def get_job_service():
    """Dependency function to get the singleton JobService instance."""
    return job_service_instance 