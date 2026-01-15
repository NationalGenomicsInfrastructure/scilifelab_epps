#!/usr/bin/env python
DESC = """Module called by other EPP scripts to write notes to couchdb
"""
import os
import sys
from typing import Any

import requests
import yaml

from scilifelab_epps.utils.genstat_conn import create_jwt_token, email_responsible


def write_note_to_couch(pid: str, note: dict[str, Any], lims: str) -> None:
    """
    Write a running note to CouchDB via the genomics status API.

    Args:
        pid: Project ID
        note: Dictionary containing note data with keys:
                - note: The content of the note
                - email: Email of the user adding the note
                - categories: List of note categories for the note
                - note_type: Type of the note (e.g., 'project')
        lims: LIMS system identifier

    Raises:
        SystemExit: If configuration is invalid or API call fails
    """
    config_genstat = "~/config/genstat-conf.yaml"
    with open(os.path.expanduser(config_genstat)) as config_file:
        config: dict[str, Any] = yaml.safe_load(config_file)
    if not config["rn_key"]:
        email_responsible(
            f"Genomics status token credentials not found in {lims}\n ",
            "genomics-bioinfo@scilifelab.se",
        )
        email_responsible(
            f"Running note save for {pid} failed on LIMS! Please contact genomics-bioinfo@scilifelab.se to resolve the issue!",
            note["email"],
        )
        sys.exit(2)

    signed_jwt: str = create_jwt_token(config["rn_key"])
    url = f"{config['genomics-status-url']}/api/v1/running_notes/{pid}"
    result: requests.Response = requests.post(
        url,
        headers={"Authorization": f"Bearer {signed_jwt}"},
        json=note,
    )
    if result.status_code != 201:
        msg = f"Running note save failed from {lims} to {config['genomics-status-url']} for {pid}"
        msg += f"\nStatus code: {result.status_code}\nResponse: {result.text}\n"
        for user_email in ["genomics-bioinfo@scilifelab.se", note["email"]]:
            email_responsible(msg, user_email)
