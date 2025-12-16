#!/usr/bin/env python
DESC = """Module called by other EPP scripts to write notes to couchdb
"""
import datetime
import json
import os
import smtplib
import sys
from email.message import Message
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from typing import Any

import requests
import yaml
from jwcrypto import jwk, jwt


def create_jwt_token(key_config: dict[str, str]) -> str:
    """
    Create a JWT token for API authentication.

    Args:
        key_config: Dictionary containing key configuration with keys:
                   - key_path: Path to private key file
                   - owner: Token owner/subject
                   - algorithm: Signing algorithm (e.g., 'ES256')
                   - key_id: Key identifier

    Returns:
        Serialized JWT token as string
    """
    with open(key_config["key_path"], "rb") as f:
        private_key_pem: bytes = f.read()

    private_key = jwk.JWK.from_pem(private_key_pem)
    now = datetime.datetime.now(datetime.timezone.utc)
    expires_at = now + datetime.timedelta(minutes=3)
    claims = {
        "sub": key_config["owner"],
        "iat": int(now.timestamp()),
        "exp": int(expires_at.timestamp()),
    }
    token = jwt.JWT(
        header={"alg": key_config["algorithm"], "kid": key_config["key_id"]},
        claims=json.dumps(claims),
    )
    token.make_signed_token(private_key)
    signed_jwt: str = token.serialize()
    return signed_jwt


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
    if not config["key"]:
        email_responsible(
            f"Genomics status token credentials not found in {lims}\n ",
            "genomics-bioinfo@scilifelab.se",
        )
        email_responsible(
            f"Running note save for {pid} failed on LIMS! Please contact genomics-bioinfo@scilifelab.se to resolve the issue!",
            note["email"],
        )
        sys.exit(2)

    signed_jwt: str = create_jwt_token(config["key"])
    url = f"{config['genomics-status-url']}/api/v1/running_notes/{pid}"
    result: requests.Response = requests.post(
        url,
        headers={"Authorization": f"Bearer {signed_jwt}"},
        json=note,
    )
    if result.status_code != 201:
        msg = f"Running note save failed from {lims} to {config['genomics-status-url']} for {pid}"
        for user_email in ["genomics-bioinfo@scilifelab.se", note["email"]]:
            email_responsible(msg, user_email)


def email_responsible(
    message: str,
    resp_email: str,
    error: bool = True,
    subject: str | None = None,
    html: str | None = None,
) -> None:
    msg: Message
    if error:
        body = "Error: " + message
        body += "\n\n--\nThis is an automatically generated error notification"
        msg = MIMEText(body)
        msg["Subject"] = "[Error] Running note sync error from LIMS to Genomics Status"
    else:
        msg = MIMEMultipart("alternative")
        msg["Subject"] = (
            subject
            if subject
            else "Running note sync info from LIMS to Genomics Status"
        )

        msg.attach(MIMEText(message, "plain"))
        if html:
            msg.attach(MIMEText(html, "html"))

    msg["From"] = "Lims_monitor"
    msg["To"] = resp_email

    with smtplib.SMTP("localhost") as s:
        s.sendmail("genologics-lims@scilifelab.se", msg["To"], msg.as_string())
