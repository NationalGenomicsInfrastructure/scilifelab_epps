import datetime
import json
import smtplib
from email.message import Message
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

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
        msg["Subject"] = "[Error] Sync error from LIMS to Genomics Status"
    else:
        msg = MIMEMultipart("alternative")
        msg["Subject"] = (
            subject if subject else "Sync info from LIMS to Genomics Status"
        )

        msg.attach(MIMEText(message, "plain"))
        if html:
            msg.attach(MIMEText(html, "html"))

    msg["From"] = "Lims_monitor"
    msg["To"] = resp_email

    with smtplib.SMTP("localhost") as s:
        s.sendmail("genologics-lims@scilifelab.se", msg["To"], msg.as_string())
