#!/usr/bin/env python

import psycopg2
import yaml
from genologics.entities import Lims, Researcher


def get_epp_user(
    lims: Lims, procid: str | None = None, project_id: str | None = None
) -> Researcher:
    """
    Get the EPP user from the database.
    Returns:
        Researcher object for whomever launched the EPP
    """
    # Setup DB connection
    with open("/opt/gls/clarity/users/glsai/config/genosqlrc.yaml") as f:
        config = yaml.safe_load(f)

    # When an epp is launched, 2 events are created in auditeventlog, one for event (i.e. advancing in a step,
    # clicking the button and so on) and one for consumption of the next EPP request by the API. Both of them map to the
    # the same rowpk in auditchangelog. This rowpk maps to the externalprogramstatusid in the table externalprogramstatus.
    # So by checking for 'EPP_CONSUME_NEXT_REQUEST' and the process thats currently running we should be able to get
    # the user that launched the EPP.
    query = """--sql
        SELECT ps.researcherid
        FROM principals ps
        JOIN externalprogramstatus eps ON ps.principalid = eps.ownerid
        JOIN auditchangelog acl ON acl.rowpk = CAST(eps.externalprogramstatusid AS TEXT)
        AND acl.tablename = 'externalprogramstatus'
        JOIN auditeventlog ael ON ael.eventid = acl.eventid
    """
    if project_id:
        query += """--sql
        JOIN project_programstatus_map ppm ON ppm.externalprogramstatusid = eps.externalprogramstatusid
        """
    query += """--sql
        WHERE ael.eventtype ='EPP_CONSUME_NEXT_REQUEST'
        AND eps.status = 'RUNNING'
     """
    if procid:
        query += f"""--sql
        AND eps.processid = {procid.split("-")[1]};
        """
    if project_id:
        query += f"""--sql
        AND ppm.projectid = {project_id.split("P")[1]};
        """

    with psycopg2.connect(
        user=config["username"],
        host=config["url"],
        database=config["db"],
        password=config["password"],
    ) as connection:
        with connection.cursor() as cursor:
            cursor.execute(query)
            query_output = cursor.fetchone()
            if query_output is None:
                raise ValueError(f"No EPP user found for process ID {procid}")

    technician_id = query_output[0]
    return Researcher(lims, id=str(technician_id))
