#!/usr/bin/env python

import psycopg2
import yaml
from genologics.entities import Researcher


def get_epp_user(lims, procid):
    """
    Get the EPP user from the database.
    Returns:
        str: The EPP user.
    """
    # Setup DB connection
    with open("/opt/gls/clarity/users/glsai/config/genosqlrc.yaml") as f:
        config = yaml.safe_load(f)

    query = (
        "SELECT ps.researcherid "
        "FROM principals ps "
        "JOIN externalprogramstatus eps ON ps.principalid = eps.ownerid "
        "JOIN auditchangelog acl ON acl.rowpk = CAST(eps.externalprogramstatusid AS TEXT) "
        "AND acl.tablename = 'externalprogramstatus' "
        "JOIN auditeventlog ael ON ael.eventid = acl.eventid "
        "WHERE ael.eventtype = 'EPP_EXECUTE' "
        f"AND eps.processid = {procid.split('-')[1]} "
        "AND eps.status = 'RUNNING';"
    )

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
