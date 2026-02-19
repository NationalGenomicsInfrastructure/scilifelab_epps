#!/usr/bin/env python

import json
import os
import re
import sys
from argparse import ArgumentParser
from datetime import datetime
from io import StringIO
from typing import Any

import pandas as pd
import requests
import yaml
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from data.Chromium_10X_indexes import Chromium_10X_indexes
from scilifelab_epps.utils.genstat_conn import create_jwt_token, email_responsible

# Load SS3 indexes
SMARTSEQ3_indexes_json = (
    "/opt/gls/clarity/users/glsai/repos/scilifelab_epps/data/SMARTSEQ3_indexes.json"
)
with open(SMARTSEQ3_indexes_json) as file:
    SMARTSEQ3_indexes = json.loads(file.read())

DESC = """EPP used to create samplesheets for Illumina sequencing platforms"""

# Pre-compile regexes in global scope:
IDX_PAT = re.compile("([ATCG]{4,}N*)-?([ATCG]*)")
TENX_SINGLE_PAT = re.compile("SI-(?:GA|NA)-[A-H][1-9][0-2]?")
TENX_DUAL_PAT = re.compile("SI-(?:TT|NT|NN|TN|TS)-[A-H][1-9][0-2]?")
SMARTSEQ_PAT = re.compile("SMARTSEQ[1-9]?-[1-9][0-9]?[A-P]")
NGISAMPLE_PAT = re.compile("P[0-9]+_[0-9]+")
SEQSETUP_PAT = re.compile("[0-9]+-[0-9A-z]+-[0-9A-z]+-[0-9]+")

compl = {"A": "T", "C": "G", "G": "C", "T": "A"}


def check_index_distance(data, log):
    lanes = {x["lane"] for x in data}
    for l in lanes:
        indexes = [
            x.get("index_1", "") + x.get("index_2", "") for x in data if x["lane"] == l
        ]
        if not indexes or len(indexes) == 1:
            return None
        for i, b in enumerate(indexes[:-1]):
            start = i + 1
            for b2 in indexes[start:]:
                d = my_distance(b, b2)
                if d < 2 and not is_special_idx(b) and not is_special_idx(b2):
                    log.append(
                        f"Found indexes {b} and {b2} in lane {l}, indexes are too close"
                    )


def is_special_idx(idx_name):
    if (
        TENX_DUAL_PAT.findall(idx_name)
        or TENX_SINGLE_PAT.findall(idx_name)
        or SMARTSEQ_PAT.findall(idx_name)
        or idx_name == "NoIndex"
    ):
        return True
    else:
        return False


def my_distance(idx1, idx2):
    short = min((idx1, idx2), key=len)
    lon = idx1 if short == idx2 else idx2

    diffs = 0
    for i, c in enumerate(short):
        if c != lon[i]:
            diffs += 1
    return diffs


def gen_NovaSeqXPlus_lane_data(pro):
    data = []
    lanes = set()
    header_ar = [
        "FCID",
        "Lane",
        "Sample_ID",
        "Sample_Name",
        "Sample_Ref",
        "index",
        "index2",
        "Description",
        "Control",
        "Recipe",
        "Operator",
        "Sample_Project",
    ]
    for out in pro.all_outputs():
        if out.type == "Analyte":
            for sample in out.samples:
                sample_idxs = set()
                find_barcode(sample_idxs, sample, pro)
                for idxs in sample_idxs:
                    sp_obj = {}
                    sp_obj["lane"] = out.location[1].split(":")[0].replace(",", "")
                    lanes.add(sp_obj["lane"])
                    if NGISAMPLE_PAT.findall(sample.name):
                        sp_obj["sample_id"] = f"Sample_{sample.name}".replace(",", "")
                        sp_obj["sample_name"] = sample.name.replace(",", "")
                        sp_obj["description"] = sample.project.name.replace(
                            ".", "__"
                        ).replace(",", "")
                        sp_obj["sample_ref"] = sample.project.udf.get(
                            "Reference genome", ""
                        ).replace(",", "")
                        seq_setup = sample.project.udf.get("Sequencing setup", "")
                        if SEQSETUP_PAT.findall(seq_setup):
                            sp_obj["rc"] = "{}-{}".format(
                                seq_setup.split("-")[0], seq_setup.split("-")[3]
                            )
                            sp_obj["recipe"] = seq_setup
                        else:
                            sp_obj["rc"] = "0-0"
                            sp_obj["recipe"] = "0-0-0-0"
                    else:
                        sp_obj["sample_id"] = (
                            f"Sample_{sample.name}".replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sample_name"] = (
                            sample.name.replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["description"] = "Control"
                        sp_obj["sample_ref"] = "Control"
                        sp_obj["rc"] = "0-0"
                        sp_obj["recipe"] = "0-0-0-0"
                    sp_obj["control"] = "N"
                    sp_obj["operator"] = pro.technician.name.replace(" ", "_").replace(
                        ",", ""
                    )
                    sp_obj["flowcell_id"] = (
                        out.location[0].name.replace(",", "").upper()
                    )
                    sp_obj["sw"] = out.location[1].replace(",", "")
                    sp_obj["index_1"] = idxs[0].replace(",", "").upper()
                    if idxs[1]:
                        sp_obj["index_2"] = idxs[1].replace(",", "").upper()
                    else:
                        sp_obj["index_2"] = ""
                    data.append(sp_obj)
    header = "{}\n".format(",".join(header_ar))
    str_data = ""
    for line in sorted(data, key=lambda x: x["lane"]):
        l_data = [
            line["flowcell_id"],
            line["lane"],
            line["sample_id"],
            line["sample_name"],
            line["sample_ref"],
            line["index_1"],
            line["index_2"],
            line["description"],
            line["control"],
            line["rc"],
            line["operator"],
            line["description"],
        ]
        str_data = str_data + ",".join(l_data) + "\n"

    content = f"{header}{str_data}"
    df = pd.read_csv(StringIO(content))
    df = df.sort_values(["Lane", "Sample_ID"])
    content = df.to_csv(index=False)

    return (content, data, len(lanes))


def gen_Miseq_header(pro, chem):
    project_name = pro.all_inputs()[0].samples[0].project.name
    header = "[Header]\nInvestigator Name,{inn}\nProject Name,{pn}\nExperiment Name,{en}\nDate,{dt}\nWorkflow,{wf}\nModule,{mod}\nAssay,{ass}\nDescription,{dsc}\nChemistry,{chem}\n".format(
        inn=pro.technician.name,
        pn=project_name,
        en=pro.udf["Flowcell ID"],
        dt=datetime.now().strftime("%Y-%m-%d"),
        wf=pro.udf["Workflow"],
        mod=pro.udf["Module"],
        ass="null",
        dsc=pro.udf["Description"],
        chem=chem,
    )
    return header


def gen_Miseq_reads(pro):
    reads = "[Reads]\n"
    if pro.udf["Read 1 Cycles"]:
        reads = reads + "{}\n".format(pro.udf["Read 1 Cycles"])
    if pro.udf.get("Read 2 Cycles"):
        reads = reads + "{}\n".format(pro.udf["Read 2 Cycles"])
    else:
        reads = reads + "0\n"
    return reads


def gen_Miseq_settings(pro):
    ogf = 1 if pro.udf["OnlyGenerateFASTQ"] else 0
    fpdcrd = 1 if pro.udf["FilterPCRDuplicates"] else 0
    custom_r1_primer = (
        "CustomRead1PrimerMix,C1\n" if pro.udf["CustomRead1PrimerMix"] else ""
    )
    custom_index_primer = (
        "CustomIndexPrimerMix,C2\n" if pro.udf["CustomIndexPrimerMix"] else ""
    )
    custom_r2_primer = (
        "CustomRead2PrimerMix,C3\n" if pro.udf["CustomRead2PrimerMix"] else ""
    )
    settings = f"[Settings]\nOnlyGenerateFASTQ,{ogf}\nFilterPCRDuplicates,{fpdcrd}\n{custom_r1_primer}{custom_index_primer}{custom_r2_primer}"
    return settings


def is_key_empty_in_all_dicts(key, list_of_dicts):
    for dictionary in list_of_dicts:
        if key not in dictionary or dictionary[key] != "":
            return False
    return True


def gen_Miseq_data(pro):
    chem = "amplicon"
    data = []
    lanes = set()
    header_ar = [
        "FCID",
        "Lane",
        "Sample_ID",
        "Sample_Name",
        "Sample_Ref",
        "index",
        "index2",
        "Description",
        "Control",
        "Recipe",
        "Operator",
        "Sample_Project",
    ]
    key_order = [
        "flowcell_id",
        "lane",
        "sample_name",
        "sample_name",
        "ref",
        "index_1",
        "index_2",
        "description",
        "control",
        "rc",
        "operator",
        "description",
    ]
    for out in pro.all_outputs():
        if out.type == "Analyte":
            for sample in out.samples:
                sample_idxs = set()
                find_barcode(sample_idxs, sample, pro)
                for idxs in sample_idxs:
                    sp_obj = {}
                    sp_obj["lane"] = "1"
                    lanes.add(sp_obj["lane"])
                    if NGISAMPLE_PAT.findall(sample.name):
                        sp_obj["sample_id"] = f"Sample_{sample.name}".replace(",", "")
                        sp_obj["sample_name"] = sample.name.replace(",", "")
                        sp_obj["description"] = sample.project.name.replace(
                            ".", "__"
                        ).replace(",", "")
                        sp_obj["sample_ref"] = sample.project.udf.get(
                            "Reference genome", ""
                        ).replace(",", "")
                        seq_setup = sample.project.udf.get("Sequencing setup", "")
                        pj_type = (
                            "by user"
                            if sample.project.udf["Library construction method"]
                            == "Finished library (by user)"
                            else "inhouse"
                        )
                        if SEQSETUP_PAT.findall(seq_setup):
                            sp_obj["rc"] = "{}-{}".format(
                                seq_setup.split("-")[0], seq_setup.split("-")[3]
                            )
                            sp_obj["recipe"] = seq_setup
                        else:
                            sp_obj["rc"] = "0-0"
                            sp_obj["recipe"] = "0-0-0-0"
                    else:
                        sp_obj["sid"] = (
                            f"Sample_{sample.name}".replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sample_name"] = (
                            sample.name.replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["description"] = "Control"
                        sp_obj["sample_ref"] = "Control"
                        sp_obj["rc"] = "0-0"
                        sp_obj["recipe"] = "0-0-0-0"
                        pj_type = "Control"
                    sp_obj["control"] = "N"
                    sp_obj["operator"] = pro.technician.name.replace(" ", "_").replace(
                        ",", ""
                    )
                    sp_obj["flowcell_id"] = (
                        out.location[0].name.replace(",", "").upper()
                    )
                    sp_obj["sw"] = out.location[1].replace(",", "")

                    # Expand 10X single indexes
                    if TENX_SINGLE_PAT.findall(idxs[0]):
                        for tenXidx in Chromium_10X_indexes[
                            TENX_SINGLE_PAT.findall(idxs[0])[0]
                        ]:
                            sp_obj_sub = {}
                            sp_obj_sub["lane"] = sp_obj["lane"]
                            sp_obj_sub["sample_id"] = sp_obj["sample_id"]
                            sp_obj_sub["sample_name"] = sp_obj["sample_name"]
                            sp_obj_sub["description"] = sp_obj["description"]
                            sp_obj_sub["sample_ref"] = sp_obj["sample_ref"]
                            sp_obj_sub["rc"] = sp_obj["rc"]
                            sp_obj_sub["control"] = sp_obj["control"]
                            sp_obj_sub["operator"] = sp_obj["operator"]
                            sp_obj_sub["flowcell_id"] = sp_obj["flowcell_id"]
                            sp_obj_sub["sw"] = sp_obj["sw"]
                            sp_obj_sub["index_1"] = tenXidx.replace(",", "")
                            sp_obj_sub["index_2"] = ""
                            data.append(sp_obj_sub)
                    # Case of 10X dual indexes
                    elif TENX_DUAL_PAT.findall(idxs[0]):
                        sp_obj["index_1"] = Chromium_10X_indexes[
                            TENX_DUAL_PAT.findall(idxs[0])[0]
                        ][0].replace(",", "")
                        sp_obj["index_2"] = "".join(
                            reversed(
                                [
                                    compl.get(b, b)
                                    for b in Chromium_10X_indexes[
                                        TENX_DUAL_PAT.findall(idxs[0])[0]
                                    ][1]
                                    .replace(",", "")
                                    .upper()
                                ]
                            )
                        )
                        data.append(sp_obj)
                    # Case of SS3 indexes
                    elif SMARTSEQ_PAT.findall(idxs[0]):
                        for i7_idx in SMARTSEQ3_indexes[idxs[0]][0]:
                            for i5_idx in SMARTSEQ3_indexes[idxs[0]][1]:
                                sp_obj_sub = {}
                                sp_obj_sub["lane"] = sp_obj["lane"]
                                sp_obj_sub["sample_id"] = sp_obj["sample_id"]
                                sp_obj_sub["sample_name"] = sp_obj["sample_name"]
                                sp_obj_sub["description"] = sp_obj["description"]
                                sp_obj_sub["sample_ref"] = sp_obj["sample_ref"]
                                sp_obj_sub["rc"] = sp_obj["rc"]
                                sp_obj_sub["control"] = sp_obj["control"]
                                sp_obj_sub["operator"] = sp_obj["operator"]
                                sp_obj_sub["flowcell_id"] = sp_obj["flowcell_id"]
                                sp_obj_sub["sw"] = sp_obj["sw"]
                                sp_obj_sub["index_1"] = i7_idx
                                sp_obj_sub["index_2"] = "".join(
                                    reversed(
                                        [
                                            compl.get(b, b)
                                            for b in i5_idx.replace(",", "").upper()
                                        ]
                                    )
                                )
                                data.append(sp_obj_sub)
                    # NoIndex cases
                    elif idxs[0].replace(",", "").upper() == "NOINDEX" or (
                        idxs[0].replace(",", "").upper() == ""
                        and idxs[1].replace(",", "").upper() == ""
                    ):
                        sp_obj["index_1"] = ""
                        sp_obj["index_2"] = ""
                        data.append(sp_obj)
                    # Ordinary indexes
                    else:
                        sp_obj["index_1"] = idxs[0].replace(",", "").upper()
                        if idxs[1]:
                            if pj_type == "by user":
                                sp_obj["index_2"] = idxs[1].replace(",", "").upper()
                            else:
                                sp_obj["index_2"] = "".join(
                                    reversed(
                                        [
                                            compl.get(b, b)
                                            for b in idxs[1].replace(",", "").upper()
                                        ]
                                    )
                                )
                        else:
                            sp_obj["index_2"] = ""
                        data.append(sp_obj)

    if is_key_empty_in_all_dicts("index_1", data):
        header_ar.remove("index")
        key_order.remove("index_1")
        chem = "Default"
    if is_key_empty_in_all_dicts("index_2", data):
        header_ar.remove("index2")
        key_order.remove("index_2")
        chem = "Default"

    header = "{}\n".format(",".join(header_ar))
    str_data = ""
    for line in sorted(data, key=lambda x: x["lane"]):
        l_data = []
        for key in key_order:
            l_data.append(line[key])
        str_data = str_data + ",".join(l_data) + "\n"

    content = f"{header}{str_data}"
    df = pd.read_csv(StringIO(content))
    df = df.sort_values(["Lane", "Sample_ID"])
    content = df.to_csv(index=False)
    content = f"[Data]\n{content}\n"

    return (content, data, chem, len(lanes))


def gen_Nextseq_lane_data(pro):
    data = []
    lanes = set()
    header_ar = [
        "FCID",
        "Lane",
        "Sample_ID",
        "Sample_Name",
        "Sample_Ref",
        "index",
        "index2",
        "Description",
        "Control",
        "Recipe",
        "Operator",
        "Sample_Project",
    ]
    for out in pro.all_outputs():
        if out.type == "Analyte":
            for sample in out.samples:
                sample_idxs = set()
                find_barcode(sample_idxs, sample, pro)
                for idxs in sample_idxs:
                    sp_obj = {}
                    sp_obj["lane"] = out.location[1].split(":")[0].replace(",", "")
                    lanes.add(sp_obj["lane"])
                    if NGISAMPLE_PAT.findall(sample.name):
                        sp_obj["sample_id"] = f"Sample_{sample.name}".replace(",", "")
                        sp_obj["sample_name"] = sample.name.replace(",", "")
                        sp_obj["sample_project"] = sample.project.name.replace(
                            ".", "__"
                        ).replace(",", "")
                        sp_obj["sample_ref"] = sample.project.udf.get(
                            "Reference genome", ""
                        ).replace(",", "")
                        seq_setup = sample.project.udf.get("Sequencing setup", "")
                        if SEQSETUP_PAT.findall(seq_setup):
                            sp_obj["rc"] = "{}-{}".format(
                                seq_setup.split("-")[0], seq_setup.split("-")[3]
                            )
                            sp_obj["recipe"] = seq_setup
                        else:
                            sp_obj["rc"] = "0-0"
                            sp_obj["recipe"] = "0-0-0-0"
                    else:
                        sp_obj["sample_id"] = (
                            f"Sample_{sample.name}".replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sample_name"] = (
                            sample.name.replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sample_project"] = "Control"
                        sp_obj["sample_ref"] = "Control"
                        sp_obj["rc"] = "0-0"
                        sp_obj["recipe"] = "0-0-0-0"
                    sp_obj["control"] = "N"
                    sp_obj["operator"] = pro.technician.name.replace(" ", "_").replace(
                        ",", ""
                    )
                    sp_obj["flowcell_id"] = (
                        out.location[0].name.replace(",", "").upper()
                    )
                    sp_obj["sw"] = out.location[1].replace(",", "")
                    sp_obj["index_1"] = idxs[0].replace(",", "")
                    if idxs[1]:
                        sp_obj["index_2"] = idxs[1].replace(",", "").upper()
                    else:
                        sp_obj["index_2"] = ""
                    data.append(sp_obj)
    header = "{}\n".format(",".join(header_ar))
    str_data = ""
    for line in sorted(data, key=lambda x: x["lane"]):
        l_data = [
            line["flowcell_id"],
            line["lane"],
            line["sample_name"],
            line["sample_name"],
            line["sample_ref"],
            line["index_1"],
            line["index_2"],
            line["sample_project"],
            line["control"],
            line["rc"],
            line["operator"],
            line["sample_project"],
        ]
        str_data = str_data + ",".join(l_data) + "\n"

    content = f"{header}{str_data}"
    df = pd.read_csv(StringIO(content))
    df = df.sort_values(["Lane", "Sample_ID"])
    content = df.to_csv(index=False)

    return (content, data, len(lanes))


def find_barcode(sample_idxs, sample, process):
    # print "trying to find {} barcode in {}".format(sample.name, process.name)
    for art in process.all_inputs():
        if sample in art.samples:
            if len(art.samples) == 1 and art.reagent_labels:
                # In rare cases we have a pool containing the same sample with different labels
                for reagent_label in art.reagent_labels:
                    reagent_label_name = reagent_label.upper().replace(" ", "")
                    idxs = (
                        TENX_SINGLE_PAT.findall(reagent_label_name)
                        or TENX_DUAL_PAT.findall(reagent_label_name)
                        or SMARTSEQ_PAT.findall(reagent_label_name)
                    )
                    if idxs:
                        # Put in tuple with empty string as second index to
                        # match expected type:
                        sample_idxs.add((idxs[0], ""))
                    else:
                        try:
                            idxs = IDX_PAT.findall(reagent_label_name)[0]
                            sample_idxs.add(idxs)
                        except IndexError:
                            try:
                                # we only have the reagent label name.
                                rt = lims.get_reagent_types(name=reagent_label_name)[0]
                                idxs = IDX_PAT.findall(rt.sequence)[0]
                                sample_idxs.add(idxs)
                            except:
                                sample_idxs.add(("NoIndex", ""))
            else:
                if art == sample.artifact or not art.parent_process:
                    pass
                else:
                    find_barcode(sample_idxs, sample, art.parent_process)


def upload_to_genstat(data, metadata, fc_name):
    config_genstat = "~/config/genstat-conf.yaml"
    with open(os.path.expanduser(config_genstat)) as config_file:
        config: dict[str, Any] = yaml.safe_load(config_file)
    if not config["samplesheet_key"]:
        email_responsible(
            f"Genomics status token credentials not found in {lims}\n Samplesheet upload for {fc_name} failed on LIMS! Please contact genomics-bioinfo@scilifelab.se to resolve the issue!",
            "genomics-bioinfo@scilifelab.se",
        )
        sys.exit(2)

    signed_jwt: str = create_jwt_token(config["samplesheet_key"])
    genstat_sample_info_url = (
        f"{config['genomics-status-url']}/api/v1/demux_sample_info/{fc_name}"
    )
    result: requests.Response = requests.post(
        genstat_sample_info_url,
        headers={"Authorization": f"Bearer {signed_jwt}"},
        json={"uploaded_lims_info": data, "metadata": metadata},
    )

    if result.status_code != 201:
        msg = f"Samplesheet upload failed from {lims} to {config['genomics-status-url']} for {fc_name}"
        msg += f"\nStatus code: {result.status_code}\nResponse: {result.text}\n"
        email_responsible(msg, "genomics-bioinfo@scilifelab.se")


def test():
    log = []
    d = [
        {"lane": 1, "index_1": "ATTT", "index_2": ""},
        {"lane": 1, "index_1": "ATCTATCG", "index_2": ""},
        {"lane": 1, "index_1": "ATCG", "index_2": "ATCG"},
    ]
    check_index_distance(d, log)
    print(log)


def main(lims, args):
    log = []
    num_lanes = 0
    thisyear = datetime.now().year
    content = None
    if args.mytest:
        test()
    else:
        process = Process(lims, id=args.pid)

        if "Load to Flowcell (NovaSeqXPlus)" in process.type.name:
            (content, obj, num_lanes) = gen_NovaSeqXPlus_lane_data(process)
            check_index_distance(obj, log)
            if os.path.exists(f"/srv/ngi-nas-ns/samplesheets/NovaSeqXPlus/{thisyear}"):
                try:
                    with open(
                        "/srv/ngi-nas-ns/samplesheets/NovaSeqXPlus/{}/{}.csv".format(
                            thisyear, obj[0]["flowcell_id"]
                        ),
                        "w",
                    ) as sf:
                        sf.write(content)
                except Exception as e:
                    log.append(str(e))

        elif process.type.name == "Denature, Dilute and Load Sample (MiSeq) 4.0":
            reads = gen_Miseq_reads(process)
            settings = gen_Miseq_settings(process)
            (content, obj, chem, num_lanes) = gen_Miseq_data(process)
            header = gen_Miseq_header(process, chem)
            check_index_distance(obj, log)
            content = f"{header}{reads}{settings}{content}"

        elif process.type.name == "Load to Flowcell (NextSeq v1.0)":
            (content, obj, num_lanes) = gen_Nextseq_lane_data(process)
            check_index_distance(obj, log)
            nextseq_fc = (
                process.udf["Flowcell Series Number"]
                if process.udf["Flowcell Series Number"]
                else obj[0]["flowcell_id"]
            )
            if os.path.exists(f"/srv/ngi-nas-ns/samplesheets/nextseq/{thisyear}"):
                try:
                    with open(
                        f"/srv/ngi-nas-ns/samplesheets/nextseq/{thisyear}/{nextseq_fc}.csv",
                        "w",
                    ) as sf:
                        sf.write(content)
                except Exception as e:
                    log.append(str(e))

        elif process.type.name == "Load to Flowcell (MiSeq i100) v1.0":
            (content, obj, num_lanes) = gen_Nextseq_lane_data(process)
            check_index_distance(obj, log)
            miseqi100_fc = (
                process.udf["Flowcell Series Number"]
                if process.udf["Flowcell Series Number"]
                else obj[0]["fc"]
            )
            if os.path.exists(f"/srv/ngi-nas-ns/samplesheets/MiSeqi100/{thisyear}"):
                try:
                    with open(
                        f"/srv/ngi-nas-ns/samplesheets/MiSeqi100/{thisyear}/{miseqi100_fc}.csv",
                        "w",
                    ) as sf:
                        sf.write(content)
                except Exception as e:
                    log.append(str(e))

        if not args.test:
            for out in process.all_outputs():
                if out.name == "Scilifelab SampleSheet":
                    ss_art = out
                elif out.name == "Scilifelab Log":
                    log_art = out
                elif out.type == "Analyte":
                    if process.type.name == "Load to Flowcell (NextSeq v1.0)":
                        fc_name = (
                            process.udf["Flowcell Series Number"].upper()
                            if process.udf["Flowcell Series Number"]
                            else out.location[0].name.upper()
                        )
                    else:
                        fc_name = out.location[0].name.upper()
                elif process.type.name in [
                    "MinION QC",
                    "Load Sample and Sequencing (MinION) 1.0",
                ]:
                    run_type = "QC" if process.type.name == "MinION QC" else "DELIVERY"
                    fc_name = (
                        run_type
                        + "_"
                        + process.udf["Nanopore Kit"]
                        + "_"
                        + process.udf["Flowcell ID"].upper()
                        + "_"
                        + "Samplesheet"
                        + "_"
                        + process.id
                    )
                else:
                    fc_name = "Samplesheet" + "_" + process.id

            with open(f"{fc_name}.csv", "w", 0o664) as f:
                f.write(content)
            os.chmod(f"{fc_name}.csv", 0o664)
            # Upload samplesheet to CouchDB through Genstat
            run_setup = f"{process.udf.get('Read 1 Cycles')}_{process.udf.get('Index Read 1', 'x')}_{process.udf.get('Index Read 2', 'x')}_{process.udf.get('Read 2 Cycles')}"

            # Determine instrument type from process name
            instrument_type_mapping = {
                "Load to Flowcell (NovaSeqXPlus)": "NovaSeqXPlus",
                "Denature, Dilute and Load Sample (MiSeq) 4.0": "MiSeq",
                "Load to Flowcell (NextSeq v1.0)": "NextSeq",
                "Load to Flowcell (MiSeq i100) v1.0": "MiSeq i100",
            }
            instrument_type = instrument_type_mapping.get(process.type.name, "")

            metadata = {
                "num_lanes": num_lanes,
                "run_setup": run_setup,
                "setup_lims_step_id": process.id,
                "instrument_type": instrument_type,
                "run_mode": process.udf.get("Run Mode", ""),
            }
            # Check that content exists to upload obj
            if content:
                upload_to_genstat(obj, metadata, fc_name)
            for f in ss_art.files:
                lims.request_session.delete(f.uri)
            lims.upload_new_file(ss_art, f"{fc_name}.csv")
            if log:
                with open(f"{log_art.id}_{fc_name}_Error.log", "w") as f:
                    f.write("\n".join(log))
                # Upload log to file slot
                lims.upload_new_file(log_art, f"{log_art.id}_{fc_name}_Error.log")

                sys.stderr.write("Errors were met, check the log.")
                sys.exit(2)

        else:
            print(content)
            print(log)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    parser.add_argument(
        "--test", action="store_true", help="do not upload the samplesheet"
    )
    parser.add_argument("--mytest", action="store_true", help="mytest")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
