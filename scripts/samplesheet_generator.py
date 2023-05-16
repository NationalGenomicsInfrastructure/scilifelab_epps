#!/usr/bin/env python

import re
import os
import sys
import pandas as pd

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

from argparse import ArgumentParser
from datetime import datetime
from genologics.lims import Lims
from genologics.entities import Process
from genologics.config import BASEURI, USERNAME, PASSWORD

from data.Chromium_10X_indexes import Chromium_10X_indexes

DESC = """EPP used to create samplesheets for Illumina sequencing platforms"""

# Pre-compile regexes in global scope:
IDX_PAT = re.compile("([ATCG]{4,}N*)-?([ATCG]*)")
TENX_SINGLE_PAT = re.compile("SI-(?:GA|NA)-[A-H][1-9][0-2]?")
TENX_DUAL_PAT = re.compile("SI-(?:TT|NT|NN|TN|TS)-[A-H][1-9][0-2]?")
SMARTSEQ_PAT = re.compile("SMARTSEQ[1-9]?-[1-9][0-9]?[A-P]")
NGISAMPLE_PAT = re.compile("P[0-9]+_[0-9]+")
SEQSETUP_PAT = re.compile("[0-9]+-[0-9A-z]+-[0-9A-z]+-[0-9]+")


def check_index_distance(data, log):
    lanes = set([x["lane"] for x in data])
    for l in lanes:
        indexes = [
            x.get("idx1", "") + x.get("idx2", "") for x in data if x["lane"] == l
        ]
        if not indexes or len(indexes) == 1:
            return None
        for i, b in enumerate(indexes[:-1]):
            start = i + 1
            for b2 in indexes[start:]:
                d = my_distance(b, b2)
                if d < 2:
                    log.append(
                        "Found indexes {} and {} in lane {}, indexes are too close".format(
                            b, b2, l
                        )
                    )


def my_distance(idx1, idx2):
    short = min((idx1, idx2), key=len)
    lon = idx1 if short == idx2 else idx2

    diffs = 0
    for i, c in enumerate(short):
        if c != lon[i]:
            diffs += 1
    return diffs


def gen_Novaseq_lane_data(pro):
    data = []
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
                    if NGISAMPLE_PAT.findall(sample.name):
                        sp_obj["sid"] = "Sample_{}".format(sample.name).replace(",", "")
                        sp_obj["sn"] = sample.name.replace(",", "")
                        sp_obj["pj"] = sample.project.name.replace(".", "__").replace(
                            ",", ""
                        )
                        sp_obj["ref"] = sample.project.udf.get(
                            "Reference genome", ""
                        ).replace(",", "")
                        seq_setup = sample.project.udf.get("Sequencing setup", "")
                        if SEQSETUP_PAT.findall(seq_setup):
                            sp_obj["rc"] = "{}-{}".format(
                                seq_setup.split("-")[0], seq_setup.split("-")[3]
                            )
                        else:
                            sp_obj["rc"] = "0-0"
                    else:
                        sp_obj["sid"] = (
                            "Sample_{}".format(sample.name)
                            .replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sn"] = (
                            sample.name.replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["pj"] = "Control"
                        sp_obj["ref"] = "Control"
                        sp_obj["rc"] = "0-0"
                    sp_obj["ct"] = "N"
                    sp_obj["op"] = pro.technician.name.replace(" ", "_").replace(
                        ",", ""
                    )
                    sp_obj["fc"] = out.location[0].name.replace(",", "")
                    sp_obj["sw"] = out.location[1].replace(",", "")
                    sp_obj["idx1"] = idxs[0].replace(",", "").upper()
                    if idxs[1]:
                        if pro.udf["Reagent Version"] == "v1.5":
                            sp_obj["idx2"] = idxs[1].replace(",", "").upper()
                        elif pro.udf["Reagent Version"] == "v1.0":
                            compl = {"A": "T", "C": "G", "G": "C", "T": "A"}
                            sp_obj["idx2"] = "".join(
                                reversed(
                                    [
                                        compl.get(b, b)
                                        for b in idxs[1].replace(",", "").upper()
                                    ]
                                )
                            )
                    else:
                        sp_obj["idx2"] = ""
                    data.append(sp_obj)
    header = "{}\n".format(",".join(header_ar))
    str_data = ""
    for line in sorted(data, key=lambda x: x["lane"]):
        l_data = [
            line["fc"],
            line["lane"],
            line["sn"],
            line["sn"],
            line["ref"],
            line["idx1"],
            line["idx2"],
            line["pj"],
            line["ct"],
            line["rc"],
            line["op"],
            line["pj"],
        ]
        str_data = str_data + ",".join(l_data) + "\n"

    content = "{}{}".format(header, str_data)
    df = pd.read_csv(StringIO(content))
    df = df.sort_values(["Lane", "Sample_ID"])
    content = df.to_csv(index=False)

    return (content, data)


def gen_Miseq_header(pro):
    project_name = pro.all_inputs()[0].samples[0].project.name
    chem = "Default"
    for io in pro.input_output_maps:
        sample_idxs = set()
        find_barcode(sample_idxs, io[1]["uri"].samples[0], pro)
        idxs = list(sample_idxs)[0]
        if len(idxs) == 2:
            chem = "amplicon"

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
    if pro.udf["Read 2 Cycles"]:
        reads = reads + "{}\n".format(pro.udf["Read 2 Cycles"])
    return reads


def gen_Miseq_settings(pro):
    ogf = 1 if pro.udf["OnlyGenerateFASTQ"] else 0
    fpdcrd = 1 if pro.udf["FilterPCRDuplicates"] else 0
    settings = (
        "[Settings]\nOnlyGenerateFASTQ,{ogf}\nFilterPCRDuplicates,{fpdcrd}\n".format(
            ogf=ogf, fpdcrd=fpdcrd
        )
    )
    return settings


def gen_Miseq_data(pro):
    data = []
    dualindex = False
    noindex = False
    header_ar = [
        "Sample_ID",
        "Sample_Name",
        "Sample_Plate",
        "Sample_Well",
        "Sample_Project",
        "index",
        "I7_Index_ID",
        "index2",
        "I5_Index_ID",
        "Description",
        "GenomeFolder",
    ]
    for io in pro.input_output_maps:
        out = io[1]["uri"]
        if out.type != "Analyte":
            continue
        for sample in out.samples:
            sample_idxs = set()
            find_barcode(sample_idxs, sample, pro)
            if not sample_idxs:
                noindex = True
                header_ar.remove("index")
                header_ar.remove("I7_Index_ID")
                header_ar.remove("index2")
                header_ar.remove("I5_Index_ID")
                sp_obj = {}
                pj_type = ""
                sp_obj["lane"] = "1"
                if NGISAMPLE_PAT.findall(sample.name):
                    sp_obj["sid"] = "Sample_{}".format(sample.name).replace(",", "")
                    sp_obj["sn"] = sample.name.replace(",", "")
                    sp_obj["pj"] = sample.project.name.replace(".", "_").replace(
                        ",", ""
                    )
                    pj_type = (
                        "by user"
                        if sample.project.udf["Library construction method"]
                        == "Finished library (by user)"
                        else "inhouse"
                    )
                else:
                    sp_obj["sid"] = (
                        "Sample_{}".format(sample.name)
                        .replace("(", "")
                        .replace(")", "")
                        .replace(".", "")
                        .replace(" ", "_")
                    )
                    sp_obj["sn"] = (
                        sample.name.replace("(", "")
                        .replace(")", "")
                        .replace(".", "")
                        .replace(" ", "_")
                    )
                    sp_obj["pj"] = "Control"
                    pj_type = "Control"
                sp_obj["fc"] = "{}-{}".format(
                    io[0]["uri"].location[0].name.replace(",", ""),
                    out.location[1].replace(":", ""),
                )
                sp_obj["sw"] = "A1"
                sp_obj["gf"] = pro.udf["GenomeFolder"].replace(",", "")
                data.append(sp_obj)
            else:
                for idxs in sample_idxs:
                    sp_obj = {}
                    pj_type = ""
                    sp_obj["lane"] = "1"
                    if NGISAMPLE_PAT.findall(sample.name):
                        sp_obj["sid"] = "Sample_{}".format(sample.name).replace(",", "")
                        sp_obj["sn"] = sample.name.replace(",", "")
                        sp_obj["pj"] = sample.project.name.replace(".", "_").replace(
                            ",", ""
                        )
                        pj_type = (
                            "by user"
                            if sample.project.udf["Library construction method"]
                            == "Finished library (by user)"
                            else "inhouse"
                        )
                    else:
                        sp_obj["sid"] = (
                            "Sample_{}".format(sample.name)
                            .replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sn"] = (
                            sample.name.replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["pj"] = "Control"
                        pj_type = "Control"
                    sp_obj["fc"] = "{}-{}".format(
                        io[0]["uri"].location[0].name.replace(",", ""),
                        out.location[1].replace(":", ""),
                    )
                    sp_obj["sw"] = "A1"
                    sp_obj["gf"] = pro.udf["GenomeFolder"].replace(",", "")

                    if TENX_DUAL_PAT.findall(idxs[0]):
                        dualindex = True
                        sp_obj["idx1"] = Chromium_10X_indexes[
                            TENX_DUAL_PAT.findall(idxs[0])[0]
                        ][0].replace(",", "")
                        sp_obj["idx1ref"] = Chromium_10X_indexes[
                            TENX_DUAL_PAT.findall(idxs[0])[0]
                        ][0].replace(",", "")
                        compl = {"A": "T", "C": "G", "G": "C", "T": "A"}
                        sp_obj["idx2"] = "".join(
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
                        sp_obj["idx2ref"] = "".join(
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
                    elif TENX_SINGLE_PAT.findall(idxs[0]):
                        if "index2" in header_ar and "I5_Index_ID" in header_ar:
                            header_ar.remove("index2")
                            header_ar.remove("I5_Index_ID")
                        for tenXidx in Chromium_10X_indexes[
                            TENX_SINGLE_PAT.findall(idxs[0])[0]
                        ]:
                            sp_obj_sub = {}
                            sp_obj_sub["lane"] = sp_obj["lane"]
                            sp_obj_sub["sid"] = sp_obj["sid"]
                            sp_obj_sub["sn"] = sp_obj["sn"]
                            sp_obj_sub["fc"] = sp_obj["fc"]
                            sp_obj_sub["sw"] = sp_obj["sw"]
                            sp_obj_sub["gf"] = sp_obj["gf"]
                            try:
                                sp_obj_sub["pj"] = sp_obj["pj"]
                            except:
                                continue
                            sp_obj_sub["idx1"] = tenXidx.replace(",", "")
                            sp_obj_sub["idx1ref"] = tenXidx.replace(",", "")
                            data.append(sp_obj_sub)
                    else:
                        sp_obj["idx1"] = idxs[0].replace(",", "").upper()
                        sp_obj["idx1ref"] = idxs[0].replace(",", "").upper()
                        if len(idxs) == 2:
                            dualindex = True
                            if pj_type != "by user":
                                compl = {"A": "T", "C": "G", "G": "C", "T": "A"}
                                sp_obj["idx2"] = "".join(
                                    reversed(
                                        [
                                            compl.get(b, b)
                                            for b in idxs[1].replace(",", "").upper()
                                        ]
                                    )
                                )
                                sp_obj["idx2ref"] = "".join(
                                    reversed(
                                        [
                                            compl.get(b, b)
                                            for b in idxs[1].replace(",", "").upper()
                                        ]
                                    )
                                )
                            else:
                                sp_obj["idx2"] = idxs[1].replace(",", "").upper()
                                sp_obj["idx2ref"] = idxs[1].replace(",", "").upper()
                        else:
                            header_ar.remove("index2")
                            header_ar.remove("I5_Index_ID")
                        data.append(sp_obj)
    header = "[Data]\n{}\n".format(",".join(header_ar))
    str_data = ""
    for line in data:
        if noindex:
            l_data = [
                line["sn"],
                line["sn"],
                line["fc"],
                line["sw"],
                line["pj"],
                pro.udf["Description"].replace(".", "_"),
                line["gf"],
            ]
        elif dualindex:
            l_data = [
                line["sn"],
                line["sn"],
                line["fc"],
                line["sw"],
                line["pj"],
                line["idx1"],
                line["idx1ref"],
                line["idx2"],
                line["idx2ref"],
                pro.udf["Description"].replace(".", "_"),
                line["gf"],
            ]
        else:
            l_data = [
                line["sn"],
                line["sn"],
                line["fc"],
                line["sw"],
                line["pj"],
                line["idx1"],
                line["idx1ref"],
                pro.udf["Description"].replace(".", "_"),
                line["gf"],
            ]
        str_data = str_data + ",".join(l_data) + "\n"

    content = "{}{}".format(header, str_data)
    df = pd.read_csv(StringIO(content), skiprows=1)
    df = df.sort_values(["Sample_ID"])
    content = df.to_csv(index=False)
    content = "[Data]\n{}\n".format(content)

    return (content, data)


def gen_Nextseq_lane_data(pro):
    data = []
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
                    if NGISAMPLE_PAT.findall(sample.name):
                        sp_obj["sid"] = "Sample_{}".format(sample.name).replace(",", "")
                        sp_obj["sn"] = sample.name.replace(",", "")
                        sp_obj["pj"] = sample.project.name.replace(".", "__").replace(
                            ",", ""
                        )
                        sp_obj["ref"] = sample.project.udf.get(
                            "Reference genome", ""
                        ).replace(",", "")
                        seq_setup = sample.project.udf.get("Sequencing setup", "")
                        if SEQSETUP_PAT.findall(seq_setup):
                            sp_obj["rc"] = "{}-{}".format(
                                seq_setup.split("-")[0], seq_setup.split("-")[3]
                            )
                        else:
                            sp_obj["rc"] = "0-0"
                    else:
                        sp_obj["sid"] = (
                            "Sample_{}".format(sample.name)
                            .replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["sn"] = (
                            sample.name.replace("(", "")
                            .replace(")", "")
                            .replace(".", "")
                            .replace(" ", "_")
                        )
                        sp_obj["pj"] = "Control"
                        sp_obj["ref"] = "Control"
                        sp_obj["rc"] = "0-0"
                    sp_obj["ct"] = "N"
                    sp_obj["op"] = pro.technician.name.replace(" ", "_").replace(
                        ",", ""
                    )
                    sp_obj["fc"] = out.location[0].name.replace(",", "")
                    sp_obj["sw"] = out.location[1].replace(",", "")
                    sp_obj["idx1"] = idxs[0].replace(",", "")
                    if idxs[1]:
                        sp_obj["idx2"] = idxs[1].replace(",", "").upper()
                    else:
                        sp_obj["idx2"] = ""
                    data.append(sp_obj)
    header = "{}\n".format(",".join(header_ar))
    str_data = ""
    for line in sorted(data, key=lambda x: x["lane"]):
        l_data = [
            line["fc"],
            line["lane"],
            line["sn"],
            line["sn"],
            line["ref"],
            line["idx1"],
            line["idx2"],
            line["pj"],
            line["ct"],
            line["rc"],
            line["op"],
            line["pj"],
        ]
        str_data = str_data + ",".join(l_data) + "\n"

    content = "{}{}".format(header, str_data)
    df = pd.read_csv(StringIO(content))
    df = df.sort_values(["Lane", "Sample_ID"])
    content = df.to_csv(index=False)

    return (content, data)


def gen_MinION_QC_data(pro):
    keep_idx_flag = True if pro.type.name == "MinION QC" else False
    data = []
    for out in pro.all_outputs():
        if NGISAMPLE_PAT.findall(out.name):
            nanopore_barcode_name = (
                out.udf["Nanopore Barcode"].split("_")[0]
                if out.udf["Nanopore Barcode"] != "None"
                else ""
            )
            nanopore_barcode_seq = (
                out.udf["Nanopore Barcode"].split("_")[1]
                if out.udf["Nanopore Barcode"] != "None"
                else ""
            )
            sample_name = out.name
            idxs = out.reagent_labels[0]

            sp_obj = {}
            sp_obj["sn"] = sample_name
            sp_obj["npbs"] = nanopore_barcode_seq

            # Case of 10X indexes
            if TENX_SINGLE_PAT.findall(idxs):
                for tenXidx in Chromium_10X_indexes[TENX_SINGLE_PAT.findall(idxs)[0]]:
                    tenXidx_no = (
                        Chromium_10X_indexes[TENX_SINGLE_PAT.findall(idxs)[0]].index(
                            tenXidx
                        )
                        + 1
                    )
                    sp_obj_sub = {}
                    sp_obj_sub["sn"] = sp_obj["sn"] + "_" + str(tenXidx_no)
                    sp_obj_sub["npbs"] = sp_obj["npbs"]
                    sp_obj_sub["idxt"] = "truseq"
                    sp_obj_sub["idx"] = tenXidx.replace(",", "")
                    data.append(sp_obj_sub)
            # Case of 10X dual indexes
            elif TENX_DUAL_PAT.findall(idxs):
                sp_obj["idxt"] = "truseq_dual"
                sp_obj["idx"] = (
                    Chromium_10X_indexes[TENX_DUAL_PAT.findall(idxs)[0]][0]
                    + "-"
                    + Chromium_10X_indexes[TENX_DUAL_PAT.findall(idxs)[0]][1]
                )
                data.append(sp_obj)
            # Case of NoIndex
            elif idxs == "NoIndex" or idxs == "" or not idxs:
                sp_obj["idxt"] = "truseq"
                sp_obj["idx"] = ""
                data.append(sp_obj)
            # Case of index sequences between brackets
            elif re.findall("\((.*?)\)", idxs):
                idxs = re.findall("\((.*?)\)", idxs)[0]
                if "-" not in idxs:
                    sp_obj["idxt"] = "truseq"
                    sp_obj["idx"] = idxs
                    data.append(sp_obj)
                else:
                    sp_obj["idxt"] = "truseq_dual"
                    sp_obj["idx"] = idxs
                    data.append(sp_obj)
            # Case of single index
            elif "-" not in idxs:
                sp_obj["idxt"] = "truseq"
                sp_obj["idx"] = idxs
                data.append(sp_obj)
            # Case of dual index
            else:
                sp_obj["idxt"] = "truseq_dual"
                sp_obj["idx"] = idxs
                data.append(sp_obj)
    str_data = ""
    for line in sorted(data, key=lambda x: x["sn"]):
        if keep_idx_flag:
            l_data = [line["sn"], line["npbs"], line["idxt"], line["idx"]]
        else:
            l_data = [line["sn"], line["npbs"], "", ""]
        str_data = str_data + ",".join(l_data) + "\n"

    return str_data


def find_barcode(sample_idxs, sample, process):
    # print "trying to find {} barcode in {}".format(sample.name, process.name)
    for art in process.all_inputs():
        if sample in art.samples:
            if len(art.samples) == 1 and art.reagent_labels:
                reagent_label_name = art.reagent_labels[0].upper().replace(" ", "")
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


def test():
    log = []
    d = [
        {"lane": 1, "idx1": "ATTT", "idx2": ""},
        {"lane": 1, "idx1": "ATCTATCG", "idx2": ""},
        {"lane": 1, "idx1": "ATCG", "idx2": "ATCG"},
    ]
    check_index_distance(d, log)
    print(log)


def main(lims, args):
    log = []
    thisyear = datetime.now().year
    content = None
    if args.mytest:
        test()
    else:
        process = Process(lims, id=args.pid)

        if process.type.name == "Load to Flowcell (NovaSeq 6000 v2.0)":
            (content, obj) = gen_Novaseq_lane_data(process)
            check_index_distance(obj, log)
            if os.path.exists("/srv/mfs/samplesheets/novaseq/{}".format(thisyear)):
                try:
                    with open(
                        "/srv/mfs/samplesheets/novaseq/{}/{}.csv".format(
                            thisyear, obj[0]["fc"]
                        ),
                        "w",
                    ) as sf:
                        sf.write(content)
                except Exception as e:
                    log.append(str(e))

        elif process.type.name == "Denature, Dilute and Load Sample (MiSeq) 4.0":
            header = gen_Miseq_header(process)
            reads = gen_Miseq_reads(process)
            settings = gen_Miseq_settings(process)
            (data, obj) = gen_Miseq_data(process)
            check_index_distance(obj, log)
            content = "{}{}{}{}".format(header, reads, settings, data)

        elif process.type.name == "Load to Flowcell (NextSeq v1.0)":
            (content, obj) = gen_Nextseq_lane_data(process)
            check_index_distance(obj, log)
            nextseq_fc = (
                process.udf["Flowcell Series Number"]
                if process.udf["Flowcell Series Number"]
                else obj[0]["fc"]
            )
            if os.path.exists("/srv/mfs/samplesheets/nextseq/{}".format(thisyear)):
                try:
                    with open(
                        "/srv/mfs/samplesheets/nextseq/{}/{}.csv".format(
                            thisyear, nextseq_fc
                        ),
                        "w",
                    ) as sf:
                        sf.write(content)
                except Exception as e:
                    log.append(str(e))

        elif process.type.name in [
            "MinION QC",
            "Load Sample and Sequencing (MinION) 1.0",
        ]:
            content = gen_MinION_QC_data(process)
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
            if os.path.exists("/srv/mfs/samplesheets/nanopore/{}".format(thisyear)):
                try:
                    with open(
                        "/srv/mfs/samplesheets/nanopore/{}/{}.csv".format(
                            thisyear, fc_name
                        ),
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
                    log_id = out.id
                elif out.type == "Analyte":
                    if process.type.name == "Load to Flowcell (NextSeq v1.0)":
                        fc_name = (
                            process.udf["Flowcell Series Number"]
                            if process.udf["Flowcell Series Number"]
                            else out.location[0].name
                        )
                    else:
                        fc_name = out.location[0].name
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

            with open("{}.csv".format(fc_name), "w", 0o664) as f:
                f.write(content)
            os.chmod("{}.csv".format(fc_name), 0o664)
            for f in ss_art.files:
                lims.request_session.delete(f.uri)
            lims.upload_new_file(ss_art, "{}.csv".format(fc_name))
            if log:
                with open("{}_{}_Error.log".format(log_id, fc_name), "w") as f:
                    f.write("\n".join(log))

                sys.stderr.write("Errors were met, check the log.")
                sys.exit(1)

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
