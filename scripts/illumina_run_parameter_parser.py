#!/usr/bin/env python

import os
import sys
import glob
import math

from argparse import ArgumentParser
from datetime import datetime
from genologics.lims import Lims
from genologics.entities import Process
from genologics.config import BASEURI, USERNAME, PASSWORD
from flowcell_parser.classes import RunParser, RunParametersParser
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary


DESC = """EPP for parsing run paramters for Illumina MiSeq, NextSeq and NovaSeq runs
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""


def fetch_fc(process):
    fc_id = ''
    if process.parent_processes()[0].type.name == 'Load to Flowcell (NextSeq v1.0)':
        if process.parent_processes()[0].udf['Flowcell Series Number']:
            fc_id = process.parent_processes()[0].udf['Flowcell Series Number'].upper()
        else:
            sys.stderr.write("Flowcell Series Number is empty in the associated Load to Flowcell (NextSeq v1.0) step.")
            sys.exit(2)
    elif process.parent_processes()[0].type.name == 'Denature, Dilute and Load Sample (MiSeq) 4.0':
        if process.parent_processes()[0].udf['Flowcell ID']:
            fc_id = process.parent_processes()[0].udf['Flowcell ID'].upper()
        else:
            sys.stderr.write("Flowcell ID is empty in the associated Denature, Dilute and Load Sample (MiSeq) 4.0 step.")
            sys.exit(2)
    elif process.parent_processes()[0].type.name == 'Load to Flowcell (NovaSeq 6000 v2.0)':
        fc_id = process.parent_processes()[0].output_containers()[0].name
    elif "Load to Flowcell (NovaSeqXPlus)" in process.parent_processes()[0].type.name:
        fc_id = process.parent_processes()[0].output_containers()[0].name
    else:
        sys.stderr.write("No associated parent step can be found.")
        sys.exit(2)
    return fc_id


def fetch_rundir(fc_id, run_type):
    run_dir = ''
    if run_type == 'nextseq':
        data_dir = 'NextSeq_data'
    elif run_type == 'miseq':
        data_dir = 'miseq_data'
    elif run_type == 'novaseq':
        data_dir = 'NovaSeq_data'
    elif run_type == "NovaSeqXPlus":
        data_dir = "NovaSeqXPlus_data"

    metadata_dir = "ngi-nas-ns"
    run_dir_path = os.path.join(
        os.sep, "srv", metadata_dir, data_dir, "*{}".format(fc_id)
    )

    if len(glob.glob(run_dir_path)) == 1:
        run_dir = glob.glob(run_dir_path)[0]
    elif len(glob.glob(run_dir_path)) == 0:
        sys.stderr.write("No run dir can be found for FC {}".format(fc_id))
        sys.exit(2)
    else:
        sys.stderr.write("Multiple run dirs found for FC {}".format(fc_id))
        sys.exit(2)
    return run_dir


def parse_run(run_dir):
    runParserObj = RunParser(run_dir)
    if os.path.exists("{}/RunParameters.xml".format(run_dir)):
        RunParametersParserObj = RunParametersParser("{}/RunParameters.xml".format(run_dir))
    elif os.path.exists("{}/runParameters.xml".format(run_dir)):
        RunParametersParserObj = RunParametersParser("{}/runParameters.xml".format(run_dir))
    else:
        sys.stderr.write("No RunParameters.xml found in path {}".format(run_dir))
        sys.exit(2)
    return runParserObj, RunParametersParserObj


def attach_xml(process, run_dir):
    for outart in process.all_outputs():
        if outart.type == 'ResultFile' and outart.name == 'Run Info':
            try:
                lims.upload_new_file(outart, "{}/RunInfo.xml".format(run_dir))
            except IOError:
                try:
                    lims.upload_new_file(outart, "{}/runInfo.xml".format(run_dir))
                except IOError:
                    sys.stderr.write("No RunInfo.xml found")
                    sys.exit(2)
        elif outart.type == 'ResultFile' and outart.name == 'Run Parameters':
            try:
                lims.upload_new_file(outart, "{}/RunParameters.xml".format(run_dir))
            except IOError:
                try:
                    lims.upload_new_file(outart, "{}/runParameters.xml".format(run_dir))
                except IOError:
                    sys.stderr.write("No RunParameters.xml found")
                    sys.exit(2)


def parse_illumina_interop(run_dir):
    # Read interop
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    try:
        run_metrics.read(run_dir, valid_to_load)
    except Exception:
        sys.stderr.write("Cannot parse information in InterOp")
        sys.exit(2)
    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)
    lanes = summary.lane_count()
    reads = summary.size()
    # Parse the interop stats lane by lane for non-index reads
    run_stats_summary = dict()
    for lane in range(lanes):
        lane_nbr = getattr(summary.at(0).at(lane), "lane")()
        for read in range(reads):
            if not summary.at(read).read().is_index():
                stats = dict()
                stats.update({"density"               :  getattr(summary.at(read).at(lane), "density")().mean()/1000})
                stats.update({"error_rate"            :  getattr(summary.at(read).at(lane), "error_rate")().mean()})
                stats.update({"first_cycle_intensity" :  getattr(summary.at(read).at(lane), "first_cycle_intensity")().mean()})
                stats.update({"percent_aligned"       :  getattr(summary.at(read).at(lane), "percent_aligned")().mean()})
                stats.update({"percent_gt_q30"        :  getattr(summary.at(read).at(lane), "percent_gt_q30")()})
                stats.update({"percent_pf"            :  getattr(summary.at(read).at(lane), "percent_pf")().mean()})
                stats.update({"phasing"               :  getattr(summary.at(read).at(lane), "phasing")().mean()})
                stats.update({"prephasing"            :  getattr(summary.at(read).at(lane), "prephasing")().mean()})
                stats.update({"reads_pf"              :  getattr(summary.at(read).at(lane), "reads_pf")()})
                stats.update({"yield_g"               :  getattr(summary.at(read).at(lane), "yield_g")()})
                if lane_nbr not in list(run_stats_summary.keys()):
                    run_stats_summary.update({lane_nbr: {read: stats}})
                else:
                    run_stats_summary[lane_nbr].update({read : stats})
    return run_stats_summary


def set_run_stats_in_lims(process, run_stats_summary):
    for art in process.all_outputs():
        if 'Lane' in art.name:
            lane_nbr = int(art.name.split(' ')[1])
            read = 1
            for i in list(run_stats_summary[lane_nbr].keys()):
                lane_stats_for_read = run_stats_summary[lane_nbr][i]
                if not math.isnan(lane_stats_for_read['density']):
                    art.udf["Cluster Density (K/mm^2) R{}".format(read)]  =  lane_stats_for_read['density']
                if not math.isnan(lane_stats_for_read['error_rate']):
                    art.udf["% Error Rate R{}".format(read)]              =  lane_stats_for_read['error_rate']
                if not math.isnan(lane_stats_for_read['first_cycle_intensity']):
                    art.udf["Intensity Cycle 1 R{}".format(read)]         =  lane_stats_for_read['first_cycle_intensity']
                if not math.isnan(lane_stats_for_read['percent_aligned']):
                    art.udf["% Aligned R{}".format(read)]                 =  lane_stats_for_read['percent_aligned']
                if not math.isnan(lane_stats_for_read['percent_gt_q30']):
                    art.udf["% Bases >=Q30 R{}".format(read)]             =  lane_stats_for_read['percent_gt_q30']
                if not math.isnan(lane_stats_for_read['percent_pf']):
                    art.udf["%PF R{}".format(read)]                       =  lane_stats_for_read['percent_pf']
                if not math.isnan(lane_stats_for_read['phasing']):
                    art.udf["% Phasing R{}".format(read)]                 =  lane_stats_for_read['phasing']
                if not math.isnan(lane_stats_for_read['prephasing']):
                    art.udf["% Prephasing R{}".format(read)]              =  lane_stats_for_read['prephasing']
                if not math.isnan(lane_stats_for_read['reads_pf']):
                    art.udf["Reads PF (M) R{}".format(read)]              =  lane_stats_for_read['reads_pf']/1000000
                if not math.isnan(lane_stats_for_read['yield_g']):
                    art.udf["Yield PF (Gb) R{}".format(read)]             =  lane_stats_for_read['yield_g']
                read += 1
            art.put()
    process.put()


def set_run_stats_in_lims_miseq(process, run_stats_summary):
    art = process.input_output_maps[0][0]['uri']
    lane_nbr = 1
    read = 1
    for i in list(run_stats_summary[lane_nbr].keys()):
        lane_stats_for_read = run_stats_summary[lane_nbr][i]
        if not math.isnan(lane_stats_for_read['density']):
            art.udf["Cluster Density (K/mm^2) R{}".format(read)]  =  lane_stats_for_read['density']
        if not math.isnan(lane_stats_for_read['error_rate']):
            art.udf["% Error Rate R{}".format(read)]              =  lane_stats_for_read['error_rate']
        if not math.isnan(lane_stats_for_read['first_cycle_intensity']):
            art.udf["Intensity Cycle 1 R{}".format(read)]         =  lane_stats_for_read['first_cycle_intensity']
        if not math.isnan(lane_stats_for_read['percent_aligned']):
            art.udf["% Aligned R{}".format(read)]                 =  lane_stats_for_read['percent_aligned']
        if not math.isnan(lane_stats_for_read['percent_gt_q30']):
            art.udf["% Bases >=Q30 R{}".format(read)]             =  lane_stats_for_read['percent_gt_q30']
        if not math.isnan(lane_stats_for_read['percent_pf']):
            art.udf["%PF R{}".format(read)]                       =  lane_stats_for_read['percent_pf']
        if not math.isnan(lane_stats_for_read['phasing']):
            art.udf["% Phasing R{}".format(read)]                 =  lane_stats_for_read['phasing']
        if not math.isnan(lane_stats_for_read['prephasing']):
            art.udf["% Prephasing R{}".format(read)]              =  lane_stats_for_read['prephasing']
        if not math.isnan(lane_stats_for_read['reads_pf']):
            art.udf["Clusters PF R{}".format(read)]              =  lane_stats_for_read['reads_pf']
        if not math.isnan(lane_stats_for_read['yield_g']):
            art.udf["Yield PF (Gb) R{}".format(read)]             =  lane_stats_for_read['yield_g']
        read += 1
    art.put()
    process.put()


def lims_for_nextseq(process, run_dir):
    # Parse run
    runParserObj, RunParametersParserObj = parse_run(run_dir)
    # Attach RunInfo.xml and RunParamters.xml
    attach_xml(process, run_dir)
    # Set values for LIMS UDFs
    runParameters = RunParametersParserObj.data['RunParameters']
    process.udf['Finish Date'] = datetime.strptime(runParameters['RunEndTime'][:10],"%Y-%m-%d").date() if 'RunEndTime' in list(runParameters.keys()) and runParameters['RunEndTime'] != '' else datetime.now().date()
    process.udf['Run Type'] = "NextSeq 2000 {}".format(runParameters['FlowCellMode'].split(' ')[2])
    process.udf['Chemistry'] = "NextSeq 2000 {}".format(runParameters['FlowCellMode'].split(' ')[2])
    planned_cycles = sum(list(map(int, list(runParameters['PlannedCycles'].values()))))
    completed_cycles = sum(list(map(int, list(runParameters['CompletedCycles'].values()))))
    process.udf['Status'] = "Cycle {} of {}".format(completed_cycles, planned_cycles)
    process.udf['Flow Cell ID'] = runParameters['FlowCellSerialNumber']
    process.udf['Experiment Name'] = runParameters['FlowCellSerialNumber']
    process.udf['Read 1 Cycles'] = int(runParameters['PlannedCycles']['Read1'])
    process.udf['Index 1 Read Cycles'] = int(runParameters['PlannedCycles']['Index1'])
    process.udf['Index 2 Read Cycles'] = int(runParameters['PlannedCycles']['Index2'])
    process.udf['Read 2 Cycles'] = int(runParameters['PlannedCycles']['Read2'])
    process.udf['Run ID'] = runParserObj.runinfo.data['Id']
    process.udf['Reagent Cartridge ID'] = runParameters['CartridgeSerialNumber']
    # Put in LIMS
    process.put()
    # Set run stats parsed from InterOp
    run_stats_summary = parse_illumina_interop(run_dir)
    set_run_stats_in_lims(process, run_stats_summary)


def lims_for_miseq(process, run_dir):
    # Parse run
    runParserObj, RunParametersParserObj = parse_run(run_dir)
    # Attach RunInfo.xml and RunParamters.xml
    attach_xml(process, run_dir)
    # Set values for LIMS UDFs
    runParameters = RunParametersParserObj.data['RunParameters']
    process.udf['Finish Date'] = datetime.now().date()
    if runParameters['Setup']['SupportMultipleSurfacesInUI'] == 'true' and runParameters['Setup']['NumTilesPerSwath'] == '19':
        process.udf['Run Type'] = 'Version3'
    elif runParameters['Setup']['SupportMultipleSurfacesInUI'] == 'true' and runParameters['Setup']['NumTilesPerSwath'] == '14':
        process.udf['Run Type'] = 'Version2'
    elif runParameters['Setup']['SupportMultipleSurfacesInUI'] == 'false' and runParameters['Setup']['NumTilesPerSwath'] == '2':
        process.udf['Run Type'] = 'Version2Nano'
    else:
        process.udf['Run Type'] = 'null'
    # Runs with single read return a dict object
    if isinstance(runParameters['Reads']['RunInfoRead'], list):
        total_cycles = sum(list(map(int, [read['NumCycles'] for read in runParameters['Reads']['RunInfoRead']])))
        non_index_read_idx = [read['Number'] for read in runParameters['Reads']['RunInfoRead'] if read['IsIndexedRead'] == 'N']
        index_read_idx = [read['Number'] for read in runParameters['Reads']['RunInfoRead'] if read['IsIndexedRead'] == 'Y']
        process.udf['Read 1 Cycles'] = int(list([read for read in runParameters['Reads']['RunInfoRead'] if read['Number'] == str(min(list(map(int,non_index_read_idx))))])[0]['NumCycles'])
    elif isinstance(runParameters['Reads']['RunInfoRead'], dict):
        total_cycles = int(runParameters['Reads']['RunInfoRead']['NumCycles'])
        non_index_read_idx = [runParameters['Reads']['RunInfoRead']['Number']]
        index_read_idx = []
        process.udf['Read 1 Cycles'] = int(runParameters['Reads']['RunInfoRead']['NumCycles'])

    process.udf['Status'] = "Cycle {} of {}".format(total_cycles, total_cycles)
    process.udf['Flow Cell ID'] = runParameters['FlowcellRFIDTag']['SerialNumber']
    process.udf['Flow Cell Version'] = runParameters['FlowcellRFIDTag']['PartNumber']
    process.udf['Experiment Name'] = process.all_inputs()[0].name

    if len(non_index_read_idx) == 2:
        process.udf['Read 2 Cycles'] = int(list([read for read in runParameters['Reads']['RunInfoRead'] if read['Number'] == str(max(list(map(int,non_index_read_idx))))])[0]['NumCycles'])

    if len(index_read_idx) > 0:
        process.udf['Index 1 Read Cycles'] = int(list([read for read in runParameters['Reads']['RunInfoRead'] if read['Number'] == str(min(list(map(int,index_read_idx))))])[0]['NumCycles'])
    if len(index_read_idx) == 2:
        process.udf['Index 2 Read Cycles'] = int(list([read for read in runParameters['Reads']['RunInfoRead'] if read['Number'] == str(max(list(map(int,index_read_idx))))])[0]['NumCycles'])


    process.udf['Run ID'] = runParameters['RunID']
    process.udf['Output Folder'] = runParameters['OutputFolder'].replace(runParameters['RunID'], '')
    process.udf['Reagent Cartridge ID'] = runParameters['ReagentKitRFIDTag']['SerialNumber']
    process.udf['Reagent Cartridge Part #'] = runParameters['ReagentKitRFIDTag']['PartNumber']
    process.udf['PR2 Bottle ID'] = runParameters['PR2BottleRFIDTag']['SerialNumber']
    process.udf['Chemistry'] = runParameters['Chemistry']
    process.udf['Workflow'] = runParameters['Workflow']['Analysis'] if runParameters.get('Workflow') else runParameters['ModuleName']
    # Put in LIMS
    process.put()
    # Set run stats parsed from InterOp
    run_stats_summary = parse_illumina_interop(run_dir)
    set_run_stats_in_lims_miseq(process, run_stats_summary)


def lims_for_novaseq(process, run_dir):
    # Parse run
    runParserObj, RunParametersParserObj = parse_run(run_dir)
    # Set values for LIMS UDFs
    runParameters = RunParametersParserObj.data['RunParameters']
    process.udf['Flow Cell ID'] = runParameters['RfidsInfo']['FlowCellSerialBarcode']
    process.udf['Flow Cell Part Number'] = runParameters['RfidsInfo']['FlowCellPartNumber']
    process.udf['Flow Cell Lot Number'] = runParameters['RfidsInfo']['FlowCellLotNumber']
    process.udf['Flow Cell Expiration Date'] = datetime.strptime(runParameters['RfidsInfo']['FlowCellExpirationdate'], "%m/%d/%Y %H:%M:%S").date()
    process.udf['Flow Cell Mode'] = runParameters['RfidsInfo']['FlowCellMode']
    process.udf['Run ID'] = runParameters['RunId']
    process.udf['Read 1 Cycles'] = int(runParameters['Read1NumberOfCycles'])
    process.udf['Read 2 Cycles'] = int(runParameters['Read2NumberOfCycles'])
    process.udf['Index Read 1'] = int(runParameters['IndexRead1NumberOfCycles'])
    process.udf['Index Read 2'] = int(runParameters['IndexRead2NumberOfCycles'])
    process.udf['PE Serial Barcode'] = runParameters['RfidsInfo']['ClusterSerialBarcode']
    process.udf['PE Part Number'] = runParameters['RfidsInfo']['ClusterPartNumber']
    process.udf['PE Lot Number'] = runParameters['RfidsInfo']['ClusterLotNumber']
    process.udf['PE Expiration Date'] = datetime.strptime(runParameters['RfidsInfo']['ClusterExpirationdate'], "%m/%d/%Y %H:%M:%S").date()
    process.udf['PE Cycle Kit'] = runParameters['RfidsInfo']['ClusterCycleKit']
    process.udf['SBS Serial Barcode'] = runParameters['RfidsInfo']['SbsSerialBarcode']
    process.udf['SBS Part Number'] = runParameters['RfidsInfo']['SbsPartNumber']
    process.udf['SBS Lot Number'] = runParameters['RfidsInfo']['SbsLotNumber']
    process.udf['SBS Expiration Date'] = datetime.strptime(runParameters['RfidsInfo']['SbsExpirationdate'], "%m/%d/%Y %H:%M:%S").date()
    process.udf['SBS Cycle Kit'] = runParameters['RfidsInfo']['SbsCycleKit']
    process.udf['Buffer Serial Barcode'] = runParameters['RfidsInfo']['BufferSerialBarcode']
    process.udf['Buffer Part Number'] = runParameters['RfidsInfo']['BufferPartNumber']
    process.udf['Buffer Lot Number'] = runParameters['RfidsInfo']['BufferLotNumber']
    process.udf['Buffer Expiration Date'] = datetime.strptime(runParameters['RfidsInfo']['BufferExpirationdate'], "%m/%d/%Y %H:%M:%S").date()
    process.udf['Output Folder'] = runParameters['OutputRunFolder']
    process.udf['Loading Workflow Type'] = runParameters['WorkflowType']
    # Put in LIMS
    process.put()
    # Set run stats parsed from InterOp
    run_stats_summary = parse_illumina_interop(run_dir)
    set_run_stats_in_lims(process, run_stats_summary)


def lims_for_NovaSeqXPlus(process, run_dir):
    # Parse run
    runParserObj, RunParametersParserObj = parse_run(run_dir)

    # Subset parsed data
    runParameters = RunParametersParserObj.data["RunParameters"]
    consumables = runParameters["ConsumableInfo"]["ConsumableInfo"]
    reads = runParameters["PlannedReads"]["Read"]

    # Set values for LIMS UDFs
    process.udf["Run ID"] = runParameters["RunId"]
    process.udf["Output Folder"] = runParameters["OutputFolder"]

    for read in reads:
        if read["ReadName"] == "Read1":
            process.udf["Read 1 Cycles"] = int(read["Cycles"])
        elif read["ReadName"] == "Index1":
            process.udf["Index Read 1"] = int(read["Cycles"])
        elif read["ReadName"] == "Index2":
            process.udf["Index Read 2"] = int(read["Cycles"])
        elif read["ReadName"] == "Read2":
            process.udf["Read 2 Cycles"] = int(read["Cycles"])

    for consumable in consumables:
        if consumable["Type"] == "FlowCell":
            process.udf["Flow Cell Mode"] = consumable["Name"] if consumable.get("Name", "") else consumable["Mode"]

            process.udf["Flow Cell ID"] = consumable["SerialNumber"]
            process.udf["Flow Cell Part Number"] = consumable["PartNumber"]
            process.udf["Flow Cell Lot Number"] = consumable["LotNumber"]
            process.udf["Flow Cell Expiration Date"] = datetime.strptime(
                consumable["ExpirationDate"][0:10], "%Y-%m-%d"
            ).date()

        elif consumable["Type"] == "Reagent":
            process.udf["Reagent Serial Barcode"] = consumable["SerialNumber"]
            process.udf["Reagent Part Number"] = consumable["PartNumber"]
            process.udf["Reagent Lot Number"] = consumable["LotNumber"]
            process.udf["Reagent Expiration Date"] = datetime.strptime(
                consumable["ExpirationDate"][0:10], "%Y-%m-%d"
            ).date()

        elif consumable["Type"] == "Buffer":
            process.udf["Buffer Serial Barcode"] = consumable["SerialNumber"]
            process.udf["Buffer Part Number"] = consumable["PartNumber"]
            process.udf["Buffer Lot Number"] = consumable["LotNumber"]
            process.udf["Buffer Expiration Date"] = datetime.strptime(
                consumable["ExpirationDate"][0:10], "%Y-%m-%d"
            ).date()

        elif consumable["Type"] == "SampleTube":
            process.udf["SampleTube Serial Barcode"] = consumable["SerialNumber"]
            process.udf["SampleTube Part Number"] = consumable["PartNumber"]
            process.udf["SampleTube Lot Number"] = consumable["LotNumber"]
            process.udf["SampleTube Expiration Date"] = datetime.strptime(
                consumable["ExpirationDate"][0:10], "%Y-%m-%d"
            ).date()

        elif consumable["Type"] == "Lyo":
            process.udf["Lyo Serial Barcode"] = consumable["SerialNumber"]
            process.udf["Lyo Part Number"] = consumable["PartNumber"]
            process.udf["Lyo Lot Number"] = consumable["LotNumber"]
            process.udf["Lyo Expiration Date"] = datetime.strptime(
                consumable["ExpirationDate"][0:10], "%Y-%m-%d"
            ).date()

    # Put in LIMS
    process.put()

    # Set run stats parsed from InterOp to measurement UDFs
    run_stats_summary = parse_illumina_interop(run_dir)
    set_run_stats_in_lims(process, run_stats_summary)


def main(lims, args):

    process = Process(lims, id=args.pid)

    if process.type.name == "Illumina Sequencing (NextSeq) v1.0":
        run_type = "nextseq"
    elif process.type.name == "MiSeq Run (MiSeq) 4.0":
        run_type = "miseq"
    elif "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)" in process.type.name:
        run_type = "novaseq"
    elif "NovaSeqXPlus Run" in process.type.name:
        run_type = "NovaSeqXPlus"

    # Fetch FC ID
    fc_id = fetch_fc(process)

    # Fetch run dir
    run_dir = fetch_rundir(fc_id, run_type)

    # Fill info in LIMS
    if run_type == 'nextseq':
        lims_for_nextseq(process, run_dir)
    elif run_type == 'miseq':
        lims_for_miseq(process, run_dir)
    elif run_type == 'novaseq':
        lims_for_novaseq(process, run_dir)
    elif run_type == "NovaSeqXPlus":
        lims_for_NovaSeqXPlus(process, run_dir)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument('--pid',
                        help='Lims id for current Process')
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
