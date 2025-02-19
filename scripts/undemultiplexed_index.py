#!/usr/bin/env python

DESC = """This EPP script reads demultiplex end undemultiplexed yields from file
system and then does the following

1)  Sets the output artifact udfs based demultiplex resultfile.

2)  Generated qc tresholds based on sequencing run setings. These can also be
    set by user.

3)  Loggs the tresholds to a qc log file

4)  Sets the output artifact qc-flaggs based on the threshold

5)  Warns if anny unexpected index has high yield ore if a lane has high yield
    of un expected indexes

6)  Loads a result file with demultiplex end undemultiplexed yields. This should
    be checked if warnings are given.

Reads from:
    --files--
    Demultiplex_Stats.htm                           in ngi-nas-ns file system
    Undemultiplexed_stats.metrics                   in ngi-nas-ns file system

Logging:
    The script outputs a regular log file with regular execution information.

Written by Maya Brandi
"""
import csv
import glob
import logging
import sys
from argparse import ArgumentParser

import numpy as np
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from scilifelab_parsers.qc.qc import FlowcellRunMetricsParser

from scilifelab_epps.epp import EppLogger, set_field

# from qc_parsers import FlowcellRunMetricsParser


class RunQC:
    def __init__(self, process):
        ## Processes, artifacts and udfs
        self.process = process
        self.input_pools = process.all_inputs()
        self.seq_run = None
        self.run_udfs = {}
        self.user_def_tresh = dict(list(self.process.udf.items()))

        ##  Stuff for logging
        self.qc_log_file = None
        self.nr_lane_samps_updat = 0
        self.nr_lane_samps_tot = 0
        self.abstract = []
        self.QC_fail = []
        self.high_index_yield = []
        self.high_lane_yield = []
        self.html_file_error = False  ## Never used????

        ##  Demultiplexing result files
        self.file_path = None
        self.dem_stat = None
        self.undem_stat = None

        ##  Other variables
        self.single = True
        self.read_length = None
        self.Q30_treshold = None

    def make_qc_log_file(self, qc_log_file):
        """File to logg qc tresholds."""
        self.qc_log_file = open(qc_log_file, "a")

    def get_run_info(self):
        cont_name = self._get_container()
        self._get_run(cont_name)
        self.run_udfs = dict(list(self.seq_run.udf.items()))
        self._get_run_id()
        self._get_cycles()
        self._get_file_path(cont_name)
        self._get_demultiplex_files()
        self._get_threshold_Q30()

    def _get_container(self):
        """Need container name to fetch parrent sequencing process"""
        try:
            return self.process.all_inputs()[0].location[0].name
        except:
            sys.exit("Could not find container name.")

    def _get_run(self, cont_name):
        """Getting parrent sequencing process and process type"""
        miseq = lims.get_processes(
            udf={"Reagent Cartridge ID": cont_name}, type="MiSeq Run (MiSeq) 4.0"
        )
        hiseq = lims.get_processes(
            udf={"Flow Cell ID": cont_name},
            type="Illumina Sequencing (Illumina SBS) 4.0",
        )
        hiseq_X10 = lims.get_processes(
            udf={"Flow Cell ID": cont_name}, type="Illumina Sequencing (HiSeq X) 1.0"
        )
        if miseq:
            self.seq_run = miseq[0]
            self.run_type = "MiSeq"
        elif hiseq:
            self.seq_run = hiseq[0]
            try:
                self.run_type = self.seq_run.udf["Flow Cell Version"]
            except:
                sys.exit("Missing field 'Flow Cell Version' in sequencing process")
        elif hiseq_X10:
            self.seq_run = hiseq_X10[0]
            self.run_type = "HiSeqX10"
        else:
            sys.exit("run not found")

    def _get_run_id(self):
        if "Run ID" in self.run_udfs:
            self.process.udf["Run ID"] = self.run_udfs["Run ID"]
            set_field(self.process)

    def _get_cycles(self):
        """To find out if it was a single end or paired end run."""
        if "Read 1 Cycles" in self.run_udfs:
            self.read_length = self.run_udfs["Read 1 Cycles"]
        else:
            sys.exit("Could not get 'Read 1 Cycles' from the sequencing step.")
        if "Read 2 Cycles" in self.run_udfs:
            self.single = False

    def _get_file_path(self, cont_name):
        """File path for HiSeq, Miseq and Xten"""
        if self.run_type == "MiSeq":
            path_id = self.run_udfs["Flow Cell ID"]
            data_folder = "miseq_data"
        elif self.run_type == "HiSeqX10":
            path_id = cont_name
            data_folder = "HiSeq_X_data"
        else:
            data_folder = "hiseq_data"
            path_id = cont_name
        try:
            self.file_path = glob.glob(
                f"/srv/ngi-nas-ns/{data_folder}/*{path_id}/Unaligned/Basecall_Stats_*/"
            )[0]
        except:
            sys.exit("Failed to get file path")

    def _get_demultiplex_files(self):
        """Files are read from the file msf system. Path hard coded."""
        FRMP = FlowcellRunMetricsParser()
        try:
            fp_dem = self.file_path + "Demultiplex_Stats.htm"
            self.dem_stat = FRMP.parse_demultiplex_stats_htm(fp_dem)
            logging.info(f"Parsed file {fp_dem}")
        except:
            sys.exit("Failed to find or parse Demultiplex_Stats.htm.")
        try:
            fp_und = self.file_path + "Undemultiplexed_stats.metrics"
            self.undem_stat = FRMP.parse_undemultiplexed_barcode_metrics(fp_und)
            logging.info(f"Parsed file {fp_und}")
        except:
            sys.exit("Failed to find or parse Undemultiplexed_stats.metrics")

    def _get_threshold_Q30(self):
        if "Threshold for % bases >= Q30" in self.user_def_tresh:
            self.Q30_treshold = self.user_def_tresh["Threshold for % bases >= Q30"]
            qc_logg = f"THRESHOLD FOR %Q30 was set by user to {self.Q30_treshold}."
            print(qc_logg, file=self.qc_log_file)
        else:
            warning = (
                f"Un recognized read length: {self.read_length}. Report this to "
                "developers! set Threshold for % bases >= Q30 if you want to "
                "run bcl conversion and demultiplexing anyway."
            )
            if self.run_type == "MiSeq":
                if self.read_length < 101:
                    Q30_threshold = 80
                elif self.read_length == 101:
                    Q30_threshold = 75
                elif self.read_length == 151:
                    Q30_threshold = 70
                elif self.read_length >= 251:
                    Q30_threshold = 60
                else:
                    sys.exit(warning)
            else:
                if self.read_length == 51:
                    Q30_threshold = 85
                elif self.read_length == 101:
                    Q30_threshold = 80
                elif self.read_length == 126:
                    Q30_threshold = 80
                elif self.read_length == 151:
                    Q30_threshold = 75
                else:
                    sys.exit(warning)
            qc_logg = (
                f"THRESHOLD FOR %Q30 was set to {Q30_threshold}. Value based on read"
                f" length: {self.read_length}, and run type {self.run_type}."
            )
            print(qc_logg, file=self.qc_log_file)
            self.Q30_treshold = Q30_threshold

    def run_QC(self):
        for pool in self.input_pools:
            outarts_per_lane = self.process.outputs_per_input(pool.id, ResultFile=True)
            lane_number = "1" if self.run_type == "MiSeq" else pool.location[1][0]
            LQC = LaneQC(
                lane_number,
                outarts_per_lane,
                self.run_type,
                self.undem_stat,
                self.dem_stat,
                self.single,
                self.Q30_treshold,
                self.qc_log_file,
                self.user_def_tresh,
                self.read_length,
            )
            LQC.set_and_log_tresholds()
            LQC.lane_QC()
            self.nr_lane_samps_tot += LQC.nr_lane_samps
            self.nr_lane_samps_updat += LQC.nr_samps_updat
            self.QC_fail += LQC.QC_fail
            if LQC.high_lane_yield:
                self.high_lane_yield.append(LQC.lane)
            if LQC.high_index_yield:
                self.high_index_yield.append(LQC.lane)

        ## moove this part to logging???-->>
        if self.high_index_yield or self.high_lane_yield:
            warn = "WARNING: "
            if self.high_index_yield:
                self.high_index_yield = ", ".join(list(set(self.high_index_yield)))
                warn = (
                    f"{warn} High yield of unexpected index on lane(s): {self.high_index_yield} ."
                    ""
                )
            if self.high_lane_yield:
                self.high_lane_yield = ", ".join(list(set(self.high_lane_yield)))
                warn = (
                    f"{warn} High total yield of unexpected index on lane(s): "
                    f"{self.high_lane_yield}."
                )
            warn = warn + "Please check the Metrics file!"
            self.abstract.insert(0, warn)

    def make_demultiplexed_counts_file(self, demuxfile):
        """Reformats the content of the demultiplex and undemultiplexed files
        to be more easy to read."""

        demuxfile = demuxfile + ".csv"
        keys = [
            "Project",
            "Sample ID",
            "Lane",
            "# Reads",
            "Index",
            "Index name",
            "% of >= Q30 Bases (PF)",
        ]
        toCSV = []
        for pool in self.input_pools:
            if self.run_type == "MiSeq":
                lane = "1"
            else:
                lane = pool.location[1][0]
            for row in self.dem_stat["Barcode_lane_statistics"]:
                if row["Lane"] == lane:
                    row_dict = {x: row[x] for x in keys if x in row}
                    row_dict["Index name"] = ""
                    toCSV.append(row_dict)
            if lane in list(self.undem_stat.keys()):
                undet_per_lane = self.undem_stat[lane]["undemultiplexed_barcodes"]
                nr_undet = len(undet_per_lane["count"])
                for row in range(nr_undet):
                    row_dict = {x: "" for x in keys}
                    row_dict["# Reads"] = undet_per_lane["count"][row]
                    row_dict["Index"] = undet_per_lane["sequence"][row]
                    row_dict["Index name"] = undet_per_lane["index_name"][row]
                    row_dict["Lane"] = undet_per_lane["lane"][row]
                    toCSV.append(row_dict)
        try:
            f = open(demuxfile, "w")
            dict_writer = csv.DictWriter(f, keys, dialect="excel")
            dict_writer.writer.writerow(keys)
            dict_writer.writerows(toCSV)
            f.close
            self.abstract.append(
                "INFO: A Metrics file has been created with "
                "demultiplexed and undemultiplexed counts for debugging."
            )
        except:
            self.abstract.append(
                "WARNING: Could not generate a Metrics file "
                "with demultiplexed and undemultiplexed counts."
            )

    def logging(self):
        """Collects and prints logging info."""
        self.abstract.append(
            f"INFO: QC-data found and QC-flags uploaded for {self.nr_lane_samps_updat}"
            f" out of {self.nr_lane_samps_tot} analytes."
        )
        if self.QC_fail:
            self.abstract.append(
                "Failed to make qc for samples: {}".format(
                    ", ".join(list(set(self.QC_fail)))
                )
            )
        if "WARNING" in " ".join(self.abstract):
            sys.exit(" ".join(self.abstract))
        else:
            print(" ".join(self.abstract), file=sys.stderr)


class LaneQC:
    def __init__(
        self,
        lane_number,
        out_arts,
        run_type,
        undem_stat,
        dem_stat,
        single,
        Q30_treshold,
        qc_log_file,
        user_def_tresholds,
        read_length,
    ):
        ##  Output artifacts and user defined tresholds
        self.out_arts = out_arts
        self.user_def_tresh = user_def_tresholds

        ##  Info from files in file system
        self.counts = undem_stat[lane_number]["undemultiplexed_barcodes"]["count"]
        self.BLS = dem_stat["Barcode_lane_statistics"]

        ##  Tresholds
        self.exp_lane_clust = None
        self.reads_threshold = None
        self.un_exp_lane = None
        self.thres_un_exp_ind = None
        self.Q30_treshold = Q30_treshold

        ##  Stuff for logging
        self.qc_log_file = qc_log_file
        self.high_lane_yield = False
        self.high_index_yield = False
        self.nr_samps_updat = 0
        self.html_file_error = False
        self.QC_fail = []

        ##  Other variables
        self.single = single
        self.read_length = read_length
        self.run_type = run_type
        self.lane = lane_number
        self.nr_lane_samps = len(out_arts)

    def set_and_log_tresholds(self):
        """Generating tresholds and writing the tresholds to log file."""
        print("", file=self.qc_log_file)
        print(f"TRESHOLDS - LANE {self.lane}:", file=self.qc_log_file)
        self._get_exp_lane_and_ind_clust()
        self._set_reads_threshold()
        self._set_tresh_un_exp_lane()
        self._set_tresh_un_exp_ind()

    def _get_exp_lane_and_ind_clust(self):
        """The expected number of lane clusters depends on run type and run mode.
        The expected number of clusters per index (sample) is then:
        (expected nr clusters on the lane)/(nr indexes on the lane)"""
        if self.run_type == "MiSeq":
            if self.read_length in [76, 301]:
                self.exp_lane_clust = 18000000
            else:
                self.exp_lane_clust = 10000000
        elif self.run_type in [
            "HiSeq Rapid Flow Cell v1",
            "HiSeq Rapid Flow Cell v2",
            "TruSeq Rapid Flow Cell v2",
            "TruSeq Rapid Flow Cell v1",
        ]:
            self.exp_lane_clust = 114000000
        elif self.run_type == "HiSeq Flow Cell v3":
            self.exp_lane_clust = 143000000
        elif self.run_type == "HiSeq Flow Cell v4":
            self.exp_lane_clust = 188000000
        elif self.run_type == "HiSeqX10":
            self.exp_lane_clust = 250000000
        else:
            sys.exit(
                f"Unrecognized run type: {self.run_type}. Report to developer! Set "
                "Threshold for # Reads if you want to run bcl conversion "
                "and demultiplexing again"
            )
        self.exp_samp_clust = np.true_divide(self.exp_lane_clust, self.nr_lane_samps)

    def _set_reads_threshold(self):
        """Treshold for nr reads per index:
        exp_samp_clust*0.5"""

        if "Threshold for # Reads" in self.user_def_tresh:
            self.reads_threshold = self.user_def_tresh["Threshold for # Reads"]
            qc_logg = f"Index yield - expected index: {self.reads_threshold}"
            print(qc_logg, file=self.qc_log_file)
        else:
            self.reads_threshold = int(self.exp_samp_clust * 0.5)
            qc_logg = (
                f"Index yield - expected index: {self.reads_threshold}. Value based on nr of "
                f"sampels in the lane: {self.nr_lane_samps}, and run type {self.run_type}."
            )
            print(qc_logg, file=self.qc_log_file)

    def _set_tresh_un_exp_lane(self):
        """Treshold for un expected indexes on a hole lane:
        If single end: exp_lane_clust*0.05
        If paired end: exp_lane_clust*0.1"""

        if self.single:
            self.un_exp_lane = int(self.exp_lane_clust * 0.05)
        else:
            self.un_exp_lane = int(self.exp_lane_clust * 0.1)
        qc_logg = (
            "Lane yield - un expected index: {}. Value based on run "
            "type {}, and run setings: {}".format(
                self.un_exp_lane,
                self.run_type,
                "Single End" if self.single else "Paired End",
            )
        )
        print(qc_logg, file=self.qc_log_file)

    def _set_tresh_un_exp_ind(self):
        """Threshold for un expected index:
        exp_samp_clust*0.1"""

        if "Threshold for Undemultiplexed Index Yield" in self.user_def_tresh:
            self.thres_un_exp_ind = self.user_def_tresh[
                "Threshold for Undemultiplexed Index Yield"
            ]
            qc_logg = (
                f"Index yield - un expected index: {self.thres_un_exp_ind}. Value set by user."
                ""
            )
            print(qc_logg, file=self.qc_log_file)
        else:
            self.thres_un_exp_ind = int(self.exp_samp_clust * 0.1)
            qc_logg = (
                f"Index yield - un expected index: {self.thres_un_exp_ind}. Value set to 10% "
                "of expected index yield"
            )
            print(qc_logg, file=self.qc_log_file)

    def lane_QC(self):
        for target_file in self.out_arts:
            samp_name = target_file.samples[0].name
            for lane_samp in self.BLS:
                if self.lane == lane_samp["Lane"]:
                    samp = lane_samp["Sample ID"]
                    if samp == samp_name:
                        IQC = IndexQC(target_file, lane_samp)
                        IQC.set_target_file_udfs()
                        IQC.set_read_pairs(self.single)
                        try:
                            IQC.lane_index_QC(self.reads_threshold, self.Q30_treshold)
                            if IQC.html_file_error:
                                self.html_file_error = IQC.html_file_error
                            set_field(IQC.t_file)
                            self.nr_samps_updat += 1
                        except:
                            self.QC_fail.append(samp)
        self._check_un_exp_lane_yield()
        for index_count in self.counts:
            self._check_un_exp_ind_yield(index_count)

    def _check_un_exp_lane_yield(self):
        unexp_lane_yield = sum([int(x) for x in self.counts])
        if unexp_lane_yield > self.un_exp_lane:
            self.high_lane_yield = True

    def _check_un_exp_ind_yield(self, index_count):
        if int(index_count) > self.thres_un_exp_ind:
            self.high_index_yield = True


class IndexQC:
    def __init__(self, t_file, stats):
        self.t_file = t_file
        self.stats = stats
        self.samp_udfs = {}
        self.html_file_error = False

    def set_target_file_udfs(self):
        """Populates the target file udfs (run lane index resolution) with
        Barcode lane statistics fetched from demultiplexed stats file"""

        self.samp_udfs = {
            "%PF": self.stats["% PF"],
            "% One Mismatch Reads (Index)": self.stats["% One Mismatch Reads (Index)"],
            "% of Raw Clusters Per Lane": self.stats["% of raw clusters per lane"],
            "Ave Q Score": self.stats["Mean Quality Score (PF)"],
            "% Perfect Index Read": self.stats["% Perfect Index Reads"],
        }
        self._make_float()
        self._set_Yield_Mbases()
        self._set_Q30()
        self._set_reads()

    def _make_float(self):
        for key, val in self.samp_udfs.items():
            try:
                self.t_file.udf[key] = float(val)
            except:
                self.html_file_error = True

    def _set_Yield_Mbases(self):
        try:
            yMb = self.stats["Yield (Mbases)"].replace(",", "")
            self.t_file.udf["Yield PF (Gb)"] = np.true_divide(float(yMb), 1000)
        except:
            self.html_file_error = True

    def _set_Q30(self):
        if "% Bases >=Q30" not in dict(list(self.t_file.udf.items())):
            self.t_file.udf["% Bases >=Q30"] = self.stats["% of >= Q30 Bases (PF)"]

    def _set_reads(self):
        if "# Reads" not in dict(list(self.t_file.udf.items())):
            try:
                self.t_file.udf["# Reads"] = float(
                    self.stats["# Reads"].replace(",", "")
                )
            except:
                self.html_file_error = True

    def set_read_pairs(self, single):
        if single and self.t_file.udf["# Reads"]:
            self.t_file.udf["# Read Pairs"] = int(self.t_file.udf["# Reads"])
        elif self.t_file.udf["# Reads"]:
            reads = float(self.t_file.udf["# Reads"])
            self.t_file.udf["# Read Pairs"] = np.true_divide(reads, 2)

    def lane_index_QC(self, reads_threshold, Q30_treshold):
        """Index-lane level QC based on derieved QC-tresholds. OBS: Fetches info
        from target file udf if they are already set. Otherwise from file
        system!!! This is to take into account yield and quality after quality
        filtering if performed."""

        Q30 = float(self.stats["% of >= Q30 Bases (PF)"])
        QC2 = Q30 >= Q30_treshold
        QC3 = self.t_file.udf["# Read Pairs"] >= reads_threshold
        if QC2 and QC3:
            self.t_file.udf["Include reads"] = "YES"
            self.t_file.qc_flag = "PASSED"
        else:
            self.t_file.udf["Include reads"] = "NO"
            self.t_file.qc_flag = "FAILED"


def main(lims, pid, epp_logger, demuxfile, qc_log_file):
    process = Process(lims, id=pid)
    RQC = RunQC(process)
    RQC.make_qc_log_file(qc_log_file)
    RQC.get_run_info()
    RQC.run_QC()
    RQC.make_demultiplexed_counts_file(demuxfile)
    RQC.logging()


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid", default=None, dest="pid", help="Lims id for current Process"
    )
    parser.add_argument(
        "--log",
        dest="log",
        help=("File name for standard log file, for runtime information and problems."),
    )
    parser.add_argument(
        "--file",
        dest="file",
        default="demux",
        help=("File path to demultiplexed metrics file"),
    )
    parser.add_argument(
        "--qc_log_file",
        dest="qc_log_file",
        default="qc_log_file",
        help=("File path to qc logfile file"),
    )
    args = parser.parse_args()
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()

    with EppLogger(log_file=args.log, lims=lims, prepend=True) as epp_logger:
        main(lims, args.pid, epp_logger, args.file, args.qc_log_file)
