"""Written by Isak Sylvin. isak.sylvin@scilifelab.se"""

import logging
import sys


class Thresholds:
    def __init__(self, instrument, chemistry, paired, read_length):
        self.logger = logging.getLogger("demux_logger.thresholds")
        self.Q30 = None
        self.exp_lane_clust = None
        self.undet_indexes_perc = None
        self.correction_factor_for_sample_in_pool = 0.75

        # Checks that only supported values are entered
        self.valid_instruments = [
            "miseq",
            "NovaSeq",
            "NextSeq",
            "NovaSeqXPlus",
            "Aviti"
            "MiSeqi100",
        ]
        self.valid_chemistry = [
            "MiSeq",
            "Version3",
            "Version2",
            "Version2Nano",
            "Version2Micro",
            "SP",
            "S1",
            "S2",
            "S4",
            "NextSeq Mid",
            "NextSeq High",
            "NextSeq 2000 P1",
            "NextSeq 2000 P2",
            "NextSeq 2000 P3",
            "10B",
            "1.5B",
            "25B",
            "AVITI High",
            "AVITI Med",
            "AVITI Low",
            "5M"
            "25M"
            "50M"
            "100M"
        ]

        if (
            instrument not in self.valid_instruments
            or chemistry not in self.valid_chemistry
        ):
            self.problem_handler(
                "exit",
                "Detected instrument and chemistry combination are not classed as valid in manage_demux_stats_thresholds.py",
            )
        else:
            self.instrument = instrument
            self.chemistry = chemistry
            self.paired = paired
            self.read_length = read_length

    def problem_handler(self, type, message):
        if type == "exit":
            self.logger.error(message)
            sys.exit(message)
        elif type == "warning":
            self.logger.warning(message)
            sys.stderr.write(message)
        else:
            self.logger.info(message)

    def set_undet_indexes_perc(self):
        if self.instrument == "miseq" and self.paired:
            self.undet_indexes_perc = 20
        else:
            self.undet_indexes_perc = 10

    """Q30 values are derived from governing document 1618"""

    def set_Q30(self):
        if self.instrument == "miseq":
            if self.read_length >= 250:
                self.Q30 = 60
            elif self.read_length >= 150:
                self.Q30 = 70
            elif self.read_length >= 100:
                self.Q30 = 75
            elif self.read_length < 100:
                self.Q30 = 80
        # Preliminary values for MiSeqi100
        elif self.instrument == "MiSeqi100":
            if self.read_length >= 250:
                self.Q30 = 85
            elif self.read_length >= 150:
                self.Q30 = 90
            elif self.read_length >= 100:
                self.Q30 = 90
            elif self.read_length < 100:
                self.Q30 = 90

        elif self.instrument == "NovaSeq":
            if self.read_length >= 150:
                self.Q30 = 75
            elif self.read_length >= 100:
                self.Q30 = 80
            elif self.read_length < 100:
                self.Q30 = 85

        elif self.instrument == "NovaSeqXPlus":
            if self.read_length >= 150:
                self.Q30 = 75
            elif self.read_length >= 100:
                self.Q30 = 80
            elif self.read_length < 100:
                self.Q30 = 85

        elif self.instrument == "NextSeq":
            if self.read_length >= 150:
                self.Q30 = 75
            elif self.read_length >= 100:
                self.Q30 = 80
            elif self.read_length < 100:
                self.Q30 = 85

        # Preliminary values for AVITI
        elif self.instrument == "Aviti":
            if self.read_length >= 250:
                self.Q30 = 85
            else:
                self.Q30 = 90

        if not self.Q30:
            self.problem_handler(
                "exit",
                f"No predefined Q30 threshold (see doc 1618). Instrument: {self.instrument}, Chemistry: {self.chemistry}, Read Length: {self.read_length}",
            )

    def set_exp_lane_clust(self):
        """Expected lanes per cluster are derived from undemultiplex_index.py"""
        if self.instrument == "miseq":
            if self.chemistry == "Version3":
                self.exp_lane_clust = 18e6
            elif self.chemistry == "Version2":
                self.exp_lane_clust = 10e6
            elif self.chemistry == "Version2Nano":
                self.exp_lane_clust = 750000
            elif self.chemistry == "Version2Micro":
                self.exp_lane_clust = 3000000
            else:
                if self.read_length >= 76 and self.read_length <= 301:
                    self.exp_lane_clust = 18e6
                else:
                    self.exp_lane_clust = 10e6
        elif self.instrument == "NovaSeq":
            if self.chemistry == "SP":
                self.exp_lane_clust = 325e6
            elif self.chemistry == "S1":
                self.exp_lane_clust = 650e6
            elif self.chemistry == "S2":
                self.exp_lane_clust = 1650e6
            elif self.chemistry == "S4":
                self.exp_lane_clust = 2000e6
        elif self.instrument == "NovaSeqXPlus":
            if self.chemistry == "10B":
                self.exp_lane_clust = 1000e6
            elif self.chemistry == "1.5B":
                self.exp_lane_clust = 750e6
            elif self.chemistry == "25B":
                self.exp_lane_clust = 3000e6
        # Preliminary values for MiSeqi100
        elif self.instrument == "MiSeqi100":
            if self.chemistry == "5M":
                self.exp_lane_clust = 5e6
            elif self.chemistry == "25M":
                self.exp_lane_clust = 25e6
            elif self.chemistry == "50M":
                self.exp_lane_clust = 70e6
            elif self.chemistry == "100M":
                self.exp_lane_clust = 100e6
        elif self.instrument == "NextSeq":
            if self.chemistry == "NextSeq Mid":
                self.exp_lane_clust = 25e6
            elif self.chemistry == "NextSeq High":
                self.exp_lane_clust = 75e6
            elif self.chemistry == "NextSeq 2000 P1":
                self.exp_lane_clust = 100e6
            elif self.chemistry == "NextSeq 2000 P2":
                self.exp_lane_clust = 400e6
            elif self.chemistry == "NextSeq 2000 P3":
                self.exp_lane_clust = 550e6
        # Preliminary values for AVITI
        elif self.instrument == "Aviti":
            if self.chemistry == "AVITI High":
                if self.read_length >= 250:
                    self.exp_lane_clust = 300e6
                else:
                    self.exp_lane_clust = 500e6
            elif self.chemistry == "AVITI Med":
                if self.read_length >= 250:
                    self.exp_lane_clust = 100e6
                else:
                    self.exp_lane_clust = 250e6
            elif self.chemistry == "AVITI Low":
                self.exp_lane_clust = 125e6
        else:
            self.problem_handler("exit", "Unknown run type!")
        if not self.exp_lane_clust:
            self.problem_handler(
                "exit",
                f"No predefined clusters per lane threshold. Instrument: {self.instrument}, Chemistry: {self.chemistry}, Read Length: {self.read_length}",
            )
