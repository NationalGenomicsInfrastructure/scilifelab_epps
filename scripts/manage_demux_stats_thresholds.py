"""Written by Isak Sylvin. isak.sylvin@scilifelab.se"""

import sys
import logging

class Thresholds():
    def __init__(self, instrument, chemistry, paired, read_length):
        self.logger = logging.getLogger('demux_logger.thresholds')
        self.Q30 = None
        self.exp_lane_clust = None
        self.undet_indexes_perc = None
        
        #Checks that only supported values are entered
        self.valid_instruments = ["miseq", "hiseq", "HiSeq_X"]
        self.valid_chemistry = ["MiSeq", "HiSeq Rapid Flow Cell v1","HiSeq Rapid Flow Cell v2", 
                             "TruSeq Rapid Flow Cell v2", "TruSeq Rapid Flow Cell v3", "HiSeq Flow Cell v4", "HiSeqX v2.5"]
        
        if not instrument in self.valid_instruments or not chemistry in self.valid_chemistry:
            self.problem_handler("exit", "Detected instrument and chemistry combination are not classed as valid in manage_demux_stats_thresholds.py")
        else:
            self.instrument = instrument
            self.chemistry = chemistry
            self.paired = paired
            self.read_length = read_length
            self.set_undet_indexes_perc()
        
    def problem_handler(self, type, message):
        if type == "exit":
            self.logger.error(message)
            sys.exit(message)
        elif type == "warning":
            self.logger.warning(message)
            sys.stderr.write(message)
        else:
            self.logger.info(message)
    
    """Q30 values are derived from governing document 1244:4"""
    def set_Q30(self):
        if self.instrument == "miseq":
            if self.chemistry == "MiSeq":
                if self.read_length >= 250:
                    self.Q30 = 60
                elif self.read_length >= 150:
                    self.Q30 = 70
                elif self.read_length >= 100:
                    self.Q30 = 75
                elif self.read_length < 100:
                    self.Q30 = 80
        elif self.instrument == "hiseq":
            #Rapid run flowcell
            if self.chemistry in ["HiSeq Rapid Flow Cell v1","HiSeq Rapid Flow Cell v2", "TruSeq Rapid Flow Cell v2", "TruSeq Rapid Flow Cell v1"] :
                if self.read_length >= 150:
                    self.Q30 = 75
                elif self.read_length >= 100:
                    self.Q30 = 80
                elif self.read_length < 100:
                    self.Q30 = 85
            #v3
            elif self.chemistry == "HiSeq Flow Cell v3":
                if self.read_length >= 100 and self.paired:
                    self.Q30 = 80
                elif self.read_length >= 50 and self.paired:
                    self.Q30 = 85
            #v4
            elif self.chemistry == "HiSeq Flow Cell v4":
                if self.read_length >= 125 and self.paired:
                    self.Q30 = 80
                elif self.read_length >= 50:
                    self.Q30 = 85
                        
        elif self.instrument == "HiSeq_X":
            if self.chemistry == "HiSeqX v2.5":
                if self.read_length >= 150 and self.paired:
                    self.Q30 = 75
        if not self.Q30:
            self.problem_handler("exit", "No predefined Q30 threshold (see doc 1244). Instrument: {}, Chemistry: {}, Read Length: {}".\
                                 format(self.instrument, self.chemistry, self.read_length))
    
    """Expected lanes per cluster are derived from undemultiplex_index.py"""
    def set_exp_lane_clust(self):
        if self.instrument == "miseq":
            if self.chemistry == "MiSeq":
                if self.read_length >= 76 and self.read_length <= 301:  
                    self.exp_lane_clust = 18000000
                else:                               
                    self.exp_lane_clust = 10000000
        elif self.instrument == "hiseq":
            #Rapid run flowcell
            if self.chemistry in ["HiSeq Rapid Flow Cell v1","HiSeq Rapid Flow Cell v2", "TruSeq Rapid Flow Cell v2", "TruSeq Rapid Flow Cell v1"] :
                self.exp_lane_clust = 114000000
            #v3
            elif self.chemistry == "HiSeq Flow Cell v3":
               self.exp_lane_clust = 143000000
            #v4
            elif self.chemistry == "HiSeq Flow Cell v4":
                self.exp_lane_clust = 188000000
        elif self.instrument == "HiSeq_X":
            #HiSeqX runs are always paired!
            if self.paired:
                #X v2.5 (common)
                if self.chemistry == "HiSeqX v2.5":
                    self.exp_lane_clust = 320000000
                #X v2.0 (rare)
                elif self.chemistry == "HiSeqX v2.0":
                    self.exp_lane_clust = 305000000
            else:
                self.problem_handler("exit", "HiSeqX runs should always be paired but script has detected otherwise. Something has gone terribly wrong.")
        if not self.exp_lane_clust:
            self.problem_handler("exit", "No predefined clusters per lane threshold. Instrument: {}, Chemistry: {}, Read Length: {}".\
                                 format(self.instrument, self.chemistry, self.read_length))
            
    """Maximum undetermined per lane are derived from the current taca.yaml"""      
    def set_undet_indexes_perc(self):
        self.undet_indexes_perc = 5
    
    
            
            
