# A full list of LIMS steps with involved instrument and details for logging
lims_process_record = {
    "Adapter ligation and reverse transcription": {
        "lims_instrument": {"dest_file": "PCR"}
    },
    "Adapter ligation and reverse transcription (TruSeq small RNA) 1.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Adapter Ligation and 1st Amplification (SMARTer Pico) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Aliquot Samples for Caliper/Bioanalyzer": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "Aliquot Samples for Qubit/Bioanalyzer": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "Aliquot Libraries for Hybridization (SS XT)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "default": {"dest_file": "Speedvac"},
    },
    "Amplification and Purification": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Amplify Adapter-Ligated Library (SS XT) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Amplify Captured Libraries to Add Index Tags (SS XT) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Amplify by PCR and Add Index Tags": {"lims_instrument": {"dest_file": "PCR"}},
    "Amplify by PCR and Add Index Tags (TruSeq small RNA) 1.0": {
        "lims_instrument": {"dest_file": "PCR"}
    },
    "AMPure Size Selection": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
        "udf_Instrument Used": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        },
    },
    "Automated Quant-iT QC (DNA) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_Spectrometry": {
            "dest_file": "Tecan",
            "details": ["Assay type", "Lot no: Quant-iT reagent kit"],
        },
    },
    "Automated Quant-iT QC (Library Validation) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_Spectrometry": {
            "dest_file": "Tecan",
            "details": ["Assay type", "Lot no: Quant-iT reagent kit"],
        },
    },
    "Automated Quant-iT QC (RNA) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_Spectrometry": {
            "dest_file": "Tecan",
            "details": ["Assay type", "Lot no: Quant-iT reagent kit"],
        },
    },
    "Bioanalyzer Fragmentation QC (TruSeq DNA) 4.0": {
        "lims_instrument": {
            "dest_file": "Bioanalyzer",
            "details": ["Lot no: Chip", "Lot no: Reagent Kit"],
        }
    },
    "Bioanalyzer QC (Library Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "Bioanalyzer",
            "details": ["Lot no: Chip", "Lot no: Reagent Kit"],
        }
    },
    "Bioanalyzer QC (DNA) 4.0": {
        "lims_instrument": {
            "dest_file": "Bioanalyzer",
            "details": ["Lot no: Chip", "Lot no: Reagent Kit"],
        }
    },
    "Bioanalyzer QC (RNA) 4.0": {
        "lims_instrument": {
            "dest_file": "Bioanalyzer",
            "details": ["Lot no: Chip", "Lot no: Reagent kit", "Lot no: Ladder"],
        }
    },
    "CA Purification": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "CaliperGX QC (DNA)": {
        "lims_instrument": {
            "dest_file": "Caliper",
            "details": ["Lot no: Chip", "Lot no: Reagent Kit"],
        }
    },
    "CaliperGX QC (RNA)": {
        "lims_instrument": {
            "dest_file": "Caliper",
            "details": ["Lot no: Chip", "Lot no: Reagent Kit", "Lot no: RNA ladder"],
        }
    },
    "Capture And Wash (SS XT) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "CytAssist Probe release and Extension": {
        "lims_instrument": {"dest_file": "CytAssist", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Decrosslinking and/or Destaining": {
        "lims_instrument": {"dest_file": "PCR"},
    },
    "Denature, Dilute and Load Sample (MiSeq) 4.0": {
        "udf_Instrument Used": {
            "dest_file": "MiSeq",
            "details": ["Flowcell ID", "RGT#s"],
        }
    },
    "Deparafinization, H&E, Tissue Imaging": {
        "lims_instrument": {"dest_file": "PCR"},
    },
    "Diluting Samples": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "End repair, adapter ligation, ligation capture and Index PCR (HiC)": {
        "lims_instrument": {"dest_file": "PCR"}
    },
    "End Repair, A-Tailing and Adapter Ligation (SS XT) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "End repair, size selection, A-tailing and adapter ligation (Lucigen NxSeq DNA) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "End repair, size selection, A-tailing and adapter ligation (TruSeq DNA Nano) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "End repair, size selection, A-tailing and adapter ligation (TruSeq PCR-free DNA) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Biomek"],
            "details": ["Processname"],
        },
        "udf_PCR machine": {"dest_file": "PCR"},
    },
    "End repair, A-tailing and adapter ligation (TruSeq RNA) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Enrich DNA fragments (Nextera) 4.0": {"lims_instrument": {"dest_file": "PCR"}},
    "Enrich DNA fragments (TruSeq DNA) 4.0": {
        "lims_instrument": {"dest_file": "PCR"},
        "udf_Bravo": {"dest_file": "Bravo", "details": ["Processname"]},
    },
    "Enrich DNA fragments (TruSeq RNA) 4.0": {"lims_instrument": {"dest_file": "PCR"}},
    "Fixation, H&E, Tissue Imaging": {
        "lims_instrument": {"dest_file": "PCR"},
    },
    "Fragment Analyzer QC (DNA) 4.0": {
        "lims_instrument": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        }
    },
    "Fragment Analyzer QC (Library Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        }
    },
    "Fragment Analyzer QC (RNA) 4.0": {
        "lims_instrument": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        }
    },
    "Fragment DNA (ThruPlex)": {
        "lims_instrument": {"dest_file": "Covaris", "details": ["Lot no: Covaris tube"]}
    },
    "Fragment DNA (TruSeq DNA) 4.0": {
        "lims_instrument": {"dest_file": "Covaris", "details": ["Lot no: Covaris tube"]}
    },
    "Fragmentation & cDNA synthesis (SMARTer Pico) 4.0": {
        "udf_PCR Machine": {"dest_file": "PCR"}
    },
    "Fragmentation & cDNA synthesis (TruSeq RNA) 4.0": {
        "udf_PCR Machine": {"dest_file": "PCR"}
    },
    "g-Tube Fragmentation": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
        "udf_Instrument Used": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        },
    },
    "GEM Generation (Chromium Genome v2)": {
        "lims_instrument": {"dest_file": "Chromium"}
    },
    "GEM Incubation (Chromium Genome v2)": {"lims_instrument": {"dest_file": "PCR"}},
    "HiC Intermediate QC": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        }
    },
    "Hybridize Library  (SS XT) 4.0": {
        "lims_instrument": {"dest_file": "PCR"},
        "udf_Instrument Used": {"dest_file": "Bravo", "details": ["Processname"]},
    },
    "Illumina DNA PCR-free Library Construction": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "Illumina DNA No-QC Library Construction": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Illumina DNA No-QC Library Pooling": {
        "lims_instrument": {"dest_file": "Mosquito", "details": ["Processname"]}
    },
    "Intermediate QC": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        }
    },
    "Library Normalization": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (AVITI) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (HiSeq X) 1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (Illumina SBS) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (MiSeq) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (NextSeq) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (NovaSeq) v2.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Normalization (NovaSeqXPlus) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Pooling (HiSeq X) 1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Pooling (Illumina SBS) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Pooling (MiSeq) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Pooling (NextSeq) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Pooling (RAD-seq) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Pooling (TruSeq Small RNA) 1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Library Preparation & Amplification": {"lims_instrument": {"dest_file": "PCR"}},
    "Library preparation (Chromium Genome v2)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Linear DNA digestion, Circularized DNA shearing and Streptavidin Bead Binding": {
        "lims_instrument": {"dest_file": "Covaris", "details": ["Lot no: Covaris tube"]}
    },
    "mRNA Purification, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "ONT Barcoding": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
        "udf_Instrument Used": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        },
    },
    "ONT End-Prep v2.0": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
        "udf_Instrument Used": {
            "dest_file": "FragmentAnalyzer",
            "details": ["Lot no: Fragment Analyzer Reagents"],
        },
    },
    "ONT Pooling v2.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        },
    },
    "ONT Reverse Transcription": {
        "lims_instrument": {
            "dest_file": "PCR",
        },
    },
    "ONT Start Sequencing v3.0": {
        "lims_instrument": {
            "dest_file": "MinION-PromethION",
            "details": [
                "Processname",
                "LIMS ID (Process)",
            ],
        },
    },
    "ONT QC Adapter Ligation v1.2": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
    },
    "ONT QC Barcoding v1.2": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
    },
    "ONT QC End-prep v1.2": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        },
    },
    "ONT QC Start Sequencing": {
        "lims_instrument": {
            "dest_file": "MinION-PromethION",
            "details": [
                "Processname",
                "LIMS ID (Process)",
            ],
        },
    },
    "PCR1 (Amplicon)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR machine": {"dest_file": "PCR"},
    },
    "PCR2 (Amplicon)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR machine": {"dest_file": "PCR"},
    },
    "Post-GEM Cleanup (Chromium Genome v2)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "Pre-Pooling": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Pre-Pooling (AVITI) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Pre-Pooling (Illumina SBS) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Pre-Pooling (MiSeq) 4.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Pre-Pooling (NextSeq) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Pre-Pooling (NovaSeq) v2.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Pre-Pooling (NovaSeqXPlus) v1.0": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Purification": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "Purification (ThruPlex)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "qPCR QC (Dilution Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "CFX",
            "details": ["Lot no. qPCR reagent kit", "Lot no. Standard"],
        },
        "udf_Instrument Used": {"dest_file": "Bravo", "details": ["Processname"]},
    },
    "qPCR QC (Library Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "CFX",
            "details": ["Lot no. qPCR reagent kit", "Lot no. Standard"],
        },
        "udf_Instrument Used": {"dest_file": "Bravo", "details": ["Processname"]},
    },
    "Quant-iT QC (DNA) 4.0": {
        "lims_instrument": {
            "dest_file": "CFX",
            "details": ["Assay type", "Lot no: Quant-iT reagent kit"],
        }
    },
    "Quant-iT QC (Library Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "CFX",
            "details": ["Assay type", "Lot no: Quant-iT reagent kit"],
        }
    },
    "Quant-iT QC (RNA) 4.0": {
        "lims_instrument": {
            "dest_file": "CFX",
            "details": ["Assay type", "Lot no: Quant-iT reagent kit"],
        }
    },
    "Qubit QC (DNA) 4.0": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        }
    },
    "Qubit QC (RNA) 4.0": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        }
    },
    "Qubit QC (Dilution Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        }
    },
    "Qubit QC (Library Validation) 4.0": {
        "lims_instrument": {
            "dest_file": "Qubit",
            "details": ["Assay", "Lot no: Qubit kit"],
        }
    },
    "RAD-seq Library Indexing v1.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Ribosomal cDNA Depletion and 2nd Amplification (SMARTer Pico) 4.0": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "RiboZero depletion": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Selection, cDNA Synthesis and Library Construction": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Biomek"],
            "details": ["Processname"],
        }
    },
    "Setup Workset/Plate": {
        "lims_instrument": {
            "dest_file": ["Bravo", "Mosquito"],
            "details": ["Processname"],
        }
    },
    "Size Selection (Pippin)": {
        "default": {
            "dest_file": "Pippin",
            "details": [
                "Type: Gel Cassette",
                "Lot no: Gel Cassette",
                "Type: Marker",
                "Lot no: Marker",
                "Lot no: Electrophoresis Buffer",
            ],
        }
    },
    "Size Selection (Robocut)": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]}
    },
    "Shear DNA (SS XT) 4.0": {
        "lims_instrument": {
            "dest_file": "Covaris",
            "details": ["Lot no: Covaris tube"],
        },
        "udf_Instrument Used": {"dest_file": "Bravo", "details": ["Processname"]},
    },
    "ThruPlex library amplification": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "ThruPlex template preparation and synthesis": {
        "lims_instrument": {"dest_file": "Bravo", "details": ["Processname"]},
        "udf_PCR Machine": {"dest_file": "PCR"},
    },
    "Volume Measurement QC": {
        "lims_instrument": {"dest_file": "VC100"},
    },
}

# A full list of GDoc electronic logooks of instruments
GDoc_logbook = {
    "Bioanalyzer": {"File": "1m2_kGf-vTi-XP8RCnuxZ3yAR1xIw_ySCpGMuRlb2gFs"},
    "Biomek": {"File": "1-0KQUfcnS1ekVVWQyLurudqVBnMu1ekpYguzdQTtnRM"},
    "Bravo": {"File": "1Di5uRlEI7zlQ7DgvQuEYpwDslf8VMrRnhruk5GWNtIo"},
    "Caliper": {"File": "1x3w-0s1-xENQORthMSF1GLTezfXQcsaVEiAfJRRjsoc"},
    "CFX": {"File": "19LKni8LO-Dzkvs7gkVHLEcTqzqX2zOerT3SOF9GPHNQ"},
    "Chromium": {"File": "1PDegRtvYhUVPJPt5ROS4Ar661QSMhRWWTA8Zo3JponU"},
    "Covaris": {"File": "1wpSzEdiZcRWk1YFo59Pzt4y-AVSb0Fi9AOg2VFrQOVE"},
    "CytAssist": {"File": "1sV0qefMAixVlnxGzB2q-YfcSl8a5Y69mAl555KfhU7Q"},
    "FragmentAnalyzer": {"File": "1T4Cy3ywZvl0-kQR-QbtXzu_sErPaYymXeGMf81fqK8k"},
    "MiSeq": {"File": "1ThnEbahwm3InlF_tUJ0riyT3RImVKQINfMD4rB6VThU"},
    "Mosquito": {"File": "1ssFoSdcWV-CRK5TR--hObNkM42zJ8X3ED5q_YmU-m_o"},
    "PCR": {"File": "1YE_M4ywhr5HuQEV2DhO0oVLDPRkThhuAytAEawcdTZM"},
    "Pippin": {"File": "1cJd2Wo9GMVq0HjXrVahxF2o_I_LqIipAreWOXeWwObM"},
    "Qubit": {"File": "1-sByQA6XVrbli0V24n4CxdxogLUlRlGvykkxOpBG-_U"},
    "Speedvac": {"File": "1Dk7qPJeNmzKtHWEdNkZ4yLB0FycjREqqIhNhInZ8G9g"},
    "Tecan": {"File": "1DUBEL8DBf0lnXJjIIjowQf2PrftMo9ECXpeNrDodM4s"},
    "VC100": {"File": "1FEH2Ilyaz2FThJZ_jjlDna7s4aTuRgfxvLFpYgOPHBw"},
    "MinION-PromethION": {"File": "1OazL3MNZImTJBHK-QVy-jbxL9nN1uCoOmOMSiXclRac"},
}
