#!/usr/bin/env python
DESC = """EPP script to copy the "Comments" field to the projects running notes on process termination

Denis Moreno, Science for Life Laboratory, Stockholm, Sweden
"""
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from write_notes_to_couchdb import write_note_to_couch


def categorization(process_name):
    decision = {
        "10X Chromium Generic: cDNA QC": "Workset",
        "10X Chromium Generic: Finish GEM Processing": "Workset",
        "10X Chromium Generic: Finish Library Processing": "Workset",
        "10X Chromium Generic: GEM Cleanup": "Workset",
        "10X Chromium Generic: GEM Generation": "Workset",
        "10X Chromium Generic: Library Construction": "Workset",
        "10X Chromium Generic: Library Purification": "Workset",
        "10X Chromium Generic: Pre-GEM Processing": "Workset",
        "10X Chromium Generic: Pre-Library Processing": "Workset",
        "Adapter ligation and reverse transcription": "Workset",
        "Aggregate QC (DNA) 4.0": "Workset",
        "Aggregate QC (Library Validation) 4.0": "Workset",
        "Aggregate QC (RNA) 4.0": "Workset",
        "Aliquot Libraries for Hybridization (SS XT)": "",
        "Aliquot Libraries for Pooling (Small RNA)": "",
        "Aliquot Samples for Caliper/Bioanalyzer": "Workset",
        "Aliquot Samples for Qubit/Bioanalyzer": "Workset",
        "Amplification and Purification": "Workset",
        "Amplify Adapter-Ligated Library (SS XT) 4.0": "",
        "Amplify by PCR and Add Index Tags (TruSeq small RNA) 1.0": "Workset",
        "Amplify by PCR and Add Index Tags": "Workset",
        "Amplify Captured Libraries to Add Index Tags (SS XT) 4.0": "",
        "AMPure Size Selection": "Workset",
        "Applications Generic Process": "Workset",
        "Applications Indexing": "Workset",
        "Applications Pre-Pooling": "",
        "Automated Quant-iT QC (DNA) 4.0": "Workset",
        "Automated Quant-iT QC (Library Validation) 4.0": "Workset",
        "Bcl Conversion & Demultiplexing (Illumina SBS) 4.0": "Bioinformatics",
        "Bioanalyzer Fragmentation QC (TruSeq DNA) 4.0": "Workset",
        "Bioanalyzer QC (DNA) 4.0": "Workset",
        "Bioanalyzer QC (Library Validation) 4.0": "Workset",
        "Bioanalyzer QC (RNA) 4.0": "Workset",
        "CA Purification": "Workset",
        "CaliperGX QC (DNA)": "",
        "CaliperGX QC (RNA)": "",
        "Capture And Wash (SS XT) 4.0": "Workset",
        "cDNA QC": "",
        "Chromatin capture, digestion, end ligation and crosslink reversal (HiC) 1.0": "",
        "Circularization": "Workset",
        "Cluster Generation (HiSeq X) 1.0": "Flowcell",
        "Cluster Generation (Illumina SBS) 4.0": "Flowcell",
        "Crosslinking & Digestion": "",
        "Customer Gel QC": "Workset",
        "Denature, Dilute and Load Sample (MiSeq) 4.0": "Flowcell",
        "End repair, A-tailing and adapter ligation (Nextera) 4.0": "Workset",
        "End Repair, A-Tailing and Adapter Ligation (SS XT) 4.0": "Workset",
        "End repair, A-tailing and adapter ligation (TruSeq RNA) 4.0": "Workset",
        "End repair, adapter ligation, ligation capture and Index PCR (HiC)": "",
        "End repair, size selection, A-tailing and adapter ligation (TruSeq PCR-free DNA) 4.0": "Workset",
        "Enrich DNA fragments (TruSeq RNA) 4.0": "Workset",
        "Fragment Analyzer QC (DNA) 4.0": "",
        "Fragment Analyzer QC (Library Validation) 4.0": "",
        "Fragment Analyzer QC (RNA) 4.0": "",
        "Fragment DNA (ThruPlex)": "Workset",
        "Fragment DNA (TruSeq DNA) 4.0": "Workset",
        "Fragmentation & cDNA synthesis (TruSeq RNA) 4.0": "Workset",
        "g-Tube Fragmentation": "Workset",
        "Generic QC": "",
        "HiC Intermediate QC": "",
        "HT-End repair, A-tailing and adapter ligation (TruSeq RNA) 4.0": "Workset",
        "Hybridize Library  (SS XT) 4.0": "",
        "Illumina DNA No-QC Library Construction": "Workset",
        "Illumina DNA No-QC Library Pooling": "Workset",
        "Illumina DNA PCR-free Library Construction": "Workset",
        "Illumina Sequencing (HiSeq X) 1.0": "Flowcell",
        "Illumina Sequencing (Illumina SBS) 4.0": "Flowcell",
        "Illumina Sequencing (NextSeq) v1.0": "Flowcell",
        "Intermediate QC": "",
        "Library Normalization (AVITI) v1.0": "",
        "Library Normalization (Illumina SBS) 4.0": "",
        "Library Normalization (MiSeq) 4.0": "",
        "Library Normalization (NextSeq) v1.0": "",
        "Library Normalization (NovaSeq) v2.0": "",
        "Library Normalization (NovaSeqXPlus) v1.0": "",
        "Library Normalization": "",
        "Library Pooling (Finished Libraries) 4.0": "",
        "Library Pooling (Illumina SBS) 4.0": "",
        "Library Pooling (MiSeq) 4.0": "",
        "Library Pooling (NextSeq) v1.0": "",
        "Library Pooling (TruSeq Small RNA) 1.0": "",
        "Library Preparation & Amplification": "",
        "Ligate 3' adapters (TruSeq small RNA) 1.0": "Workset",
        "Ligate 5' adapters (TruSeq small RNA) 1.0": "Workset",
        "Linear DNA digestion, Circularized DNA shearing and Streptavidin Bead Binding": "Workset",
        "Load to Flowcell (NextSeq v1.0)": "",
        "MinElute Purification": "",
        "MinION QC": "Workset",
        "MiSeq Run (MiSeq) 4.0": "Flowcell",
        "mRNA Purification, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0": "Workset",
        "NeoPrep Library QC v1.0": "Workset",
        "ONT Adapter ligation v2.0": "Workset",
        "ONT Barcoding": "Workset",
        "ONT End-Prep v2.0": "Workset",
        "ONT Finish Sequencing v2.1": "Flowcell",
        "ONT Pooling v2.0": "Workset",
        "ONT QC Adapter Ligation v1.2": "",
        "ONT QC Analysis": "",
        "ONT QC Barcoding v1.2": "",
        "ONT QC End-prep v1.2": "",
        "ONT QC Pooling v1.2": "",
        "ONT QC Setup Workset/Plate": "",
        "ONT QC Start Sequencing": "Flowcell",
        "ONT Sequencing and Reloading v3.1": "Flowcell",
        "ONT Start Sequencing v3.0": "Flowcell",
        "Pre-Pooling (AVITI) v1.0": "",
        "Pre-Pooling (Illumina SBS) 4.0": "",
        "Pre-Pooling (MiSeq) 4.0": "",
        "Pre-Pooling (NextSeq) v1.0": "",
        "Pre-Pooling (NovaSeq) v2.0": "",
        "Pre-Pooling (NovaSeqXPlus) v1.0": "",
        "Pre-Pooling": "",
        "Project Summary 1.3": "",
        "PromethION Sequencing": "",
        "Proximity Ligation": "",
        "Purification (ThruPlex)": "Workset",
        "Purification": "Workset",
        "qPCR QC (Dilution Validation) 4.0": "",
        "qPCR QC (Library Validation) 4.0": "",
        "Quant-iT QC (DNA) 4.0": "Workset",
        "Quant-iT QC (Library Validation) 4.0": "Workset",
        "Quant-iT QC (RNA) 4.0": "Workset",
        "Qubit QC (Dilution Validation) 4.0": "Workset",
        "Qubit QC (DNA) 4.0": "Workset",
        "Qubit QC (Library Validation) 4.0": "Workset",
        "Qubit QC (RNA) 4.0": "Workset",
        "RAD-seq Library Indexing v1.0": "",
        "Reverse Transcribe (TruSeq small RNA) 1.0": "Workset",
        "RiboZero depletion": "Workset",
        "Sample Crosslinking": "",
        "Sample Inspection": "",
        "Sample Placement (Size Selection)": "",
        "Sample Setup": "",
        "Selection, cDNA Synthesis and Library Construction": "Workset",
        "Setup Workset/Plate": "Workset",
        "Shear DNA (SS XT) 4.0": "Workset",
        "Size Selection (Caliper XT) 1.0": "Workset",
        "Size Selection (Pippin)": "Workset",
        "Size Selection (Robocut)": "Workset",
        "Sort AVITI Samples v1.0": "",
        "Sort HiSeq Samples (HiSeq) 4.0": "",
        "Sort HiSeq X Samples (HiSeq X) 1.0": "",
        "Sort MiSeq Samples (MiSeq) 4.0": "",
        "Sort NextSeq Samples (NextSeq) v1.0": "",
        "Sort NovaSeq Samples (NovaSeq) v2.0": "",
        "Sort NovaSeqXPlus Samples (NovaSeqXPlus) v1.0": "",
        "Sort Samples (ONT Pre-Prep)": "Workset",
        "Sort Samples (ONT Seq)": "Workset",
        "Sort Samples for Norm/Pooling": "",
        "Tagmentation, Strand displacement and AMPure purification": "Workset",
        "ThruPlex library amplification": "Workset",
        "ThruPlex template preparation and synthesis": "Workset",
        "Tissue Extraction": "",
        "Tissue QC": "",
        "Volume Measurement QC": "Workset",
    }

    return decision[process_name]


def main(lims, args):
    noteobj = {}
    pro = Process(lims, id=args.pid)
    if "Comments" in pro.udf and pro.udf["Comments"] != "":
        if isinstance(pro.udf["Comments"], str):
            comments = pro.udf["Comments"]
        else:
            comments = pro.udf["Comments"].encode("utf-8")
        note = "Comment from {} ({}) : \n{}".format(
            pro.type.name,
            "[LIMS]({}/clarity/work-details/{})".format(BASEURI, pro.id.split("-")[1]),
            comments,
        )
        noteobj["note"] = note
        noteobj["email"] = pro.technician.email
        noteobj["categories"] = [categorization(pro.type.name)]
        noteobj["note_type"] = "project"

        # find the correct projects.
        samples = set()
        projects = set()
        for inp in pro.all_inputs():
            # bitwise or to add inp.samples to samplesas a set
            samples |= set(inp.samples)
        for sam in samples:
            if sam.project:
                projects.add(sam.project)

        for proj in projects:
            write_note_to_couch(proj.id, noteobj, lims.get_uri())


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()

    main(lims, args)
