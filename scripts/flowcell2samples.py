#!/usr/bin/env python
import os, sys, time, argparse, re
from fnmatch import fnmatch

## python 3 : scandir
try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

ROOT = "/scicore/projects/openbis/userstore/"

METADATA = {
    "SampleName": r"\nSAMPLE_CODE\s+(GFB[\w\-_]+)\n",
    #'SampleName' : r'\nSAMPLE_CODE\s+(BSSE[\w\-_]+)\n',
    "ExternalSampleName": r"\nEXTERNAL_SAMPLE_NAME\s+([^\n]+)\n",
    "SampleType": r"\nSAMPLE_TYPE\s+([^\n]+)\n",
    #'SampleType' : r'\nSAMPLE_KIND\s+([^\n]+)\n',
    "Kit": r"\nKIT\s+([^\n]+)\n",
    "Species": r"\nNCBI_ORGANISM_TAXONOMY\s+([\w]+)\n",
    "EndType": r"\nEND_TYPE\s+([^\n]+)\n",
    "NCycles": r"\nR1_CYCLE\s+([\w]+)\n",
    #'NCycles' : r'\nCYCLES_REQUESTED_BY_CUSTOMER\s+([\w]+)\n',
    "ContactPerson": r"\nCONTACT_PERSON_NAME\s+([^\n]+)\n",
    #'FlowcellID' : r'\nSAMPLE_CODE\s+([\w]+)\n',
    "FlowcellID": r"\nFLOWCELL_NAME\s+([\w]+)\n",
    "Fastq": r"\nFASTQ_FILES\n([\w\-\_\.\t\n]+)[\n]+",
}

##################################################


def getMetaFileData(metaFile, METADATA):
    patterns = dict(
        [(key, re.compile("{}".format(value))) for key, value in METADATA.items()]
    )
    fields = dict([(key, "") for key in METADATA.keys()])
    root = os.path.dirname(metaFile)
    with open(metaFile) as fn:
        lines = fn.read()
    for key, pattern in patterns.items():
        m = pattern.search(lines)
        if m:
            if key == "Fastq":
                fqs = [
                    os.path.join(root, elm.split()[0])
                    for elm in m.group(1).splitlines()
                ]
                fields[key] = fqs
                fields["FileSize"] = {fq: os.path.getsize(fq) / 2**30 for fq in fqs}
            elif key == "SampleName":
                fields[key] = m.group(1).replace("-", "_")
            else:
                fields[key] = m.group(1)
        else:
            fields[key] = ""
    return fields


##################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "flowcell",
        help="Flowcell ID(s). Specify '*' to search for any flowcell ID.",
        nargs="+",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-g", "--group", help="Group name, including prefix, e.g. dbm_ or biozentrum_."
    )
    group.add_argument("-p", "--path", help="ROOT path")
    parser.add_argument(
        "-d",
        "--date",
        help="Minimal date. Format: day/month/year",
        default="01/01/1980",
    )
    parser.add_argument("-o", "--output", help="Output file name")

    args = parser.parse_args()

    if args.group:
        ROOT = os.path.join(ROOT, args.group)
    elif args.path:
        ROOT = args.path
    else:
        pass  ##  ROOT = ROOT

    ## format date
    startDate = time.mktime(
        time.strptime(args.date, "%d/%m/%Y")
    )  ## seconds since Epoch

    ## flowcell file patters
    #    pattern = ["BSSE*{}*metadata.tsv".format(fc) for fc in args.flowcell]
    pattern = [f"*{fc}*.fastq.gz" for fc in args.flowcell]

    ## get files
    result = []
    for root, dirs, files in walk(ROOT):
        if any([fnmatch(file, pat) for pat in pattern for file in files]):
            fpath = [
                os.path.join(root, file)
                for file in files
                if fnmatch(file, "*metadata.tsv")
            ]
            if len(fpath) != 1:
                sys.stderr.write(
                    f"WARNING: no/multiple *metadata.tsv found! dir: {root}\n"
                )
                continue  # skip directory
            fpath = fpath[0]
            if os.path.getmtime(fpath) > startDate:
                result.append(getMetaFileData(fpath, METADATA))

    if len(result) == 0:
        sys.stderr.write("ERROR: no files found!\n")
        sys.exit()

    ## order result by id
    id = [elm["SampleName"] for elm in result]
    result = [result[idx] for idx in sorted(range(len(id)), key=lambda k: id[k])]

    ## file size
    FileSize = {k: v for d in [elm["FileSize"] for elm in result] for k, v in d.items()}

    ## single or paired or mixture?
    PAIRED = False
    if "Paired End" in [elm["EndType"] for elm in result]:
        #    if "PAIRED_END" in [elm['EndType'] for elm in result]:
        head = "FileName1\tFileName2\tFileSizeGB1\tFileSizeGB2"
        PAIRED = True
        if "SINGLE_END" in [elm["EndType"] for elm in result]:
            sys.stderr.write("Warning: Single and paired-end files detected!\n")
    else:  # single-end
        head = "FileName\tFileSizeGB"

    ## print in that order
    fields = [
        "SampleName",
        "ExternalSampleName",
        "FlowcellID",
        "SampleType",
        "Kit",
        "NCycles",
        "Species",
        "ContactPerson",
    ]
    with open(args.output, "w") if args.output else sys.stdout as out:
        out.write("{}\t{}\n".format(head, "\t".join(fields)))
        for elm in result:
            fqInfo = "\t".join([elm[field] for field in fields])
            fqFiles = elm["Fastq"]
            # only keep R1/R2 FASTQ, i.e. remove dual index files (I1/I2)
            # fqFiles = [fq for fq in fqFiles if re.search(r'_R[1,2]_[\dM\_]+.fastq.gz', fq)]
            ## For AVITI this doesn't work!
            fqFiles = [
                fq for fq in fqFiles if re.search(r"_R[1,2][_\dM\_]*.fastq.gz", fq)
            ]

            if len(fqFiles) == 1:  # single end
                if PAIRED:  # mix of single+paired: add empty slots
                    line = "{0}\t\t{1:.2f}\t\t{2}".format(
                        fqFiles[0], FileSize[fqFiles[0]], fqInfo
                    )
                else:
                    line = "{0}\t{1:.2f}\t{2}".format(
                        fqFiles[0], FileSize[fqFiles[0]], fqInfo
                    )
            else:
                if PAIRED:
                    # order R1 first
                    fqPairs = [
                        [pat.replace("**", "_R1_"), pat.replace("**", "_R2_")]
                        for pat in set([re.sub("_R[1,2]_", "**", fq) for fq in fqFiles])
                    ]
                    line = "\n".join(
                        [
                            "{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4}".format(
                                fq1, fq2, FileSize[fq1], FileSize[fq2], fqInfo
                            )
                            for fq1, fq2 in fqPairs
                        ]
                    )
                else:
                    line = "\n".join(
                        [
                            "{0}\t{1:.2f}\t{2}".format(fq, FileSize[fq], fqInfo)
                            for fq in fqFiles
                        ]
                    )
            out.write(line + "\n")
