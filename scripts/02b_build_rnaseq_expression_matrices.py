"""
GEM1 - Step 4 (continued): Build gene x sample expression matrices for the two
RNA-seq cohorts (GSE126848, GSE135251), matching the format 02_fetch_geo_cohorts.py
already produced for GSE89632 (`{gse_id}_expression_matrix.csv`, genes x GSM
samples). 02_fetch_geo_cohorts.py's auto-build (gse.pivot_samples("VALUE"))
only works for series that publish a GSM-level VALUE field; both RNA-seq
cohorts here only publish supplementary count files, so this script does
that assembly by hand.

Both cohorts already report counts natively indexed by Ensembl gene ID
(ENSG...), which is Human-GEM's own gene namespace (model.genes[i].id is
ENSG, confirmed directly against models/Human-GEM.json) -- unlike GSE89632's
Illumina microarray data, no probe->Entrez->model-gene bridge is needed here.

Methodological choice (recorded per Mohammed's standing instruction to note
rather than pause on this kind of call): raw counts are converted to
counts-per-million (CPM) before saving, so samples are comparable to each
other and to GSE89632's already-normalized intensities when
group_mean_expression() averages across samples. Neither cohort has
duplicate gene rows needing a probe-style aggregation choice.

Requires: pandas   (already in gem1-main)
Run from gem1-main (no troppo/cobra dependency needed here).
"""

import gzip
import io
import re
import tarfile

import pandas as pd

DATA_DIR = "D:/mahdi/gem1/data/geo_cohorts"


def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    return counts.div(counts.sum(axis=0), axis=1) * 1e6


def build_gse126848():
    """
    GSE126848_Gene_counts_raw.txt.gz is one combined table (genes x samples),
    but its columns are internal lab codes ("Sample_description" in the GEO
    SOFT record), not GSM accessions. Build the code->GSM mapping from the
    family.soft.gz file, then rename columns to GSM so the matrix lines up
    with `sample_metadata.csv`'s `sample` column (as 03_build_context_specific_models.py
    expects).
    """
    gse_id = "GSE126848"
    soft_path = f"{DATA_DIR}/{gse_id}_family.soft.gz"
    counts_path = f"{DATA_DIR}/{gse_id}_Gene_counts_raw.txt.gz"

    gsm = None
    code_to_gsm = {}
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            m = re.match(r"\^SAMPLE = (GSM\d+)", line)
            if m:
                gsm = m.group(1)
                continue
            m = re.match(r"!Sample_description = (\S+)", line)
            if m and gsm:
                code_to_gsm[m.group(1)] = gsm

    # Sample codes are numeric lab IDs on both sides, but the counts file
    # zero-pads them to 4 digits (e.g. "0869") while the SOFT record's
    # Sample_description doesn't (e.g. "869") -- match as int, not string.
    code_to_gsm = {int(code): gsm for code, gsm in code_to_gsm.items()}

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.columns = [int(c) for c in counts.columns]
    unmapped = [c for c in counts.columns if c not in code_to_gsm]
    if unmapped:
        raise ValueError(f"{gse_id}: no GSM mapping found for sample codes {unmapped}")
    counts = counts.rename(columns=code_to_gsm)

    expr = cpm(counts)
    out_path = f"{DATA_DIR}/{gse_id}_expression_matrix.csv"
    expr.to_csv(out_path)
    print(f"{gse_id}: {expr.shape[0]} genes x {expr.shape[1]} samples (CPM) -> {out_path}")


def build_gse135251():
    """
    GSE135251_RAW.tar contains one HTSeq-count-style file per sample, named
    GSM<id>_<label>.counts.txt.gz (gene<TAB>count, no header, plus 5 trailing
    HTSeq summary rows __no_feature/__ambiguous/__too_low_aQual/__not_aligned/
    __alignment_not_unique that aren't genes and must be dropped). The GSM id
    is embedded directly in the filename, so no separate id-mapping step is
    needed here.
    """
    gse_id = "GSE135251"
    tar_path = f"{DATA_DIR}/{gse_id}_RAW.tar"
    summary_rows = {
        "__no_feature", "__ambiguous", "__too_low_aQual",
        "__not_aligned", "__alignment_not_unique",
    }

    series = {}
    with tarfile.open(tar_path) as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".counts.txt.gz")]
        for member in members:
            gsm_match = re.match(r"(GSM\d+)_", member.name)
            if not gsm_match:
                raise ValueError(f"{gse_id}: could not parse GSM id from {member.name}")
            gsm = gsm_match.group(1)
            raw = tar.extractfile(member).read()
            with gzip.open(io.BytesIO(raw), "rt") as f:
                s = pd.read_csv(f, sep="\t", header=None, index_col=0).iloc[:, 0]
            s = s[~s.index.isin(summary_rows)]
            series[gsm] = s

    counts = pd.DataFrame(series)
    expr = cpm(counts)
    out_path = f"{DATA_DIR}/{gse_id}_expression_matrix.csv"
    expr.to_csv(out_path)
    print(f"{gse_id}: {expr.shape[0]} genes x {expr.shape[1]} samples (CPM) -> {out_path}")


if __name__ == "__main__":
    build_gse126848()
    build_gse135251()
    print("\nBoth RNA-seq cohorts now have {gse_id}_expression_matrix.csv, same "
          "format as GSE89632's. Next: run 03_build_context_specific_models.py.")
