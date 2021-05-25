# -*- coding: utf-8 -*-
"""Process GTEx json data to xlsx."""
import json

import pandas as pd
from logs.get_logger import get_logger

INS = snakemake.input  # noqa: F821
LOGS = snakemake.log[0]  # noqa: F821
OUTS = snakemake.output  # noqa: F821

logger = get_logger(__name__, LOGS)

with pd.ExcelWriter(OUTS["data"], engine="openpyxl") as writer:
    for path in INS["data"]:
        with open(path, "r") as file:
            raw = json.load(file)
        data = pd.DataFrame.from_dict(raw["medianTranscriptExpression"])
        data = data.sort_values("median", ascending=False)
        gene = data["geneSymbol"].unique()[0]
        data.to_excel(writer, index=False, sheet_name=gene)
