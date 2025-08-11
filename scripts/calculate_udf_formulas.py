#!/usr/bin/env python

import logging
import re
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims

from scilifelab_epps.utils.formula import (
    fmol_to_ng,
    ng_to_fmol,
    ng_ul_to_nM,
    nM_to_ng_ul,
)
from scilifelab_epps.utils.udf_tools import fetch_last, get_art_tuples
from scilifelab_epps.wrapper import epp_decorator

DESC = """Script to perform calculations on input/output/step fields
by reading equations with special syntax from a step field.

The whole idea is to have a single calculation script whose behavior can be
customized on a step-by-step basis in the LIMS front-end configuration.

The equations are herein referred to as "formulas" and the input/output/step fields
are referenced using special syntax "placeholders".


Simple example:

    # This formula will calculate the field 'Sample volume (ul)' in the outputs of a step
    # It uses three different UDF placeholders

    outp['Diluted volume (ul)'] = outp['Sample volume (ul)'] + outp['Buffer volume (ul)']


Defintions:

    - UDF placeholders

        A special string referencing a UDF name in an input/output/step, and whether to fetch
        it recursively (indicated by a leading underscore). Can also be supplied as a
        prioritized list of UDF names.

        E.g.
            inp['foo']      step input artifact UDF 'foo'
            _outp['bar']    artifact UDF 'bar', fetched recursively from step output artifact
            step['A', 'B']  step UDF 'A', or 'B' if 'A' is not found

    - Formulas

        A string containing an equation with UDF placeholders.

        The left-hand side needs to be an isolated UDF placeholder, which is the one to
        be assigned, and the right-hand side will be evaluated using eval() after
        replacing UDF placeholders with their corresponding values.

        The separator between the left and right hand side can be either '=' or '==',
        where the former will overwrite existing UDFs and the latter will not.
        The latter case is useful for making sure multiple formulas do not overwrite
        one another, i.e. the first calculation to write to the field will be
        the one to stick.

        Right-hand side placeholders that can't be resolved will result in skipping the
        calculation. The silent skipping allows us to have multiple formulas targeting
        the same UDF that can be run in priority order without raising errors or
        warnings.

        Formulas are only allowed to contain
        - UDF placeholders
        - pure math, i.e. numbers, operators and parentheses
        - allowed functions
            - ng_ul, nM: For enforcing a unit-ambiguous concentration to be in
                ng/ul or nM, respectively.
            - fmol_to_ng, ng_to_fmol, ng_ul_to_nM, nM_to_ng_ul
                For converting between concentration units
        - allowed strings (e.g. 'ng/ul', 'nM')

Complex examples:

    1) Calculate the ng and fmol amount from a given concentration, volume and size of
       an output artifact.

        outp['Amount (ng)'] = ng_ul(outp['Concentration'], outp['Conc. Units'], outp['Size (bp)']) * outp['Volume (ul)']
        outp['Amount (fmol)'] = ng_to_fmol(outp['Amount (ng)'], outp['Size (bp)'])


    2) Prioritized calculations! Populate three UDFs, based on the supplied value of one
       of them, in priority order. Useful for having flexibility in what to base calculations on.
       It is done by having multiple non-overwriting (==) formulas trying to write to the same field
       from different inputs. The first one to resolve will be the one to stick.

        # Calculate ng amount and volume from given fmol amount
        outp['Amount for prep (ng)'] == fmol_to_ng(outp['Amount for prep (fmol)'], _outp['Size (bp)'])
        outp['Volume to take (uL)'] == fmol_to_ng(outp['Amount for prep (fmol)'], _outp['Size (bp)']) / ng_ul(inp['Concentration'], inp['Conc. Units'], _outp['Size (bp)'])

        # Calculate fmol amount and volume from given ng amount
        outp['Amount for prep (fmol)'] == ng_to_fmol(outp['Amount for prep (ng)'], _outp['Size (bp)'])
        outp['Volume to take (uL)'] == ng_to_fmol(outp['Amount for prep (ng)'], _outp['Size (bp)']) / nM(inp['Concentration'], inp['Conc. Units'], _outp['Size (bp)'])

        # Calculate fmol amount and ng amount given volume
        outp['Amount for prep (fmol)'] == outp['Volume to take (uL)'] * nM(inp['Concentration'], inp['Conc. Units'], _outp['Size (bp)'])
        outp['Amount for prep (ng)'] == outp['Volume to take (uL)'] * ng_ul(inp['Concentration'], inp['Conc. Units'], _outp['Size (bp)'])


    3) Prioritized UDFs! Fallback to using a different field, if the first one is not found.

        # Calculate amounts from input amounts by volume fraction
        outp['Amount (ng)'] = ( outp['Volume to take (uL)'] / inp['Volume (ul)', 'Total Volume (uL)'] ) * inp['Amount (ng)', 'Amount for prep (ng)']
        outp['Amount (fmol)'] = ( outp['Volume to take (uL)'] / inp['Volume (ul)', 'Total Volume (uL)'] ) * inp['Amount (fmol)', 'Amount for prep (fmol)']

        # Calculate new conc. from new volume
        outp['Concentration'] = outp['Amount (ng)'] / outp['Volume (ul)']
        outp['Conc. Units'] = 'ng/ul'
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


def ng_ul(conc: float, conc_units: str, size: float | None = None) -> float:
    """Force a concentration to be in ng/ul."""
    if conc_units == "ng/ul":
        return conc
    elif conc_units == "nM":
        assert size, "Size not provided for conversion to ng"
        return nM_to_ng_ul(conc, size)
    else:
        raise AssertionError(f"Concentration units '{conc_units}' not recognized")


def nM(conc: float, conc_units: str, size: float | None = None) -> float:
    """Force a concentration to be in nM."""
    if conc_units == "nM":
        return conc
    elif conc_units == "ng/ul":
        assert size, "Size not provided for conversion to nM"
        return ng_ul_to_nM(conc, size)
    else:
        raise AssertionError(f"Concentration units '{conc_units}' not recognized")


def assign_val_to_placeholder(
    val: str | float,
    placeholder: str,
    art_in: Artifact | None = None,
    art_out: Artifact | None = None,
    step: Process | None = None,
    overwrite: bool = True,
) -> None:
    """Assign a value to a UDF placeholder."""
    try:
        udf_name = re.search(r"\['(.*?)'\]", placeholder).groups()[0]  # type: ignore
    except (AttributeError, IndexError):
        raise AssertionError(
            f"Could not extract UDF name from placeholder '{placeholder}'"
        )

    # Where to assign UDF
    if "inp" in placeholder:
        assert art_in, "Input artifact not provided"
        obj = art_in
    elif "outp" in placeholder:
        assert art_out, "Output artifact not provided"
        obj = art_out
    elif "step" in placeholder:
        assert step, "Step not provided"
        obj = step

    if obj.udf.get(udf_name) and not overwrite:
        logging.info(
            f"UDF {placeholder} already exists and overwrite is False. Skipping calculation."
        )
        raise SkipCalculation()
    else:
        obj.udf[udf_name] = val
        obj.put()


def get_val_from_placeholder(
    placeholder: str,
    art_in: Artifact | None = None,
    art_out: Artifact | None = None,
    step: Process | None = None,
) -> str | float | int:
    """Fetch a value from a UDF placeholder."""
    recursive = True if placeholder[0] == "_" else False

    # Extract UDF names
    udf_names = re.findall(r"'(.*?)'", placeholder)
    assert udf_names, f"Could not extract UDF names from placeholder '{placeholder}'"

    # Where to fetch UDF
    if "inp" in placeholder:
        assert art_in, "Input artifact not provided"
        obj = art_in
        obj_type = "input artifact"
    elif "outp" in placeholder:
        assert art_out, "Output artifact not provided"
        obj = art_out
        obj_type = "output artifact"
    elif "step" in placeholder:
        assert step, "Step not provided"
        obj = step
        obj_type = "step"

    # Iterate across UDFs
    for i, udf_name in enumerate(udf_names):
        # How to fetch UDF
        if recursive:
            assert type(obj) is Artifact, (
                "Recursive UDF references only allowed for artifacts"
            )
            val = fetch_last(obj, udf_name, include_current=True, on_fail=None)
        else:
            val = obj.udf.get(udf_name)

        if val is None:
            # Re-try next UDF
            if i + 1 < len(udf_names):
                continue
            # All UDFs tried
            else:
                udf_names_quoted = [f"'{i}'" for i in udf_names]
                msg = (
                    "Could not resolve any of "
                    + f"UDFs {', '.join(udf_names_quoted)} "
                    + f"for {obj_type} '{obj.type.name if 'step' in placeholder else obj.name}' "
                    + f"({obj.id})."
                )
                logging.info(msg)
                val = None
        else:
            break

    # Returned values will pass through eval, so strings need
    # to be escaped to not be interpreted as variables
    if isinstance(val, str):
        val = f"'{val}'"

    if len(udf_names) > 1:
        # Only clarify which UDF was used if multiple ones were provided
        logging.info(f"Resolved {placeholder} UDF '{udf_name}' to {val}")

    assert isinstance(val, (int, str, float)) or val is None, (
        f"Unexpected type for val: {type(val)}"
    )

    return val


class SkipCalculation(Exception):
    """Custom exception for skipping calculations."""

    def __init__(self):
        super().__init__()


def parse_formula(formula: str) -> tuple[str, list[str]]:
    """From a formula string, extract UDF placeholders.

    Return the formula string with placeholders replaced by curly brackets,
    as well as a list of the placeholders.

    Also perform various assertions to sanity check the input and prevent code injection.
    """
    logging.info(f"Parsing formula: {formula}")

    # Explicate non-placeholders allowed in formula, to prevent code injection
    allowed_functions = [
        ng_ul,
        nM,
        fmol_to_ng,
        ng_to_fmol,
        ng_ul_to_nM,
        nM_to_ng_ul,
    ]
    allowed_strings = [
        "'ng/ul'",
        "'nM'",
    ]

    # Patterns
    placeholder_pattern = r"_?((inp)|(outp)|(step))\[.*?\]"
    allowed_functions_pattern = "|".join([f"({f.__name__})" for f in allowed_functions])
    allowed_strings_pattern = "|".join([f"({s})" for s in allowed_strings])
    pure_math_pattern = r"[=\d,\+\-\*\/\(\)\{\}\s]+"

    # Collect UDF references from formula
    placeholders: list[str] = [
        match[0] for match in re.finditer(placeholder_pattern, formula)
    ]

    # Replace UDF references with curly brackets
    formula_fstring: str = re.sub(placeholder_pattern, r"{}", formula)

    # Sanity check matching number of UDF references and placeholders
    assert formula_fstring.count(r"{}") == len(placeholders), (
        f"Number of extracted UDF references ({len(placeholders)})"
        + f"do not match number of format placeholders ({formula_fstring.count(r'{}}')})"
    )

    # Strip down formula pattern by pattern, to see if any disallowed content remains
    # This is to prevent code injection
    formula_stripped = formula_fstring
    formula_stripped = re.sub(allowed_functions_pattern, "", formula_stripped)
    formula_stripped = re.sub(allowed_strings_pattern, "", formula_stripped)
    formula_stripped = re.sub(pure_math_pattern, "", formula_stripped)

    # Assert contents are empty after stripping
    # If not raise the remainins for troubleshooting
    assert formula_stripped == "", (
        f"Formula '{formula}' appears to contain disallowed characters: '{formula_stripped}'"
    )

    return formula_fstring, placeholders


def eval_rh(
    formula_fstring: str,
    placeholders: list[str],
    art_in: Artifact | None = None,
    art_out: Artifact | None = None,
    step: Process | None = None,
) -> str | float:
    """Evaluate the right-hand side of the formula, with placeholders replaced by values."""
    # Translate right-hand placeholders to values
    rh_values = []
    for placeholder in placeholders[1:]:  # First placeholder to be assigned, not read
        rh_val = get_val_from_placeholder(placeholder, art_in, art_out, step)
        rh_values.append(rh_val)

    # Evaluate right-hand side of equation
    sep = "==" if "==" in formula_fstring else "="
    rh_fstring = formula_fstring.split(sep)[1].strip()
    rh_eval_string = rh_fstring.format(*rh_values)

    # Solve for x :)
    try:
        lh_val = eval(rh_eval_string)
    except Exception as e:
        logging.error(f'Could not evaluate: "{rh_eval_string}"')
        raise e
    assert type(lh_val) in [float, int, str], (
        f'Evaluation of "{rh_eval_string}" gave invalid output: "{lh_val}" of type "{type(lh_val)}"'
    )

    values = [lh_val] + rh_values
    values_2f = [
        f"{value:.2f}" if isinstance(value, float) else value for value in values
    ]

    # Print equations with placeholders and populated values
    logging.info(f"    Calculation:  {formula_fstring.format(*values_2f)}")

    return lh_val


def apply_formula(
    process: Process, step_type: str, formula_fstring: str, placeholders: list[str]
) -> None:
    """This function takes the parsed formula and applies it.

    The application will differ depending on the type of step.
    """

    overwrite = False if "==" in formula_fstring else True

    # Iterate across artifacts
    if step_type in ["Standard", "Add Labels"]:
        for art_tuple in get_art_tuples(process):
            art_in = art_tuple[0]["uri"]
            art_out = art_tuple[1]["uri"]
            logging.info(
                f"Calculations for input-output '{art_in.name}' ({art_in.id}) --> '{art_out.name}' ({art_out.id})"
            )
            try:
                val = eval_rh(
                    formula_fstring,
                    placeholders,
                    art_in=art_in,
                    art_out=art_out,
                    step=process,
                )
                assign_val_to_placeholder(
                    val, placeholders[0], art_in, art_out, process, overwrite
                )
            except SkipCalculation:
                continue
    elif step_type == "No Outputs":
        logging.info("Step type: No-output")
        for art_in in [i for i in process.all_inputs() if i.type == "Analyte"]:
            logging.info(f"Calculations for input '{art_in.name}' ({art_in.id})")
            try:
                val = eval_rh(
                    formula_fstring,
                    placeholders,
                    art_in=art_in,
                    step=process,
                )
                assign_val_to_placeholder(
                    val, placeholders[0], art_in, process, overwrite
                )
            except SkipCalculation:
                continue
    else:
        raise NotImplementedError(f"Step type '{step_type}' not implemented")


def get_formulas(step, formula_field):
    """Extract formulas from a step UDF text field. Skip empty lines and comments."""

    udf_names = [udf[0] for udf in step.udf.items()]
    formula_field_contents = [
        step.udf.get(udf_name) for udf_name in udf_names if formula_field in udf_name
    ]
    assert len(formula_field_contents) == 1, (
        f"Searching for '{formula_field}' in step UDFs gave {len(formula_field_contents)} matches."
    )
    rows = formula_field_contents[0].split("\n")

    formulas = []
    for row in rows:
        row = row.strip()
        if not row or row[0] == "#":
            continue
        formulas.append(row)

    return formulas


def get_step_type(process: Process) -> str:
    """Determine the type of step, based on the number of outputs."""

    ## Step type IDs
    # Standard / No Outputs  24-XXXXXX
    # Add Labels             151-XXXXXX
    # Pooling                122-XXXXXX

    step_prefix = process.id.split("-")[0]

    if step_prefix == "24":
        # This function should return None for 'No Outputs'-steps
        if get_art_tuples(process):
            step_type = "Standard"
        else:
            step_type = "No Outputs"
    elif step_prefix == "151":
        step_type = "Add Labels"
    else:
        raise AssertionError(f"Step prefix '{step_prefix}' not recognized")

    logging.info(f"Step type: {step_type}")
    return step_type


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    step_type = get_step_type(process)

    formulas = get_formulas(process, args.formula_field)
    for formula in formulas:
        formula_fstring, placeholders = parse_formula(formula)
        apply_formula(process, step_type, formula_fstring, placeholders)


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid",
        required=True,
        type=str,
        help="Lims ID for current Process",
    )
    parser.add_argument(
        "--log",
        required=True,
        type=str,
        help="Which log file slot to use",
    )
    parser.add_argument(
        "--formula_field",
        type=str,
        default="UDF formulas",
        help="String to identify the step field containing formulas",
    )
    args = parser.parse_args()

    main(args)
