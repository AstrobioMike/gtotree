"""
VeryFastTree is the default tree program (it was FastTreeMP): same algorithm, same
trees in testing, about a third of the time. `-T FastTreeMP` keeps the old behavior.
"""

import shutil

import pytest  # type: ignore

from gtotree.cli.parser import parser
from gtotree.tests.paths import DATA_DIR


def test_the_default_tree_program_is_veryfasttree():
    assert parser().parse_args(["-A", "x", "-H", "x"]).tree_program == "VeryFastTree"


def test_fasttreemp_can_still_be_chosen():
    assert parser().parse_args(["-A", "x", "-H", "x", "-T", "FastTreeMP"]).tree_program \
        == "FastTreeMP"


@pytest.mark.skipif(
    any(shutil.which(b) is None for b in ("muscle", "trimal", "VeryFastTree")),
    reason="preflight needs muscle, trimal and VeryFastTree")
def test_a_default_run_cites_veryfasttree_and_fasttree(tmp_path):
    from gtotree.utils.misc.preflight_checks import preflight_checks
    from gtotree.utils.misc.context import log_file_var
    token = log_file_var.set(log_file_var.get())
    try:
        _args, run_data = preflight_checks(parser().parse_args(
            ["-A", *(str(DATA_DIR / f"mock-{i}.faa") for i in range(1, 5)),
             "-H", str(DATA_DIR / "mock.hmm"), "-o", str(tmp_path / "out")]))
    finally:
        log_file_var.reset(token)
    # VeryFastTree implements FastTree's method, so both are cited
    assert run_data.tools_used.veryfasttree_used
    assert run_data.tools_used.fasttree_used
