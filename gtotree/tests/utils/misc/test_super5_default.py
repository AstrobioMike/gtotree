"""
muscle's 'super5' is the alignment default at any genome count; `-X` turns it off at
any genome count (it used to only stop an automatic switch above 500 genomes)
"""

import shutil

import pytest  # type: ignore

from gtotree.tests.paths import DATA_DIR

pytestmark = pytest.mark.skipif(
    any(shutil.which(b) is None for b in ("muscle", "trimal", "VeryFastTree")),
    reason="preflight needs muscle, trimal and VeryFastTree")


def _preflight(tmp_path, *extra):
    from gtotree.cli.parser import parser
    from gtotree.utils.misc.preflight_checks import preflight_checks
    from gtotree.utils.misc.context import log_file_var

    listing = tmp_path / "aa.txt"
    listing.write_text("\n".join(str(DATA_DIR / f"mock-{i}.faa") for i in range(1, 5)))
    argv = ["-A", str(listing), "-H", str(DATA_DIR / "mock.hmm"),
            "-o", str(tmp_path / "out"), *extra]
    token = log_file_var.set(log_file_var.get())
    try:
        return preflight_checks(parser().parse_args(argv))
    finally:
        log_file_var.reset(token)


def test_super5_is_used_by_default_even_with_few_genomes(tmp_path):
    _args, run_data = _preflight(tmp_path)
    assert run_data.use_muscle_super5


def test_no_super5_turns_it_off_at_any_genome_count(tmp_path):
    _args, run_data = _preflight(tmp_path, "-X")
    assert not run_data.use_muscle_super5


def test_muscle_is_called_with_the_chosen_algorithm(tmp_path, monkeypatch):
    import types
    from gtotree.utils.misc import seqs
    calls = []
    monkeypatch.setattr(seqs.subprocess, "run",
                        lambda cmd, **kw: calls.append(cmd))
    for flag, expected in ((True, "-super5"), (False, "-align")):
        rd = types.SimpleNamespace(use_muscle_super5=flag, num_muscle_threads=1)
        seqs.run_muscle("S", rd, "in.faa", "out.faa", str(tmp_path / "log"))
        assert calls[-1][1] == expected
