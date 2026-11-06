"""
Collecting the requested Pfam profiles out of the master Pfam-A.hmm
"""

import os
import shutil

from gtotree.tests.paths import MOCK_PFAM_HMM, MOCK_PFAM_INFO
from gtotree.utils.misc.general import RunData
from gtotree.utils.pfam.get_pfam_data import HMM_FILENAME, INFO_FILENAME
from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets


def _collect(tmp_path, monkeypatch, targets):
    pfam_dir = tmp_path / "pfam-data"
    pfam_dir.mkdir()
    shutil.copy(MOCK_PFAM_HMM, pfam_dir / HMM_FILENAME)
    shutil.copy(MOCK_PFAM_INFO, pfam_dir / INFO_FILENAME)
    monkeypatch.setenv("Pfam_data_dir", str(pfam_dir))

    target_file = tmp_path / "targets.txt"
    target_file.write_text("\n".join(targets) + "\n")
    rd = RunData()
    rd.target_pfams_file = str(target_file)
    rd.pfam_results_dir = str(tmp_path / "pfam-results")
    os.makedirs(os.path.join(rd.pfam_results_dir, "info"))
    return get_additional_pfam_targets(rd)


def test_one_pfam_listed_under_two_spellings_is_pulled_once_and_neither_fails(
        tmp_path, monkeypatch, capsys):
    # both spellings used to get their own entry, only the last was matched, and the
    # other was reported as a failed target even though the profile was pulled
    rd = _collect(tmp_path, monkeypatch, ["PF90001", "PF90001.3", "PF90002"])

    assert rd.found_pfam_targets == ["PF90001", "PF90002"]
    assert rd.failed_pfam_targets == []
    assert rd.pfam_dict == {"PF90001": "PF90001.3", "PF90001.3": "PF90001.3",
                            "PF90002": "PF90002.7"}
    assert "same Pfam (PF90001)" in capsys.readouterr().out


def test_a_genuinely_missing_pfam_is_still_reported_as_failed(tmp_path, monkeypatch):
    rd = _collect(tmp_path, monkeypatch, ["PF90001", "PF99999"])
    assert rd.found_pfam_targets == ["PF90001"]
    assert rd.failed_pfam_targets == ["PF99999"]
