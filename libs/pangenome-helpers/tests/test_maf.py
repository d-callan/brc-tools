from pathlib import Path

import pytest

from pangenome_helpers.maf import derive_multiz_order, process_maf_file

DATA = Path(__file__).parent / "data"


def test_process_maf_file_filters_and_sorts(tmp_path):
    src = DATA / "maf" / "sample.maf"
    dst = tmp_path / "out.maf"
    result = process_maf_file("REF", src, dst)
    assert result.kept == 1
    assert result.dropped == 1
    text = dst.read_text().strip().splitlines()
    assert text[0].startswith("##maf")
    assert any(line.startswith("s REF.chr1") for line in text)


def test_process_maf_file_requires_ref(tmp_path):
    src = tmp_path / "bad.maf"
    src.write_text("##maf\n")
    dst = tmp_path / "out.maf"
    result = process_maf_file("REF", src, dst)
    assert result.kept == 0
    assert result.dropped == 0


def test_derive_multiz_order_sorted_desc():
    compare = DATA / "multiz" / "compare.csv"
    ordered = derive_multiz_order(compare, "hinge", ["query1", "query2"])
    assert ordered == ["query1", "query2"]


def test_derive_multiz_order_requires_queries():
    compare = DATA / "multiz" / "compare.csv"
    with pytest.raises(ValueError):
        derive_multiz_order(compare, "hinge", [])
