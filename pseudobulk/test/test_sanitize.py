import pandas as pd
import pytest

from pseudobulk.utils import _get_map_to_sanitized, sanitize_to_ascii_underscore


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("SC.alpha", "SC_alpha"),
        ("T cell, CD4+", "T_cell_CD4+"),
        ("CD4/CD8", "CD4_CD8"),
        ("CD4\\CD8", "CD4_CD8"),
        ("a\tb\nc", "a_b_c"),
        ("a. /\\,b", "a_b"),
        ("café", "cafe"),
        ("already_clean-name", "already_clean-name"),
    ],
)
def test_sanitize_to_ascii_underscore(text: str, expected: str) -> None:
    assert sanitize_to_ascii_underscore(text) == expected


def test_get_map_to_sanitized_keeps_collisions_unique() -> None:
    cell_names = pd.Series(["SC.alpha", "SC/alpha", "SC_alpha", "SC.beta", "SC.alpha"])
    mapping = _get_map_to_sanitized(cell_names)
    assert set(mapping) == {"SC.alpha", "SC/alpha", "SC_alpha", "SC.beta"}
    assert len(set(mapping.values())) == len(mapping)
    assert mapping["SC.beta"] == "SC_beta"
    assert all(not any(char in sanitized for char in "./\\, ") for sanitized in mapping.values())
