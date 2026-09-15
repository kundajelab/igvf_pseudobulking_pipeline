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
        ("already_clean_name", "already_clean_name"),
    ],
)
@pytest.mark.parametrize("sanitize_hyphen", [True, False])
def test_sanitize_to_ascii_underscore(text: str, expected: str, sanitize_hyphen: bool) -> None:
    assert sanitize_to_ascii_underscore(text, sanitize_hyphen=sanitize_hyphen) == expected


@pytest.mark.parametrize(
    ("text", "sanitize_hyphen", "expected"),
    [
        # by default "-" is replaced, since it separates the cell name from the subsample in
        # pseudobulk IDs
        ("already_clean-name", True, "already_clean_name"),
        ("CD4-positive, alpha-beta T cell", True, "CD4_positive_alpha_beta_T_cell"),
        ("a - b", True, "a_b"),
        ("already_clean-name", False, "already_clean-name"),
        ("CD4-positive, alpha-beta T cell", False, "CD4-positive_alpha-beta_T_cell"),
        ("a - b", False, "a_-_b"),
    ],
)
def test_sanitize_to_ascii_underscore_hyphens(
    text: str, sanitize_hyphen: bool, expected: str
) -> None:
    assert sanitize_to_ascii_underscore(text, sanitize_hyphen=sanitize_hyphen) == expected


def test_get_map_to_sanitized_keeps_collisions_unique() -> None:
    cell_names = pd.Series(["SC.alpha", "SC/alpha", "SC_alpha", "SC.beta", "SC.alpha"])
    mapping = _get_map_to_sanitized(cell_names)
    assert set(mapping) == {"SC.alpha", "SC/alpha", "SC_alpha", "SC.beta"}
    assert len(set(mapping.values())) == len(mapping)
    assert mapping["SC.beta"] == "SC_beta"
    assert all(not any(char in sanitized for char in "./\\, ") for sanitized in mapping.values())


@pytest.mark.parametrize("sanitize_hyphen", [True, False])
def test_get_map_to_sanitized_hyphens(sanitize_hyphen: bool) -> None:
    """Names that differ only by "-" or "_" collide only when hyphens are sanitized."""
    mapping = _get_map_to_sanitized(pd.Series(["a-b", "a_b"]), sanitize_hyphen=sanitize_hyphen)
    if sanitize_hyphen:
        assert mapping == {"a-b": "a_b_1", "a_b": "a_b_2"}
    else:
        assert mapping == {"a-b": "a-b", "a_b": "a_b"}
