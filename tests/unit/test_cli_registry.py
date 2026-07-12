from phykit.cli_registry import (
    ALIAS_TO_HANDLER,
    COMMAND_IDENTITIES,
    PUBLIC_COMMAND_TO_HANDLER,
)


def test_every_handler_has_one_command_identity():
    handlers = {identity.handler for identity in COMMAND_IDENTITIES}
    assert handlers == set(ALIAS_TO_HANDLER.values())
    assert len(handlers) == len(COMMAND_IDENTITIES)


def test_public_command_names_are_unique_and_complete():
    expected = {
        name: identity.handler
        for identity in COMMAND_IDENTITIES
        for name in (identity.canonical, *identity.aliases)
    }
    assert PUBLIC_COMMAND_TO_HANDLER == expected


def test_canonical_entry_points_include_descriptive_legacy_commands():
    expected = {
        "degree_of_violation_of_a_molecular_clock": "dvmc",
        "long_branch_score": "lb_score",
        "relative_composition_variability": "rcv",
        "relative_composition_variability_taxon": "rcvt",
        "robinson_foulds_distance": "rf_distance",
    }
    for command, handler in expected.items():
        assert PUBLIC_COMMAND_TO_HANDLER[command] == handler


def test_every_canonical_command_has_a_standalone_entry_point():
    for identity in COMMAND_IDENTITIES:
        assert f"pk_{identity.canonical}" in identity.entry_points
