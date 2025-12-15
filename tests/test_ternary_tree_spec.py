import pytest
import random

from yaf2q.ternary_tree_spec import TernaryTreeSpec


def test_constructor():

    ttspec = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec.indices == [1, 0, 3, 2]
    assert ttspec.edges == {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')}
    assert ttspec.is_valid() is True


def test_constructor_too_large_indices():

    with pytest.raises(ValueError):
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2, 4],
            edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
        )
        ttspec


def test_constructor_duplicate_indices():

    with pytest.raises(ValueError):
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 1],
            edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
        )
        ttspec


def test_constructor_duplicate_parent():

    with pytest.raises(ValueError):
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2],
            edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Z')},
        )
        ttspec


def test_random():

    # seed = 1
    random.seed(1)
    ttspec = TernaryTreeSpec.random(4)
    assert ttspec.is_valid() is True

    # seed = 12
    random.seed(12)
    ttspec = TernaryTreeSpec.random(4)
    assert ttspec.is_valid() is True

    # seed = 123
    random.seed(123)
    ttspec = TernaryTreeSpec.random(4)
    assert ttspec.is_valid() is True


def test_to_string_rust():

    ttspec = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec._to_string_rust() == "1 0 3 2\n1 2\n2 3\n3 1"


def test_to_string():

    ttspec = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec.to_string() == "indices:[1, 0, 3, 2], edges:{1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')"


def test_swap_indices():

    ttspec = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )

    # seed = 1
    random.seed(1)
    ttspec_swapped = ttspec.swap_indices()
    ttspec_expect = TernaryTreeSpec(
        indices = [0, 1, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec_swapped == ttspec_expect

    # seed = 12
    random.seed(12)
    ttspec_swapped = ttspec.swap_indices()

    ttspec_expect = TernaryTreeSpec(
        indices = [1, 0, 2, 3],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec_swapped == ttspec_expect
    
    # seed = 123
    random.seed(123)
    ttspec_swapped = ttspec.swap_indices()
    ttspec_expect = TernaryTreeSpec(
        indices = [3, 0, 1, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec_swapped == ttspec_expect


def test_change_parents():

    ttspec = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
    )

    # seed = 1
    random.seed(1)
    ttspec_changed = ttspec.change_parent()
    ttspec_expect = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec_changed == ttspec_expect

    # seed = 12
    random.seed(12)
    ttspec_changed = ttspec.change_parent()
    ttspec_expect = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'Z'), 2: (1, 'Z'), 3: (0, 'Y')},
    )
    assert ttspec_changed == ttspec_expect
    
    # seed = 123
    random.seed(123)
    ttspec_changed = ttspec.change_parent()
    ttspec_expect = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert ttspec_changed == ttspec_expect
