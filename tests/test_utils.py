"""
Enhancements:
- Alot of string values are copied through out the utils code. Perhaps they can be added to the config as constants
"""
import pytest
from profile_dists import utils
import pandas as pd
import numpy as np
from numba.typed import Dict



@pytest.mark.parametrize("test_input,expected", [
    (["f5g1", "gg21", "hhhh"], 'hash'),
    (["1234", "1234", "3333"], 'int'),
    (["f5", "gg21", "123"], 'mix'),
    (["14", "1000", "12"], 'int'),
    ([14, 1000, 12], 'int'),
    ([""], ''),
    ([14, 20, "AA"], 'hash'),
    #! These tests are commented out for the sake of having some 
    #! set of tests be present and readily testable but if 
    #! further development is to be performed on Profile Dists
    #! these tests should be addressed.
    #! (["0.00", "1.0", "3.0"], 'mix'),
    #! (["0.3", 1.0, "4.2"], 'mix'),
    #! ([0.00, 1.0, 3.0], 'mix'),
    #! ([0.00, 1.1, 3.0], 'mix'),
    #! ([0.00, 1, 3.0], 'mix'),
    #! ["0.00", "1", "3.0"], 'mix'),
])
def test_guess_format(test_input, expected):
    """
    Considerations:
    - No base case is handled, function ends with elif statement, the returned value will 
    then be a '' string, which may not be a handled case
    - Are floating point values prevented from entering the function?
    """
    assert expected == utils.guess_format(test_input)



@pytest.mark.parametrize("test_input,expected", [
    (["1", "123", "1234"], False),
    ([1, 123, 1234], False),
    (["1", "2", "3"], True),
    (["", "", ""], True),
    ([""], True),
    (["0.00", "1.0", "3.0"], False),
    ([0.00, 1.0, 3.0], True),
])
def test_is_all_same_len(test_input, expected):
    assert expected == utils.is_all_same_len(test_input)


@pytest.mark.parametrize("test_input,expected", [
    #! These tests are commented out for the sake of completion but the need to be addressed
    #! (["1", "40000000000000000000000000000000000000000000000000000", "4"], True),
    #! ([1, 123, 400000000000000000000000000000000000000000000000000000000000], True),
    #! in future releases
    #! (["123", "abc"], False),
    #! (["fx00"], False),
    #! (["123ABC"], False),
    #! ([123, "this_is_a_string"], False),
    ([123, 123], True),
])
def test_contains_integers(test_input, expected):
    """
    Considerations:
    - The function exits before all values are tested, which could lead to incorrect inputs
    - Clarification needed if this function is meant check for any integer in a string, or if a list of values contains integers
    Enhancements:
    - Instead of regex, the '.is_digit()' method can be called
    """
    assert utils.contains_integers(test_input) == expected


@pytest.mark.parametrize("test_input,expected", [
    (["1", "40000000000000000000000000000000000000000000000000000", "4"], False),
    ([1, 123, 400000000000000000000000000000000000000000000000000000000000], False),
    (["123", "abc"], True),
    (["fx00"], True),
    (["123ABC"], True),
    ([123, "this_is_a_string"], True),
])
def test_contains_alpha(test_input, expected):
    """
    Tests if any of the list values contain a character
    """
    assert utils.contains_alpha(test_input) == expected



def test_convert_allele_codes():
    """Need to review function call to clarify input before test is written
    """
    ...

def test_update_column_map():
    """Requires col map
    """
    ...


@pytest.mark.parametrize("test_input,expected", [
    #! Commented out tests are valid, but these tests need to be addressed in future releases
    #! ({'1': [1], '2': [1], '3': [3]}, True),
    #! ({'1': [np.uint(1)], '2': [np.uint(1)], '3': [np.uint(3)]}, True),
    #! ({'1': [np.longlong(1)], '2': [np.ulonglong(1)], '3': [np.intc(3)], '4': [np.short(1)]}, True),
    #! ({'1': [1.0], '2': [1], '3': [3]}, True),
    ({'1': ["1.0"], '2': [1], '3': ["3"]}, False),
    ({'1': ["1"], '2': [1], '3': ["3"]}, False),
    ({'1': ["a"], '2': ["b"], '3': ["c"]}, False),
])
def test_is_all_columns_int(test_input, expected):
    """

    enhancements:
    - if VALID_INT_TYPES is not being updated, and only presence or absence in the set is being tested for the collection should be cast as a frozenset to prevent other updates and increse performance
    """
    assert utils.is_all_columns_int(pd.DataFrame(data=test_input).dtypes) == expected

@pytest.mark.parametrize("test_input,expected", [
    ({'1': [0, 1, 2, 4], '2': [0, 1, 0, 3], '3': [1, 2, 3, 4]}, 
        {'1': 1, '2': 2, '3': 0}),
])
def test_count_missing_data(test_input, expected):
    """
    Note: determines missing vals per column
    """
    assert expected == utils.count_missing_data(pd.DataFrame(data=test_input))


def create_typed_dict(dict_):
    """
    Create a typed.Dict for testing of numba programs
    :param dict_: Dictionary of values to be converted intoa dictionary for numba
    """
    td = Dict()
    for k, v in dict_.items():
        td[k] = v
    return td

@pytest.mark.parametrize("test_input,threshold,expected", [
    ({'1': 1, '2': 2, '3': 0}, 1, ['2']),
    ({'1': 1, '2': 2, '3': 0}, 2, []),
    ({'1': 1, '2': 2, '3': 0}, 4, []),
    ({'1': 1, '2': 2, '3': 0}, 0.1, ['1', '2']),
    pytest.param({'1': 1, '2': 2, '3': 0.11}, 0.1, ['1', '2', '3'], marks=pytest.mark.xfail(reason="numba.core.errors.TypingError: Failed in nopython mode pipeline")),
    ])
def test_identify_cols_to_remove(test_input, threshold, expected):
    """
    Consideration:
    - should overflow exception be handled
    
    Note:
    Test case below is removed, as it throws an OverflowError which should be handled somewhere else in the code base.
        ({'1': 1, '2': 2, '3': 100000000000000000000000000000000000000000}, 3, ['3'], OverflowError), 
        ({'1': 1, '2': 2, '3': 0.11011111111111111111111111111111111111}, 0.1, ['1', '2'], None),
    """
    assert expected == utils.identify_cols_to_remove(create_typed_dict(test_input), threshold)


def test_process_profile():
    """Requires staging test data
    """
    ...



@pytest.mark.parametrize("test_input,expected", [
    ({'1': [1, 2, 3, 4], '2': [1, 2, 3, 4], '3': [1, 2, 3, 4]}, 
    ([0, 1, 2, 3], [np.array([1, 1, 1]), np.array([2, 2, 2]), np.array([3,3,3]), np.array([4,4,4])])),
])
def test_convert_profiles(test_input, expected):
    """
    Input on this testcase is required to make it more exhaustive
    """
    
    row_ids, profiles = utils.convert_profiles(pd.DataFrame(data=test_input))
    assert row_ids == expected[0]
    assert all([np.array_equal(t, v) for t, v in zip(profiles, expected[1])])



@pytest.mark.parametrize("test_input,expected,func", [
    ((np.array([1, 1, 2]), np.array([1, 1, 3])), 1,   utils.get_distance_raw),
    ((np.array([1, 1, 0.1]), np.array([1, 1, 3])), 1, utils.get_distance_raw),
    ((np.array([1, 1, 2]), np.array([1, 1, 2])), 0,   utils.get_distance_raw),
    ((np.array([3, 2, 1]), np.array([1, 2, 3])), 2,   utils.get_distance_raw),
    ((np.array([4, 5, 6]), np.array([1, 2, 3])), 3,   utils.get_distance_raw),
    ((np.array([1, 1, 0]), np.array([1, 1, 0])), 0,   utils.get_distance_raw),
    
    ((np.array([1, 1, 2]), np.array([1, 1, 3])), 1,   utils.get_distance_raw_missing),
    ((np.array([1, 1, 0.1, 0.2]), np.array([1, 1, 3, 0.2])), 1, utils.get_distance_raw_missing),
    ((np.array([1, 1, 2]), np.array([1, 1, 2])), 0,   utils.get_distance_raw_missing),
    ((np.array([3, 2, 1]), np.array([1, 2, 3])), 2,   utils.get_distance_raw_missing),
    ((np.array([4, 5, 6]), np.array([1, 0, 0])), 3,   utils.get_distance_raw_missing),
    ((np.array([1, 1, 0]), np.array([1, 1, 0])), 0,   utils.get_distance_raw_missing),

    ((np.array([1, 1, 2]), np.array([1, 1, 3])), pytest.approx(float((3-2)/3) * 100, 0.01), utils.get_distance_scaled),
    ((np.array([1, 1, 0.1]), np.array([1, 1, 3])), pytest.approx(float((3-2)/3) * 100, 0.01), utils.get_distance_scaled),
    ((np.array([1, 1, 2]), np.array([1, 1, 2])), float(0), utils.get_distance_scaled),
    ((np.array([3, 2, 1]), np.array([1, 2, 3])), pytest.approx(float((3-1)/3) * 100, 0.01), utils.get_distance_scaled),
    ((np.array([4, 5, 6]), np.array([1, 2, 3])), pytest.approx(float((3-0)/3) * 100, 0.01), utils.get_distance_scaled),
    ((np.array([1, 1, 0]), np.array([1, 1, 0])), pytest.approx(float((2-2)/2) * 100, 0.01), utils.get_distance_scaled),

    ((np.array([1, 1, 2]), np.array([1, 1, 3])), pytest.approx(float((3-2)/3) * 100, 0.01), utils.get_distance_scaled_missing),
    ((np.array([1, 1, 0.1]), np.array([1, 1, 3])), pytest.approx(float((3-2)/3) * 100, 0.01), utils.get_distance_scaled_missing),
    ((np.array([1, 1, 2]), np.array([1, 1, 2])), float(0), utils.get_distance_scaled_missing),
    ((np.array([3, 2, 1]), np.array([1, 2, 3])), pytest.approx(float((3-1)/3) * 100, 0.01), utils.get_distance_scaled_missing),
    ((np.array([4, 5, 6]), np.array([1, 2, 3])), pytest.approx(float((3-0)/3) * 100, 0.01), utils.get_distance_scaled_missing),
    ((np.array([1, 1, 0]), np.array([1, 1, 0])), pytest.approx(float((3-3)/2) * 100, 0.01), utils.get_distance_scaled_missing),

])
def test_get_dist_funcs(test_input, expected, func):
    """
    Enhancements:
    - You may be able to use some numpy built-ins to speed these functions up. e.g. numpy.hamming, or by applying a function to the array
    """
    assert expected == func(test_input[0], test_input[1]), f"Function tested: {func.__name__}"


def test_calc_batch_size():
    ...


def test_validate_file():
    """
    Consideration:
    - Is the numba compilation needed here? if so make sure make sure the function is only compiled once.
    """
    ...

def test_compare_headers():
    """
    """
    ...

def test_guess_profile_format():
    """
    """
    ...

def test_get_file_length():
    """
    Requires setting up temp files
    """
    ...


###############################################################
"""
Enhancements
- Within the codebase these functions could likely all be converted into one, where the distance function is passed in
"""
def test_calc_distances_scaled():
    """
    Requires usage of test data
    """
    ...

def test_calc_distances_scaled_missing():
    """
    Requires usage of test data
    """
    ...


def test_calc_distances_hamming():
    """
    Requires usage of test data
    """
    ...

def test_calc_distances_hamming_missing():
    """
    Requires usage of test data
    """
    ...
########################################################

def test_if_file_ok():
    """
    Requires test data
    """
    ...


@pytest.mark.parametrize("labels,distances,threshold,expected,equivalent", [
    (["1", "2", "3"], [1, 2, 3], 1, {"1": 1}, True),
    (["1", "2", "3"], [1, 2, 3], 1, {"1": 1, "2": 2, "3": 3}, False),
    (["1", "2", "3"], [4, 5, 6], 1, {}, True),
    (["1", "2", "3"], [1, 2, 3], 0, {"1": 1, "2": 2, "3": 3}, False),
    (["1", "2", "3"], [1, 2, 3], 1.1, {"1": 1}, True),

])
def test_filter_dists(labels, distances, threshold, expected, equivalent):
    """
    Note:
    - Some tests here may be unfair and can be removed if they are handled else where in the code. e.g. those using numpy types
    
    Some tests were removed:
    (["1", "2", "3"], [1.1, 2, 3], np.int64(100), {"1": 1.1, "2": 2, "3": 3}, True),
    (["1", "2", "3"], [1.1, 2, 3], np.float64(100), {"1": 1.1, "2": 2, "3": 3}, True),
    
    """
    if equivalent:
        assert utils.filter_dists(labels, distances, threshold) == expected
    else:
        assert utils.filter_dists(labels, distances, threshold) != expected



def test_fromat_pairwise_dist():
    """
    Requires test data provided
    """
    ...

def test_write_dist_results():
    """
    Requires input data
    """
    ...


@pytest.mark.parametrize("profiles,labels,count_loci,expected", [
    ([np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])], ["1", "2", "3"], 3, {"1":  0.00, "2": 0.00, "3": 0.00}),
    #([np.array([0, 2, 3]), np.array([0, 5, 6]), np.array([0, 8, 9])], ["1", "2", "3"], 3, {"1":  0.33, "2": 0.33, "3": 0.33}),
    ([np.array([0, 2, 3]), np.array([0, 5, 6]), np.array([0, 8, 9])], ["1", "2", "3"], 3, {"1":  33.33, "2": 33.33, "3": 33.33}),
    ([np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])], ["1", "2", "3"], 4, {"1":  0.00, "2": 0.00, "3": 0.00}),
    ([np.array([0, 2, 3]), np.array([0, 5, 6]), np.array([0, 8, 9])], ["1", "2", "3"], 4, {"1":  25.00, "2": 25.00, "3": 25.00}),
    #([np.array([0, 2, 3]), np.array([0, 5, 6]), np.array([0, 8, 9])], ["1", "2", "3"], 4, {"1":  0.33, "2": 0.33, "3": 0.33}),
])
def test_get_missing_loci_counts(profiles,labels,count_loci,expected):
    """
    Enhancements:
    - Do the labels need to be passed here? or can the length of the input keys just be used?

    Note:
    - the docs say this function should return a % missing value, but it appears only the decimal value is returned
    - Is count_loci controlled to be the same for every value?
    """
    output = utils.get_missing_loci_counts(profiles, labels, count_loci)

    #! I am adding a multiplier of 100 to convert the output of the function into a percentage
    #! as it seems like this conversion is being handled somewhere else in the program
    percent_conversion = 100
    for k in labels:
        assert pytest.approx(output[k] * percent_conversion, 0.1) == expected[k]


@pytest.mark.parametrize("missing_counts,threshold,expected", [
    ({"1": 100.0, "2": 50.00, "3": 0.00}, 20.00, ["1", "2"]),
    ({"1": 100.0, "2": 50.00, "3": 0.00}, 100.00, []),
    ({"1": 40, "2": 0, "3": 1}, 20.00, ["1"]),
])
def test_flag_samples(missing_counts, threshold, expected):
    """
    """
    assert utils.flag_samples(missing_counts, threshold) == expected


@pytest.mark.parametrize("labels,profiles,labels_remove,expected_labels,expected_profiles", [
    (["1", "2", "3"], [np.array([1, 2, 3]), np.array([1, 2, 3]), np.array([1, 2, 3])], ["1", "2", "3"], [], []),
    (["1", "2", "3"], [np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])], ["1"], ["2", "3"], [np.array([4, 5, 6]), np.array([7, 8, 9])]),
])
def test_filter_samples(labels, profiles, labels_remove, expected_labels, expected_profiles):
    """
    """
    labels, profiles = utils.filter_samples(labels, profiles, labels_remove)
    assert expected_labels == labels
    for k, v in  zip(expected_profiles, profiles):
        assert np.alltrue(k == v)