import pytest
import prism.util as util


def test_is_fully_methylated_or_unmethylated():
    assert util.is_fully_methylated_or_unmethylated('1111') == True
    assert util.is_fully_methylated_or_unmethylated('0000') == True
    assert util.is_fully_methylated_or_unmethylated('0111') == False
    assert util.is_fully_methylated_or_unmethylated('0001') == False

def test_jaccard_similarity():
    s1 = set([1, 2, 3, 4])
    s2 = set([2, 3, 4, 5])
    assert util.jaccard_similarity(s1, s2) == 0.6

def test_region_fail_because_of_invalid_region():
    with pytest.raises(ValueError):
        util.Region('chr1', 20, 10)
    
    with pytest.raises(ValueError):
        util.Region('chr1', 10, 10)

def test_region_string():
    r = util.Region('chr1', 10, 20)
    assert str(r) == 'chr1:10-20'

def test_extract_chromosome():
    assert util.extract_chromosome('chr1:100000;chr1:100002;chr1:100004') == 'chr1'

def test_merge_two_headers():
    h1 = 'chr1:100000;chr1:100002;chr1:100004;chr1:100010'
    h2 = 'chr:100002;chr1:100004;chr1:100010;chr1:100014;chr1:100018'
    assert util.merge_two_headers(h1, h2) == 'chr1:100000;chr1:100018'

def test_merge_two_headers_when_contained():
    h1 = 'chr1:100000;chr1:100002;chr1:100004;chr1:100010'
    h2 = 'chr:99998;chr1:100000;chr1:100004;chr1:100010;chr1:100014;chr1:100018'
    assert util.merge_two_headers(h1, h2) == 'chr1:99998;chr1:100018'

def test_preset_rc():
    util.preset_rc()
    util.preset_rc(font_family='FreeSans')
