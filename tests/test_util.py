import prism.util as util


def test_is_fully_methylated_or_unmethylated():
    assert util.is_fully_methylated_or_unmethylated('1111') == True
    assert util.is_fully_methylated_or_unmethylated('0000') == True
    assert util.is_fully_methylated_or_unmethylated('0111') == False
    assert util.is_fully_methylated_or_unmethylated('0001') == False
