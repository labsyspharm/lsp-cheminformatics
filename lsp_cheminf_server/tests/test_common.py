import pytest

import lspcheminf

test_compounds = {
    "inchi": [
        "InChI=1S/C23H31N5O4/c1-15(2)14-28-21(29)19-20(25(3)23(28)30)24-22-26(10-6-11-27(19)22)12-9-16-7-8-17(31-4)18(13-16)32-5/h7-8,13,15H,6,9-12,14H2,1-5H3",
        "InChI=1S/C31H30F3N3O3S/c1-19(2)26-13-12-25(16-20(26)3)37(17-21-4-6-23(7-5-21)29(40)35-15-14-28(38)39)30-36-27(18-41-30)22-8-10-24(11-9-22)31(32,33)34/h4-13,16,18-19H,14-15,17H2,1-3H3,(H,35,40)(H,38,39)",
        "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
        "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+",
    ],
    "smiles": [
        r"COc1ccc(CCN2CCCn3c2nc4N(C)C(=O)N(CC(C)C)C(=O)c34)cc1OC",
        r"CC(C)c1ccc(cc1C)N(Cc2ccc(cc2)C(=O)NCCC(=O)O)c3nc(cs3)c4ccc(cc4)C(F)(F)F",
        r"C[C@@H]1CCN(C[C@@H]1N(C)c2ncnc3[nH]ccc23)C(=O)CC#N",
        r"Oc1ccc(\C=C\c2cc(O)cc(O)c2)cc1",
    ],
}


@pytest.fixture
def client():
    lspcheminf.app.testing = True
    with lspcheminf.app.test_client() as client:
        yield client
