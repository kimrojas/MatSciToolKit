from matscitoolkit import gnu_parser


def test_read_gnu_datafile__np():
    d = gnu_parser.read_gnu_datafile(
        "test/test_data/bandstructure_np/si_bands.dat.gnu"
    )
    assert d.shape == (8, 91, 2)