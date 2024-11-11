
import os


def test_file_stitch_1toA():
    import numpy as np

    from carmm.utils.file_stitch import file_stitch

    # 2 frames of 6 atoms. 1 dataset of 6 values broadcast to all
    mock_data = np.ones(6)

    file_stitch(path='data/H+CH4_CH3+H2_path/H+CH4_CH3+H2.xyz',
                out_fname='data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz',
                data=mock_data,
                lines_per_image=8,
                mode='1toA')

    with open('data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz', 'r') as test:
        lines = test.readlines()

    assert str(lines[0]) == '6\n'
    assert str(lines[1]) == 'Frame 1\n'
    assert str(lines[2]) == '   C   0.63304144726145573 -0.1087864399671382 -0.56719035628535563   1.0\n'
    assert str(lines[8]) == '6\n'
    assert str(lines[9]) == 'Frame 1\n'
    assert str(lines[10]) == '   C   0.88562473811384212 0.97416578530588571 1.0453399509875991   1.0\n'


def test_file_stitch_AtoA():
    import numpy as np

    from carmm.utils.file_stitch import file_stitch

    # 2 fames of six atoms. Dataset of 12 values broadcast
    mock_data = np.ones((12,3))

    file_stitch(path='data/H+CH4_CH3+H2_path/H+CH4_CH3+H2.xyz',
                out_fname='data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_AtoA.xyz',
                data=mock_data,
                lines_per_image=8,
                mode='AtoA')

    with open('data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_AtoA.xyz', 'r') as test:
        lines = test.readlines()

    assert str(lines[0]) == '6\n'
    assert str(lines[1]) == 'Frame 1\n'
    assert str(lines[2]) == '   C   0.63304144726145573 -0.1087864399671382 -0.56719035628535563   1. 1. 1.\n'
    assert str(lines[8]) == '6\n'
    assert str(lines[9]) == 'Frame 1\n'
    assert str(lines[10]) == '   C   0.88562473811384212 0.97416578530588571 1.0453399509875991   1. 1. 1.\n'


test_file_stitch_1toA()

test_file_stitch_AtoA()

# Copy for reference
os.system('cp data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_REF.xyz')
os.system('cp data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_AtoA.xyz data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_AtoA_REF.xyz')

# Remove the test file
os.system('rm data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz')
os.system('rm data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_AtoA.xyz')
