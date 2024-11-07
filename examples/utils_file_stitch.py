
import os


def test_file_stitch():
    import numpy as np

    from carmm.utils.file_stitch import file_stitch

    mock_data = np.ones(6)

    file_stitch(path='data/H+CH4_CH3+H2_path/H+CH4_CH3+H2.xyz',
                out_fname='data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz',
                data=mock_data,
                lines_per_image=8)
    with open('data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz', 'r') as test:
        lines = test.readlines()

    assert str(lines[0]) == '6\n'
    assert str(lines[1]) == 'Frame 1\n'
    assert str(lines[2]) == '   C   0.63304144726145573 -0.1087864399671382 -0.56719035628535563   1.0\n'
    assert str(lines[8]) == '6\n'
    assert str(lines[9]) == 'Frame 1\n'
    assert str(lines[10]) == '   C   0.88562473811384212 0.97416578530588571 1.0453399509875991   1.0\n'


test_file_stitch()

# Copy for reference
os.system('cp data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched_REF.xyz')

# Remove the test file
os.system('rm data/H+CH4_CH3+H2_path/H+CH4_CH3+H2_stitched.xyz')
