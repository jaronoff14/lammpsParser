
def test_import():
    import importlib
    mod = importlib.import_module("lammpsParser")
    assert hasattr(mod, "parse_frame_to_numpy")
    
def test_dump():
    import lammpsParser
    file = '/Users/14ath/Documents/Coding/C/pybind/parser/lammpsParser/tests/dump.test.lammpstrj'
    lammpsParser.parse_frame_to_numpy(file)
    file1_dump= lammpsParser.parse_frame_to_numpy(file,frame_num=1)
    box_data = file1_dump["header"]
    box_bounds = box_data["box_bounds"]
    assert(file1_dump["n_atoms"] != 0)
    assert(box_bounds[0][1] !=0 and box_bounds[0][0] !=0)
    
    
    file2= '/Users/14ath/Documents/Coding/C/pybind/parser/lammpsParser/tests/dump.big.lammpstrj'
    file2_dump = lammpsParser.parse_frame_to_numpy(file2,frame_num=785)
    assert(file2_dump["n_atoms"] != 0)
    box_data = file2_dump["header"]
    box_bounds = box_data["box_bounds"]
    assert(box_bounds[0][1] !=0 and box_bounds[0][0] !=0)
    assert(abs(box_bounds[0][1]) + abs(box_bounds[0][0]) != 0  and abs(box_bounds[1][0])+ abs(box_bounds[1][1]) != 0)
    
def test_data():
    import lammpsParser
    
    
    
    
    
    data = '/Users/14ath/Documents/Coding/C/pybind/parser/lammpsParser/tests/datafile_test.data'
    atoms_data = lammpsParser.parse_data_to_numpy(data)
    box_data = atoms_data["header"]
    assert(box_data["xhi"] != 0 and box_data["xlo"] !=0)
    assert(box_data["yhi"] != 0 and box_data["ylo"] !=0)
    
