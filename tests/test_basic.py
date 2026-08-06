def test_import():
    import importlib
    mod = importlib.import_module("lammpsParser")
    assert hasattr(mod, "parse_frame_to_numpy")
    

def test_data():
    import lammpsParser
    
    file = 'dump.test.lammpstrj'
    lammpsParser.parse_frame_to_numpy(file)
    
    lammpsParser.parse_frame_to_numpy(file,frame_num=1)
    print("Finished")