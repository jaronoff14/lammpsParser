def test_import():
    import importlib
    mod = importlib.import_module("lammpsParser")
    assert hasattr(mod, "parse_frame_to_numpy")
