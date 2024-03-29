from pathlib import Path
import shutil
import os
import sys


class TemporaryEnvironment:
    def __init__(self, *tmpname):
        # Reference to the main path (relative path handling)
        self.main_path = Path.cwd()

        # Test path handling
        self.test_path = self.main_path / Path("test")

        # Temporary path handling
        tmpname = "_".join([str(i) for i in tmpname])
        self.tmp_path = self.test_path / Path("tmp") / Path(tmpname)
        if self.tmp_path.exists():
            shutil.rmtree(self.tmp_path)
        
        self.tmp_path.mkdir(parents=True, exist_ok=True)
        os.chdir(self.tmp_path)

    def goto_tmp(self):
        os.chdir(self.tmp_path)

    def clean(self):
        shutil.rmtree(self.tmp_path)

    def return_to_main(self):
        os.chdir(self.main_path)

    def add(self, path):
        _name = Path(path).name
        link = Path(_name)
        
        link.symlink_to(path)