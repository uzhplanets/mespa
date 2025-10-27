"""
Author: Henrik Knierim
Script: install.py
Description: This script installs element_sedimentation.
Date: 2025-05-28
"""

import os
import shutil
import argparse

# dictionary with the installation information
NEW_FILE_NAME = "element_sedimentation.f90"
install_info = {
    # # star/make
    "makefile_base": {"mesa_path": "star/make", "replace": False},
    # # star/defaults
    "controls.defaults": {"mesa_path": "star/defaults", "replace": False},
    # # star/private
    "alloc.f90": {"mesa_path": "star/private", "replace": False},
    "ctrls_io.f90": {"mesa_path": "star/private", "replace": False},
    "evolve.f90": {"mesa_path": "star/private", "replace": False},
    "hydro_energy.f90": {"mesa_path": "star/private", "replace": False},
    "init.f90": {"mesa_path": "star/private", "replace": False},
    NEW_FILE_NAME: {"mesa_path": "star/private", "replace": True},
    # # star_data/private
    "star_controls.inc": {"mesa_path": "star_data/private", "replace": False},
    # # star_data/public
    "star_data_step_work.inc": {"mesa_path": "star_data/public", "replace": False},
}

verbose = False
debug = False


# modes of the installation file
DO_APPEND: int = 0  # append the new content after the anchor
DO_REPLACE: int = 1 # replace the anchor with the new content
DO_PREPEND: int = 2 # prepend the new content before the anchor
DO_APPEND_ALL: int = 3 # append the new content after all occurences of the anchor

# First, check if `$MESA_DIR` is set
if "MESA_DIR" not in os.environ:
    raise ValueError(
        "Please set the environment variable $MESA_DIR to the path of your MESA installation."
    )
else:
    mesa_dir = os.environ["MESA_DIR"]
    print(f"MESA_DIR is set to {mesa_dir}.") if debug else None

mesa_path = lambda *path: os.path.join(mesa_dir, *path)
# Define the path to the new code in ./code (relative to the script)
code_dir = "./code"
backup_dir = "./backup"

# Check if the code directory exists
if not os.path.exists(code_dir):
    raise FileNotFoundError(f"{code_dir} not found.")
else:
    print(f"Found {code_dir}.") if debug else None

# Check if the backup directory exists
if not os.path.exists(backup_dir):
    os.makedirs(backup_dir)
    print(f"Created {backup_dir}.") if verbose else None


# Define Class for the files that need to be modified
class InstallationFile:
    def __init__(
        self, name: str, replace_file: bool = False, *mesa_path_vars, **kwargs
    ) -> None:
        """Initializes the InstallationFile object.

        If replace_file is True, the file will be replaced with the new content.
        Otherwise, the new content will be inseted to the existing file.
        """
        self.name = name

        self.code_path = os.path.join(code_dir, name)

        # check if the code file exists
        if not os.path.exists(self.code_path):
            raise FileNotFoundError(f"{self.code_path} not found.")

        self.mesa_path: str = mesa_path(*mesa_path_vars)
        self.full_mesa_path: str = os.path.join(self.mesa_path, name)
        self.module_path: str = os.path.dirname(self.mesa_path)

        # check if the mesa file exists
        # skip the check for the new file
        if not os.path.exists(self.full_mesa_path) and self.name != NEW_FILE_NAME:
            raise FileNotFoundError(f"{self.full_mesa_path} not found.")

        self.backup_path: str = os.path.join(backup_dir, name + ".bak")
        self.replace_file: bool = replace_file

    def backup(self) -> None:
        """Backs up the mesa file to the backup directory."""
        # check if the backup directory exists
        if not os.path.exists(backup_dir):
            os.makedirs(backup_dir)

        # check if the backup file already exists and return if it does
        if os.path.exists(self.backup_path):
            return

        # skip backup for the new file
        if self.name == NEW_FILE_NAME:
            return

        shutil.copy(self.full_mesa_path, self.backup_path)
        print(f"Backed up {self.name} to {self.backup_path}.") if verbose else None

    def read_installation_file(self, **kwargs) -> tuple[list[str], list[list[str]], list[int]]:
        """Reads the installation file and returns the anchors and the new contents."""

        with open(self.code_path, "r") as code_file:
            # the anchors always follows the pattern ! anchor = '<anchor>'
            # first, we find all lines that start with `! anchor =`

            lines = code_file.readlines()

            # check if the file is empty
            if not lines:
                raise ValueError(f"{self.code_path} is empty.")

            # find all lines their respective indices that start with `! anchor =`
            i_anchors: list[int] = [
                i for i, line in enumerate(lines) if line.startswith("! anchor =")
            ]
            
            # check if the file contains the `replace` keyword
            # in this case, we replace the anchor instead of inserting new content after it
            insert_mode: list[int] = []
            for line in lines:
                if line.startswith("! anchor ="):
                    if line.strip().split(";")[-1].endswith("replace"):
                        insert_mode.append(DO_REPLACE)
                    elif line.strip().split(";")[-1].endswith("prepend"):
                        insert_mode.append(DO_PREPEND)
                    elif line.strip().split(";")[-1].endswith("all"):
                        insert_mode.append(DO_APPEND_ALL)
                    else:   # default is to append
                        insert_mode.append(DO_APPEND)

            # if there are no anchors, raise an error
            if not i_anchors:
                raise ValueError(f"No anchor found in {self.code_path}.")

            anchors: list[str] = [
                lines[i].strip().split("=", maxsplit=1)[1].split(";")[0].strip().strip("'\"")
                for i in i_anchors
            ]

            # the new contents are between the anchors
            new_contents: list[list[str]] = []
            for i, anchor in enumerate(anchors):
                # find the start and end index of the new content
                i_start = i_anchors[i] + 1
                i_end = i_anchors[i + 1] if i < len(i_anchors) - 1 else len(lines)

                # append the new content to the new_content list
                new_contents.append(lines[i_start:i_end])

            # if there is no '\n' at the end of new content, add it
            for new_content in new_contents:
                if not new_content[-1].endswith("\n"):
                    new_content[-1] += "\n"

            # print anchor and new content
            print("Anchors:", *anchors) if debug else None
            print("New contents:", *new_content) if debug else None
            print("insert_mode:", *insert_mode) if debug else None

        return anchors, new_contents, insert_mode

    def modify(self, new_content: list, anchor: str, insert_mode : int, **kwargs) -> None:
        """Inserts new content into the mesa file at the anchor."""
        self.backup()
        lines = []
        with open(self.full_mesa_path, "r") as file:
            lines = file.readlines()

        # find all occurences of the anchor
        anchors: list[int] = [i for i, line in enumerate(lines) if anchor in line]

        # if there are no anchors, raise an error
        if not anchors:
            raise ValueError(f"Anchor {anchor} not found in {self.name}.")

        elif len(anchors) > 1 and insert_mode != DO_APPEND_ALL:
            raise ValueError(f"Multiple anchors {anchor} found in {self.name}.")

        for i_anchor in reversed(anchors):
            if insert_mode == DO_APPEND or insert_mode == DO_APPEND_ALL:
                # new content is sandwiched between the anchor and the next line
                i_start: int = i_anchor + 1
                i_end: int = i_anchor + 1

            elif insert_mode == DO_REPLACE:
                # new content replaces the anchor
                i_start: int = i_anchor
                i_end: int = i_anchor + 1

            elif insert_mode == DO_PREPEND:
                # new content is prepended before the anchor
                i_start: int = i_anchor
                i_end: int = i_anchor

            # check if the new content is already in the file
            if new_content == lines[i_start : i_start + len(new_content)]:
                print(f"Content already in {self.name}.") if verbose else None
                return

            lines[i_start : i_end] = new_content

        with open(self.full_mesa_path, "w") as file:
            file.writelines(lines)

        print(f"Modified {self.name}.") if verbose else None

        format = kwargs.get("format", False)
        if format:
            self.format(**kwargs)

    def format(self, **kwargs) -> None:
        """Formats the MESA file."""

        indentation: int = kwargs.get("indentation", 3)
        # test if indentation is an integer
        if not isinstance(indentation, int):
            raise ValueError("Indentation must be an integer.")

        shift: int = kwargs.get("shift", 0)
        # test if shift is an integer
        if not isinstance(shift, int):
            raise ValueError("Shift must be an integer.")

        # try to format the file. If it fails, print a warning.
        try:
            # format the file and save it to a temporary file
            os.system(
                f"findent -I{shift} -i{indentation} -funix < {self.full_mesa_path} > {self.full_mesa_path}.tmp"
            )

            # move the temporary file to the original file
            os.system(f"mv {self.full_mesa_path}.tmp {self.full_mesa_path}")

            # print a success message
            print(f"Formatted {self.name}.") if verbose else None

        except:
            (
                print(f"Could not format {self.name}. Is `findent` installed?")
                if verbose
                else None
            )

    def insert(self, **kwargs) -> None:
        """Reads the new code and inserts it into the relevant MESA file."""

        # if self.replace is True, throw an error
        if self.replace_file:
            raise ValueError("Cannot insert code if replace is True.")

        # get the anchor and the new content
        anchors, new_contents, insert_modes = self.read_installation_file(**kwargs)

        # modify the mesa file
        for anchor, new_content, insert_mode in zip(anchors, new_contents, insert_modes):
            self.modify(new_content, anchor, insert_mode, **kwargs)

    def install(self, **kwargs) -> None:
        """Installs the new code in the MESA directory.

        Description:
        If replace is True, the file will be replaced with the new content.
        Otherwise, the new content will be inseted to the existing file.
        """

        # first, backup the original file
        self.backup()

        # Either replace the file or insert the new code into the file
        if self.replace_file:
            shutil.copy(self.code_path, self.full_mesa_path)

            # if the file is the new file, we are done
            if self.name == NEW_FILE_NAME:
                print(f"Copied {self.name} into {self.mesa_path}.") if verbose else None
            else:
                print(f"Replaced {self.name} with new file.") if verbose else None
        else:
            self.insert(**kwargs)
            (
                print(f"Finished inserting the new code into {self.name}.")
                if verbose
                else None
            )

    def uninstall(self) -> None:
        """Uninstalls the new code from the MESA directory."""

        # for the new file, simply remove the file
        if self.name == NEW_FILE_NAME:
            os.remove(self.full_mesa_path)
            print(f"Removed {self.name}.") if verbose else None
            return

        # check if the backup file exists
        if not os.path.exists(self.backup_path):
            raise FileNotFoundError(f"No backup found for {self.name}.")

        # restore the original file
        shutil.copy(self.backup_path, self.full_mesa_path)
        print(f"Restored {self.name}.") if verbose else None


# Init the argument parser
parser = argparse.ArgumentParser(description="Install or uninstall the software.")

# Add the arguments
parser.add_argument(
    "command",
    choices=[
        "install",
        "uninstall",
        "reinstall",
        "clean",
        "soft_install",
        "soft_uninstall",
        "soft_reinstall",
        "single_install",
        "single_uninstall",
        "single_reinstall",
    ],
    help="The command to execute.",
)

parser.add_argument('-f', '--filename', help='Filename for single operations')

# Parse the arguments
args = parser.parse_args()

print(f"Executing command: {args.command}") if debug else None

if args.command == "install":

    # install all the files in the install_info dictionary
    for key, value in install_info.items():
    
        print("-" * 80) if verbose else None
        file = InstallationFile(key, value["replace"], value["mesa_path"])
        file.install()

    # execute ./clean and ./install in the MESA directory
    os.system("cd $MESA_DIR; ./clean; ./install")


    # execute ./clean and ./mk in the $MESA_DIR/eos directory
    os.system("cd $MESA_DIR/eos; ./clean; ./mk; ./export")

    # execute ./clean and ./mk in the $MESA_DIR/star_data directory
    os.system("cd $MESA_DIR/star_data; ./clean; ./mk; ./export")

    # execute ./clean and ./mk in the $MESA_DIR/star directory
    os.system("cd $MESA_DIR/star; ./clean; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print("Gentle mixing was succesfully installed!")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "soft_install":
    # recompiling the code without running ./clean and ./install
    # install all the files in the install_info dictionary
    for key, value in install_info.items():
        print("-" * 80) if verbose else None
        file = InstallationFile(key, value["replace"], value["mesa_path"])
        file.install()

    # ask the user if they want to clean the mesa dirs
    clean: bool = (
        input("Do you want to clean the mesa directories? (y/n): ").lower() == "y"
    )
    if clean:

        # execute ./clean and ./mk in the $MESA_DIR/eos directory
        os.system("cd $MESA_DIR/eos; ./clean; ./mk; ./export")

        # execute ./clean and ./mk in the $MESA_DIR/star_data directory
        os.system("cd $MESA_DIR/star_data; ./clean; ./mk; ./export")

        # execute ./clean and ./mk in the $MESA_DIR/star directory
        os.system("cd $MESA_DIR/star; ./clean; ./mk; ./export")
    else:

        # execute ./mk in the $MESA_DIR/eos directory
        os.system("cd $MESA_DIR/eos; ./mk; ./export")

        # execute ./mk in the $MESA_DIR/star_data directory
        os.system("cd $MESA_DIR/star_data; ./mk; ./export")

        # execute ./mk in the $MESA_DIR/star directory
        os.system("cd $MESA_DIR/star; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print("Gentle mixing was succesfully compilled!")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "uninstall":

    # uninstall all the files in the install_info dictionary
    for key, value in install_info.items():
        file = InstallationFile(key, value["replace"], value["mesa_path"])
        file.uninstall()

    # execute ./clean and ./install in the MESA directory
    os.system("cd $MESA_DIR; ./clean; ./install")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print("Gentle mixing was succesfully uninstalled.")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "soft_uninstall":

    # uninstall all the files in the install_info dictionary
    for key, value in install_info.items():
        file = InstallationFile(key, value["replace"], value["mesa_path"])
        file.uninstall()

    # ask the user if they want to clean the mesa dirs
    clean: bool = (
        input("Do you want to clean the mesa directories? (y/n): ").lower() == "y"
    )
    if clean:

        # execute ./clean and ./mk in the $MESA_DIR/eos directory
        os.system("cd $MESA_DIR/eos; ./clean; ./mk; ./export")

        # execute ./clean and ./mk in the $MESA_DIR/star_data directory
        os.system("cd $MESA_DIR/star_data; ./clean; ./mk; ./export")

        # execute ./clean and ./mk in the $MESA_DIR/star directory
        os.system("cd $MESA_DIR/star; ./clean; ./mk; ./export")
    else:

        # execute ./mk in the $MESA_DIR/eos directory
        os.system("cd $MESA_DIR/eos; ./mk; ./export")

        # execute ./mk in the $MESA_DIR/star_data directory
        os.system("cd $MESA_DIR/star_data; ./mk; ./export")

        # execute ./mk in the $MESA_DIR/star directory
        os.system("cd $MESA_DIR/star; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print("Gentle mixing was succesfully recompiled to the original version.")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "soft_reinstall":

    # uninstall all the files in the install_info dictionary
    for key, value in install_info.items():
        file = InstallationFile(key, value["replace"], value["mesa_path"])
        file.uninstall()
    
    # install all the files in the install_info dictionary
    for key, value in install_info.items():
        file = InstallationFile(key, value["replace"], value["mesa_path"])
        file.install()

    # ask the user if they want to clean the mesa dirs
    clean: bool = (
        input("Do you want to clean the mesa directories? (y/n): ").lower() == "y"
    )

    if clean:

        # execute ./clean and ./mk in the $MESA_DIR/eos directory
        os.system("cd $MESA_DIR/eos; ./clean; ./mk; ./export")

        # execute ./clean and ./mk in the $MESA_DIR/star_data directory
        os.system("cd $MESA_DIR/star_data; ./clean; ./mk; ./export")

        # execute ./clean and ./mk in the $MESA_DIR/star directory
        os.system("cd $MESA_DIR/star; ./clean; ./mk; ./export")
    else:

        # execute ./mk in the $MESA_DIR/eos directory
        os.system("cd $MESA_DIR/eos; ./mk; ./export")

        # execute ./mk in the $MESA_DIR/star_data directory
        os.system("cd $MESA_DIR/star_data; ./mk; ./export")

        # execute ./mk in the $MESA_DIR/star directory
        os.system("cd $MESA_DIR/star; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print("Gentle mixing was succesfully recompiled to the original version.")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "single_install":
    # install a single file
    # get filename from -f argument or ask user
    if args.filename:
        file_name = args.filename
    else:
        file_name = input("Enter the file name: ")
    file = InstallationFile(
        file_name,
        install_info[file_name]["replace"],
        install_info[file_name]["mesa_path"],
    )
    file.install()

    # execute ./mk and ./export in the relevant directory
    os.system(f"cd {file.module_path}; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print(f"{file_name} was succesfully installed.")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "single_uninstall":
    # uninstall a single file
    # get filename from -f argument or ask user
    if args.filename:
        file_name = args.filename
    else:
        file_name = input("Enter the file name: ")
    file = InstallationFile(
        file_name,
        install_info[file_name]["replace"],
        install_info[file_name]["mesa_path"],
    )
    file.uninstall()

    # execute ./mk and ./export in the relevant directory
    os.system(f"cd {file.module_path}; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print(f"{file_name} was succesfully uninstalled.")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "single_reinstall":
    # ask the user for the file name
    # first, restore the original file
    if args.filename:
        file_name = args.filename
    else:
        file_name = input("Enter the file name: ")

    file = InstallationFile(
        file_name,
        install_info[file_name]["replace"],
        install_info[file_name]["mesa_path"],
    )
    file.uninstall()

    # install the new file
    file.install()

    # execute ./mk and ./export in the relevant directory
    os.system(f"cd {file.module_path}; ./mk; ./export")

    # remind the user to run ./clean and ./mk in their work directory
    print("-" * 80)
    print(f"{file_name} was succesfully reinstalled.")
    print("Don't forget to run ./clean and ./mk in your work directory.")
    print("-" * 80)

elif args.command == "clean":
    # remove the backup directory
    shutil.rmtree(backup_dir)
    print(f"Removed {backup_dir}.") if verbose else None

else:
    raise ValueError("Invalid command.")
