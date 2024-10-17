#!/usr/bin/env python3
"""
Compile COMPASS
===============

Usage: compile.py [options] [INSTALL_PATH]

Options:
    INSTALL_PATH          The install path. Default is $COMPASS_INSTALL_ROOT.
    --rt                  Allow parallel build on isolated cores. Default is False.
    --config              Force the configuration of the project. Default is False.

Environment variables:
    COMPASS_INSTALL_ROOT  The install path. Default is $project_dir/local.
    COMPASS_DEBUG         The build type. Default is Release.
    CUDA_SM               The CUDA architecture. Default is native.
    PYTHON_VERSION        The Python version. Default is the current version.

This script compiles and installs COMPASS. It can be used in two ways:
=
1. Without any arguments, it will use the default installation path. You can
   change it by setting the environment variable COMPASS_INSTALL_ROOT.

2. You can also provide an installation path as an argument.

The build type can be set with the environment variable COMPASS_DEBUG. The
default is Release. Set it to Debug to build a debug version.

The CUDA architecture can be set with the environment variable CUDA_SM. The
default is native. You can also set it to a specific architecture, like 75.

The Python version can be set with the environment variable PYTHON_VERSION.
The default is the current version. You can also set it to a specific version,
like 3.8. If you set it to OFF, the Python module will not be built.

The script generates a CMakePresets.json file and calls cmake to build the
project. It then installs the project using the install target.

Finally, it copies the binaries to the installation path and sets the
capabilities for the milk binary.

Example:
    $ ./compile.py /opt/compass
    $ COMPASS_INSTALL_ROOT=/opt/compass ./compile.py
    $ COMPASS_DEBUG=Debug ./compile.py /opt/compass
    $ CUDA_SM=75 ./compile.py /opt/compass
    $ PYTHON_VERSION=3.8 ./compile.py /opt/compass

Note:
    This script is intended to be used with the CI/CD pipeline. It is not
    intended to be used for a development environment. For that, you should
    use the standard cmake commands.
"""

from docopt import docopt
import json
import os
from pathlib import Path
from rich.console import Console
import subprocess
import sys


def get_install_path(args, project_dir):
    """Get the installation path."""
    # Get the installation path
    default_path = os.environ.get("COMPASS_INSTALL_ROOT", default=str(project_dir / "local"))

    if args:
        install_path = args
    else:
        install_path = default_path

    # Get the build type
    build_type = os.environ.get("COMPASS_DEBUG", default="Release")
    return install_path, build_type


def get_env_vars():
    """Get environment variables."""
    cuda_sm = os.environ.get("CUDA_SM", "native")
    build_libs = os.environ.get("BUILD_LIBS", default="on").lower() in [
        "yes",
        "true",
        "on",
    ]
    python_version_ = f"{sys.version_info.major}.{sys.version_info.minor}"
    python_version = os.environ.get("PYTHON_VERSION", default=python_version_)
    return cuda_sm, build_libs, python_version


def run_command(project_dir, cmd):
    """Print and run a command."""
    with Console() as console:
        try:
            sp_kwargs = {"shell": True, "check": True, "cwd": project_dir}
            console.print(f"Running {cmd}", style="bold green")
            subprocess.run(cmd, **sp_kwargs)
        except subprocess.CalledProcessError as e:
            console.print(
                f"Error calling {e.cmd} with return code {e.returncode}",
                style="bold red",
            )
            sys.exit(e.returncode)


def generate_CMakePresets_json(install_path, build_dir, build_type):
    """Generate the command args to install."""
    cuda_args, build_libs, python_version = get_env_vars()

    # CMakePresets.json file
    cmake_presets = {
        "version": 3,
        "configurePresets": [
            {
                "name": "default",
                "binaryDir": str(build_dir),
                "cacheVariables": {
                    "CMAKE_BUILD_TYPE": build_type,
                    "CMAKE_INSTALL_PREFIX": install_path,
                    "CMAKE_TOOLCHAIN_FILE":
                        "$env{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake",
                    "CMAKE_CUDA_ARCHITECTURES": cuda_args,
                    "CMAKE_POSITION_INDEPENDENT_CODE": "ON",
                    "build_python_module": "ON" if python_version != "OFF" else "OFF",
                    "libs_build": build_libs,
                },
            }
        ],
    }

    # Add CMAKE_EXPORT_COMPILE_COMMANDS for Debug builds
    if build_type == "Debug":
        preset = cmake_presets["configurePresets"][0]
        preset["cacheVariables"]["CMAKE_EXPORT_COMPILE_COMMANDS"] = "ON"

    # Write CMakePresets.json file
    with open("CMakePresets.json", "w") as f:
        json.dump(cmake_presets, f, indent=4)


if __name__ == "__main__":
    arguments = docopt(__doc__, version="Compile COMPASS 2.0", options_first=True)
    if arguments["--rt"]:
        import os

        pid = os.getpid()
        os.system(f"sudo chrt -r --pid 1 {pid}")
    args = arguments["INSTALL_PATH"]

    project_dir = Path(__file__).absolute().parent

    build_dir = project_dir / "build"
    install_path, build_type = get_install_path(args, project_dir)

    generate_CMakePresets_json(install_path, build_dir, build_type)

    # test if the Makefile exists in the build directory
    if not build_dir.exists() or not (build_dir / "Makefile").exists() or arguments["--config"]:
        run_command(project_dir, "cmake --preset=default")

    run_command(project_dir, f"cmake --build {build_dir} --target install --parallel")
