#!/usr/bin/env python

import os
import sys
import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


def get_install_path(args, project_dir):
    """Get the installation path."""
    default_path = os.environ.get('COMPASS_INSTALL_PATH',
                                  default=project_dir / 'local')
    build_type = os.environ.get('COMPASS_DEBUG',
                                default='Release')
    return args[0] if args else default_path, build_type


def get_env_vars():
    """Get environment variables."""
    cuda_sm = (['-o', 'cuda_sm=' + os.environ.get('CUDA_SM')]
               if os.environ.get('CUDA_SM') else [])
    do_half = (os.environ.get('COMPASS_DO_HALF', default='off').lower()
               in ['yes', 'true', 'on'])
    build_libs = (os.environ.get('BUILD_LIBS', default='on').lower()
                  in ['yes', 'true', 'on'])
    python_version_ = f'{sys.version_info.major}.{sys.version_info.minor}'
    python_version = os.environ.get('PYTHON_VERSION', default=python_version_)
    return cuda_sm, do_half, build_libs, python_version


def generate_conan_cmd(build_type):
    """Generate the command args to install compass."""
    cuda_args, do_half, build_libs, python_version = get_env_vars()
    return ' '.join(['conan', 'build', '.', '-b', 'missing',
                     '-o', f'half={do_half}',
                     '-o', f'libs={build_libs}',
                     '-o', f'python_version={python_version}',
                     '-s', f'build_type={build_type}',
                     *cuda_args
                     ])


def run_command(project_dir, cmd):
    """Print and run a command."""
    try:
        sp_kwargs = {'shell': True, 'check': True, 'cwd': project_dir}
        console.print(f'Running {cmd}', style="bold green")
        subprocess.run(cmd, **sp_kwargs)
    except subprocess.CalledProcessError as e:
        console.print(f'Error calling {e.cmd} with return code {e.returncode}',
                      style="bold red")
        sys.exit(e.returncode)


def main():
    args = sys.argv[1:]
    project_dir = Path(__file__).absolute().parent

    install_path, build_type = get_install_path(args, project_dir)
    build_dir = project_dir / 'build' / build_type
    conan_preset = f"--preset conan-{build_type.lower()}"

    if not build_dir.exists():
        conan_cmd = generate_conan_cmd(build_type)
        run_command(project_dir, conan_cmd)
        cmake_arg = f'-DCMAKE_INSTALL_PREFIX={install_path} {conan_preset} .'
        run_command(project_dir, f'cmake {cmake_arg}')

    cmake_arg = f'--build {build_dir} -t install {conan_preset} --parallel'
    run_command(project_dir, f'cmake {cmake_arg}')


if __name__ == '__main__':
    main()
