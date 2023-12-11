#!/usr/bin/env python

import os
import sys
import subprocess
from pathlib import Path

def get_install_path(args, project_dir):
    """Get the installation path."""
    default_path = os.environ.get('COMPASS_INSTALL_PATH', project_dir / 'local')
    return args[0] if args else default_path

def get_env_vars():
    """Get environment variables."""
    cuda_sm = ['-o', 'cuda_sm=' + os.environ.get('CUDA_SM')] if os.environ.get('CUDA_SM') else []
    do_half = os.environ.get('COMPASS_DO_HALF', default='off').lower() in ['yes', 'true', 'on']
    build_libs = os.environ.get('BUILD_LIBS', default='on').lower() in ['yes', 'true', 'on']
    python_version = os.environ.get('PYTHON_VERSION', default=f'{sys.version_info.major}.{sys.version_info.minor}')
    build_type = os.environ.get('COMPASS_DEBUG', default='Release')
    return cuda_sm, do_half, build_libs, python_version, build_type

def generate_conan_cmd(do_half, build_libs, python_version, build_type, cuda_args):
    """Generate the command args to install compass."""
    return ' '.join(['conan', 'build', '.', '-b', 'missing',
                     '-o', f'half={do_half}',
                     '-o', f'libs={build_libs}',
                     '-o', f'python_version={python_version}',
                     '-s', f'build_type={build_type}',
                     *cuda_args
                     ])

def run_command(project_dir, cmd):
    """Print and run a command."""
    print(f'Running {cmd}')
    try:
        sp_kwargs = {'shell': True, 'check': True, 'cwd': project_dir}
        subprocess.run(cmd, **sp_kwargs)
    except subprocess.CalledProcessError as e:
        print(f'Error calling {e.cmd} with return code {e.returncode}')
        sys.exit(e.returncode)

def main():
    args = sys.argv[1:]
    project_dir = Path(__file__).absolute().parent

    install_path = get_install_path(args, project_dir)
    cuda_args, do_half, build_libs, python_version, build_type = get_env_vars()
    
    conan_cmd = generate_conan_cmd(do_half, build_libs, python_version, build_type, cuda_args)

    build_dir = project_dir / 'build' / build_type

    if not build_dir.exists():
        run_command(project_dir, conan_cmd)
        cmake_cmd = f'cmake -DCMAKE_INSTALL_PREFIX={install_path} --preset conan-{build_type.lower()} .'
        run_command(project_dir, cmake_cmd)

    cmake_cmd = f'cmake --build {build_dir} -t install --preset conan-{build_type.lower()} --parallel'
    run_command(project_dir, cmake_cmd)

if __name__ == '__main__':
    main()