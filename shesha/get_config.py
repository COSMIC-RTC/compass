"""Check computer configuration
.
Usage:
  get_config.py [options]

Options:
  -h --help                 Show this screen.
  -o, --output=<outputfile> Filename to write [default: config.md].
"""

import shlex, sys
from subprocess import PIPE, Popen
from docopt import docopt
from glob import glob
from os.path import basename

import shesha


def launch_cmd(command_line: str):
    pipes = command_line.split('|')

    pipe_in = None
    ps = None
    while len(pipes) > 0:
        cmd_args = shlex.split(pipes[0])
        ps = Popen(cmd_args, stdout=PIPE, stdin=pipe_in, encoding="utf-8")
        pipe_in = ps.stdout
        del pipes[0]

    return ps.communicate()[0]


def check_command(file_handler, command_line: str, *, desc: str = None):
    """ Check the command output and format it for log

    Check the command output and format it for log

    Parameters:
        file_handler : (file) : file where the command output will be written

        command_line : (str) : command to launch

        desc : (str) : description of the command to launch (optionnal)
    """
    if desc is None:
        desc = command_line

    file_handler.write("\n## " + desc + "\n\n")
    file_handler.write("```bash\n" + launch_cmd(command_line) + "```\n")


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["--output"]
    if filename is None:
        filename = "config.md"

    with open(filename, 'w') as file:
        file.write("# COMPASS configuration script output\n")
        cmds = {
                "hostname": None,
                "lspci | grep -i vga": "VGA Cards",
                "cat /etc/lsb-release": "Distribution",
                "uname -a": "Kernel",
        }
        for cmd in cmds:
            check_command(file, cmd, desc=cmds[cmd])

        file.write("\n## python version\n")
        file.write("\n```bash\n" + sys.version + "\n ``` \n")

        file.write("\n## shesha version\n")
        file.write("\n```bash\n" + shesha.__version__ + "\n ``` \n")
        check_command(file, "conda list")

        if launch_cmd("conda list | grep compass | wc -l") == "0\n":
            # development version
            file.write("\n## development version\n")

            so_files = glob("local/*/*.so")
            cmds = {}
            for so_file in so_files:
                cmds["ldd -d " + so_file + " | grep -v Py"] = basename(so_file)
            for cmd in cmds:
                check_command(file, cmd, desc=cmds[cmd])
        else:
            # public version
            check_command(file, "conda list | grep compass", desc="public version")
