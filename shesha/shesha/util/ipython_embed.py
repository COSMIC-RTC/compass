## @package   shesha.util
## @brief     Shesha utilities
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @date      2022/01/24
## @copyright 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>
#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>

# COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the 
# License, or any later version.

# COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
# If not, see <https://www.gnu.org/licenses/>

# Copyright (C) 2011-2024 COSMIC Team <https//://github.com/COSMIC-RTC/compass>
try:
    from IPython.terminal.prompts import Prompts, Token
    from IPython.terminal.embed import embed as std_embed

    class CustomPrompt(Prompts):
        name = ""

        def __init__(self, shell):
            Prompts.__init__(self, shell)

        def setName(self, name):
            self.name = name

        def in_prompt_tokens(self, cli=None):
            return [
                    (Token.Prompt, self.name),
                    (Token.Prompt, ' In ['),
                    (Token.PromptNum, str(self.shell.execution_count)),
                    (Token.Prompt, ']: '),
            ]

        def out_prompt_tokens(self):
            return [
                    (Token.OutPrompt, self.name),
                    (Token.OutPrompt, ' Out['),
                    (Token.OutPromptNum, str(self.shell.execution_count)),
                    (Token.OutPrompt, ']: '),
            ]

    def embed(name: str="", loc_vars: dict=None):
        from traitlets.config import Config

        glob_vars = globals()
        if loc_vars is None:
            glob_vars.update(locals())
        else:
            glob_vars.update(loc_vars)

        cfg = Config()
        cfg.InteractiveShellApp.gui = "qt5"
        cfg.TerminalInteractiveShell.prompts_class = CustomPrompt
        CustomPrompt.name = name
        std_embed(config=cfg,
                  banner1='Dropping into IPython, type %gui qt5 to unlock GUI')

except ImportError:
    import code

    def embed(name: str="", loc_vars: dict=None):
        import sys
        sys.ps1 = name + " >>> "
        sys.ps2 = name + " ... "

        glob_vars = globals()
        if loc_vars is None:
            glob_vars.update(locals())
        else:
            glob_vars.update(loc_vars)
        shell = code.InteractiveConsole(glob_vars)
        shell.interact()
