# Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
# See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: MIT

import pathlib
import json
import subprocess
#from clippy.error import ClippyConfigurationError
from error import ClippyConfigurationError

CLIPPY_ENV = 'CLIPPY_ENV'
CLIPPY_CFG = '.clippy'
JSON_FLAG = '--clippy-help'


def _is_exe(f):
    '''
    Returns true if `f` is an executable file.
    '''
    if not f.exists():
        return False
    if not f.is_file():
        return False
    s = oct(f.stat().st_mode)[-3:]
    return '7' in s


def get_registered_commands(logger, cmd_dict=None):
    '''
    Returns a dictionary of namespaces with keys of str (representing the namespace)
    and values of dict, with keys of method names, and vals of dicts representing
    the arguments and docstrings.
    Error as appropriate.
    '''

    namespaces = {}
    for namespace, path in cmd_dict.items():
        p = pathlib.Path(path)
        if not p.is_dir():
            raise ClippyConfigurationError(f'executable directory {p} does not exist')

        namespaces[namespace] = {}
        for f in p.iterdir():
            logger.debug(f"trying {f}")
            if _is_exe(f):
                logger.debug(f'running {f} {JSON_FLAG}')
                cmd = [f, JSON_FLAG]
                exe = subprocess.run(cmd, capture_output=True)
                if exe.returncode:
                    logger.warn(f'{f} exited with error code {exe.returncode} - ignoring')
                else:
                    try:
                        j = json.loads(exe.stdout)
                        j['exe_name'] = str(f)
                        ns = namespaces[namespace]
                        ns[j['method_name']] = j
                        logger.debug(f'Adding {f} to valid commands under namespace {namespace}')
                    except json.JSONDecodeError:
                        logger.warn(f'JSON parsing error for {f} - ignoring')

    logger.debug(f"namespaces = {namespaces}")
    return namespaces
