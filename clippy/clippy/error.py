# Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
# See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: MIT

# This file contains custom Clippy errors.

class ClippyError(Exception):
    '''
    This is a top-level custom exception for Clippy.
    '''
    pass


class ClippyConfigurationError(ClippyError):
    '''
    This error represents a configuration error on user input.
    '''
    pass


class ClippyBackendError(ClippyError):
    '''
    This error should be thrown when the backend returns an abend.
    '''
    pass


class ClippyValidationError(ClippyError):
    '''
    This error represents a validation error in the inputs to a clippy job.
    '''
    pass
