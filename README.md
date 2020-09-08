# Overview

HavoqGT (Highly Asynchronous Visitor Queue Graph Toolkit) is a framework for
expressing asynchronous vertex-centric graph algorithms.  It provides a visitor
interface, where actions are defined at an individual vertex level.
This code was developed at Lawrence Livermore National Laboratory.

Built in C++, the framework provides a runtime for parallel communication and
algorithm termination detection.   V0.1 is an initial release with only MPI support.
All graph data is stored in mmaped files, using Boost.Interprocess and Memory 
Mapped (mmap) I/O.   Large graphs that cannot fit in main-memory may still be
processed using mmap as external memory.  For best results, high speed Flash 
devices are preferred for external memory storage.

For documentation, see http://havoqgt.bitbucket.org

--------------------------------------------------------------------------------
# About

## Authors

* Roger A Pearce (rpearce at llnl dot gov)
* Keita Iwabuchi (kiwabuchi at llnl dot gov)
* Tahsin A Reza (reza2 at llnl dot gov)

## License

HavoqGT is distributed under the terms of the MIT license.
All new contributions must be made under this license.

See [LICENSE](LICENSE) and [NOTICE](NOTICE) for details.

SPDX-License-Identifier: MIT

## Release

LLNL-CODE-644630